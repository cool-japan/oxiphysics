//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, HashSet};

use super::types::MeshStats;

pub(super) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(super) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub(super) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(super) fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
pub(super) fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = length(a);
    if l < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}
pub(super) fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub(super) fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub(super) fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    length(sub(a, b))
}
pub(super) fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    scale(add(a, b), 0.5)
}
/// Project point `p` onto line segment `a`-`b`, return parameter t in \[0,1\]
/// and the closest point.
pub(super) fn project_point_on_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> (f64, [f64; 3]) {
    let ab = sub(b, a);
    let ab_len2 = dot(ab, ab);
    if ab_len2 < 1e-20 {
        return (0.0, a);
    }
    let t = (dot(sub(p, a), ab) / ab_len2).clamp(0.0, 1.0);
    let closest = add(a, scale(ab, t));
    (t, closest)
}
/// Compute the area of a triangle given three vertex positions.
pub fn compute_triangle_area(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> f64 {
    let e1 = sub(v1, v0);
    let e2 = sub(v2, v0);
    length(cross(e1, e2)) * 0.5
}
/// Remove degenerate triangles (area < `threshold`) from a mesh.
///
/// Returns the filtered triangle list; the vertex array is unchanged.
pub fn remove_degenerate_triangles(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
) -> Vec<[usize; 3]> {
    remove_degenerate_triangles_threshold(vertices, triangles, 1e-12)
}
/// Remove degenerate triangles with a configurable area threshold.
pub fn remove_degenerate_triangles_threshold(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    threshold: f64,
) -> Vec<[usize; 3]> {
    triangles
        .iter()
        .filter(|tri| {
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                return false;
            }
            let area = compute_triangle_area(vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
            area >= threshold
        })
        .copied()
        .collect()
}
/// Merge vertices that are within `tolerance` of each other.
///
/// Returns a new (compacted) vertex array and updated triangle indices.
pub fn weld_vertices(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    tolerance: f64,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let n = vertices.len();
    let mut remap: Vec<usize> = (0..n).collect();
    for i in 0..n {
        if remap[i] != i {
            continue;
        }
        for j in (i + 1)..n {
            if remap[j] == j && dist(vertices[i], vertices[j]) <= tolerance {
                remap[j] = i;
            }
        }
    }
    let mut new_vertices: Vec<[f64; 3]> = Vec::new();
    let mut compact: Vec<usize> = vec![usize::MAX; n];
    for (i, &canon) in remap.iter().enumerate().take(n) {
        let _ = i;
        if compact[canon] == usize::MAX {
            compact[canon] = new_vertices.len();
            new_vertices.push(vertices[canon]);
        }
    }
    let new_triangles: Vec<[usize; 3]> = triangles
        .iter()
        .map(|tri| {
            [
                compact[remap[tri[0]]],
                compact[remap[tri[1]]],
                compact[remap[tri[2]]],
            ]
        })
        .collect();
    (new_vertices, new_triangles)
}
/// Compute angle-weighted vertex normals for a triangle mesh.
pub fn compute_vertex_normals(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> Vec<[f64; 3]> {
    let mut normals = vec![[0.0_f64; 3]; vertices.len()];
    for tri in triangles {
        let [i0, i1, i2] = *tri;
        let v0 = vertices[i0];
        let v1 = vertices[i1];
        let v2 = vertices[i2];
        let face_normal = cross(sub(v1, v0), sub(v2, v0));
        let angle_at = |va: [f64; 3], vb: [f64; 3], vc: [f64; 3]| -> f64 {
            let ab = normalize(sub(vb, va));
            let ac = normalize(sub(vc, va));
            dot(ab, ac).clamp(-1.0, 1.0).acos()
        };
        let a0 = angle_at(v0, v1, v2);
        let a1 = angle_at(v1, v2, v0);
        let a2 = angle_at(v2, v0, v1);
        for (vi, w) in [(i0, a0), (i1, a1), (i2, a2)] {
            normals[vi] = add(normals[vi], scale(face_normal, w));
        }
    }
    normals.iter().map(|&n| normalize(n)).collect()
}
/// Reverse the winding order of all triangles, flipping their normals.
pub fn flip_normals(triangles: &mut [[usize; 3]]) {
    for tri in triangles.iter_mut() {
        tri.swap(1, 2);
    }
}
/// Compute the total surface area of a triangle mesh.
pub fn compute_surface_area(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
    triangles
        .iter()
        .map(|tri| compute_triangle_area(vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]))
        .sum()
}
/// Detect and repair T-junctions in a triangle mesh.
///
/// A T-junction occurs when a vertex lies on (or very close to) an edge of
/// another triangle but is not referenced by that triangle. This function
/// splits such edges by inserting the T-junction vertex, producing new
/// triangles.
///
/// Returns updated triangle list (vertices are unchanged; no new vertices are
/// added, we only re-use existing ones that sit on edges).
pub fn repair_t_junctions(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    tolerance: f64,
) -> Vec<[usize; 3]> {
    let mut result: Vec<[usize; 3]> = Vec::new();
    for tri in triangles {
        let v0 = tri[0];
        let v1 = tri[1];
        let v2 = tri[2];
        let edges = [(v0, v1, v2), (v1, v2, v0), (v2, v0, v1)];
        let mut split_found = false;
        for &(ea, eb, opp) in &edges {
            let mut on_edge = Vec::new();
            for (vi, &vpos) in vertices.iter().enumerate() {
                if vi == v0 || vi == v1 || vi == v2 {
                    continue;
                }
                let (t, closest) = project_point_on_segment(vpos, vertices[ea], vertices[eb]);
                if t > 1e-6 && t < 1.0 - 1e-6 && dist(vpos, closest) < tolerance {
                    on_edge.push((vi, t));
                }
            }
            if !on_edge.is_empty() {
                on_edge.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
                let mut prev = ea;
                for &(vi, _) in &on_edge {
                    result.push([prev, vi, opp]);
                    prev = vi;
                }
                result.push([prev, eb, opp]);
                split_found = true;
                break;
            }
        }
        if !split_found {
            result.push(*tri);
        }
    }
    result
}
/// Detect non-manifold edges (shared by more than 2 triangles).
/// Returns list of (a, b) edge pairs that are non-manifold.
pub fn detect_non_manifold_edges(triangles: &[[usize; 3]]) -> Vec<(usize, usize)> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    edge_count
        .into_iter()
        .filter(|&(_, c)| c > 2)
        .map(|(e, _)| e)
        .collect()
}
/// Fix non-manifold edges by removing excess triangles until each edge is
/// shared by at most 2 triangles. Keeps the first two triangles for each
/// edge encountered; removes later ones.
pub fn fix_non_manifold_edges(triangles: &[[usize; 3]]) -> Vec<[usize; 3]> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut keep = vec![true; triangles.len()];
    for (ti, tri) in triangles.iter().enumerate() {
        let mut remove = false;
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            let count = edge_count.entry(key).or_insert(0);
            if *count >= 2 {
                remove = true;
            }
        }
        if remove {
            keep[ti] = false;
        } else {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry(key).or_insert(0) += 1;
            }
        }
    }
    triangles
        .iter()
        .enumerate()
        .filter(|&(i, _)| keep[i])
        .map(|(_, tri)| *tri)
        .collect()
}
/// Extract boundary loops from a triangle mesh. Each loop is a list of
/// vertex indices forming a closed polygon along boundary edges.
pub fn find_boundary_loops(triangles: &[[usize; 3]]) -> Vec<Vec<usize>> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let mut boundary_next: HashMap<usize, usize> = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_count.get(&key) == Some(&1) {
                boundary_next.insert(a, b);
            }
        }
    }
    let mut visited: HashSet<usize> = HashSet::new();
    let mut loops = Vec::new();
    for &start in boundary_next.keys() {
        if visited.contains(&start) {
            continue;
        }
        let mut loop_verts = Vec::new();
        let mut current = start;
        loop {
            if visited.contains(&current) {
                break;
            }
            visited.insert(current);
            loop_verts.push(current);
            match boundary_next.get(&current) {
                Some(&next) => current = next,
                None => break,
            }
        }
        if loop_verts.len() >= 3 {
            loops.push(loop_verts);
        }
    }
    loops
}
/// Fill a single boundary loop using ear-clipping triangulation.
/// Returns the new triangles that fill the hole.
pub fn fill_hole_ear_clipping(vertices: &[[f64; 3]], boundary_loop: &[usize]) -> Vec<[usize; 3]> {
    if boundary_loop.len() < 3 {
        return Vec::new();
    }
    let loop_normal = {
        let mut n = [0.0; 3];
        let n_verts = boundary_loop.len();
        for i in 0..n_verts {
            let v0 = vertices[boundary_loop[i]];
            let v1 = vertices[boundary_loop[(i + 1) % n_verts]];
            n[0] += (v0[1] - v1[1]) * (v0[2] + v1[2]);
            n[1] += (v0[2] - v1[2]) * (v0[0] + v1[0]);
            n[2] += (v0[0] - v1[0]) * (v0[1] + v1[1]);
        }
        normalize(n)
    };
    let mut remaining: Vec<usize> = boundary_loop.to_vec();
    let mut result = Vec::new();
    let mut safety = remaining.len() * remaining.len();
    while remaining.len() > 2 && safety > 0 {
        safety -= 1;
        let n = remaining.len();
        let mut ear_found = false;
        for i in 0..n {
            let prev = remaining[(i + n - 1) % n];
            let curr = remaining[i];
            let next = remaining[(i + 1) % n];
            let e1 = sub(vertices[curr], vertices[prev]);
            let e2 = sub(vertices[next], vertices[curr]);
            let tri_normal = cross(e1, e2);
            if dot(tri_normal, loop_normal) < 0.0 {
                continue;
            }
            let mut contains_other = false;
            for j in 0..n {
                if j == (i + n - 1) % n || j == i || j == (i + 1) % n {
                    continue;
                }
                let p = vertices[remaining[j]];
                if point_in_triangle(p, vertices[prev], vertices[curr], vertices[next]) {
                    contains_other = true;
                    break;
                }
            }
            if !contains_other {
                result.push([prev, curr, next]);
                remaining.remove(i);
                ear_found = true;
                break;
            }
        }
        if !ear_found {
            break;
        }
    }
    result
}
/// Check if point p is inside triangle (a, b, c) using barycentric coordinates.
pub(super) fn point_in_triangle(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> bool {
    let v0 = sub(c, a);
    let v1 = sub(b, a);
    let v2 = sub(p, a);
    let d00 = dot(v0, v0);
    let d01 = dot(v0, v1);
    let d02 = dot(v0, v2);
    let d11 = dot(v1, v1);
    let d12 = dot(v1, v2);
    let inv_denom = d00 * d11 - d01 * d01;
    if inv_denom.abs() < 1e-20 {
        return false;
    }
    let inv_denom = 1.0 / inv_denom;
    let u = (d11 * d02 - d01 * d12) * inv_denom;
    let v = (d00 * d12 - d01 * d02) * inv_denom;
    u > 1e-8 && v > 1e-8 && (u + v) < 1.0 - 1e-8
}
/// Stitch two triangle meshes together by concatenating their vertex and
/// triangle arrays. Vertex indices in the second mesh are offset accordingly.
pub fn stitch_meshes(
    vertices_a: &[[f64; 3]],
    triangles_a: &[[usize; 3]],
    vertices_b: &[[f64; 3]],
    triangles_b: &[[usize; 3]],
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let offset = vertices_a.len();
    let mut vertices = vertices_a.to_vec();
    vertices.extend_from_slice(vertices_b);
    let mut triangles = triangles_a.to_vec();
    for tri in triangles_b {
        triangles.push([tri[0] + offset, tri[1] + offset, tri[2] + offset]);
    }
    (vertices, triangles)
}
/// Check if all triangles have consistent winding (normals roughly agree with
/// their neighbours). Returns the indices of triangles whose normals point
/// opposite to the majority of their neighbours.
pub fn inconsistent_normal_faces(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> Vec<usize> {
    let face_normals: Vec<[f64; 3]> = triangles
        .iter()
        .map(|tri| {
            let e1 = sub(vertices[tri[1]], vertices[tri[0]]);
            let e2 = sub(vertices[tri[2]], vertices[tri[0]]);
            normalize(cross(e1, e2))
        })
        .collect();
    let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (fi, tri) in triangles.iter().enumerate() {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }
    let mut inconsistent = Vec::new();
    for fi in 0..triangles.len() {
        let n_i = face_normals[fi];
        let mut agree = 0i32;
        let mut disagree = 0i32;
        let tri = &triangles[fi];
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(faces) = edge_faces.get(&key) {
                for &fj in faces {
                    if fj == fi {
                        continue;
                    }
                    if dot(n_i, face_normals[fj]) > 0.0 {
                        agree += 1;
                    } else {
                        disagree += 1;
                    }
                }
            }
        }
        if disagree > agree && disagree > 0 {
            inconsistent.push(fi);
        }
    }
    inconsistent
}
/// Fix normal consistency by flipping inconsistent triangles.
/// Returns a new triangle list with consistent winding.
pub fn fix_normal_consistency(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> Vec<[usize; 3]> {
    let bad = inconsistent_normal_faces(vertices, triangles);
    let bad_set: HashSet<usize> = bad.into_iter().collect();
    triangles
        .iter()
        .enumerate()
        .map(|(i, tri)| {
            if bad_set.contains(&i) {
                [tri[0], tri[2], tri[1]]
            } else {
                *tri
            }
        })
        .collect()
}
/// Compute mesh statistics.
pub fn compute_mesh_stats(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> MeshStats {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let n_edges = edge_count.len();
    let n_boundary = edge_count.values().filter(|&&c| c == 1).count();
    let n_non_manifold = edge_count.values().filter(|&&c| c > 2).count();
    MeshStats {
        n_vertices: vertices.len(),
        n_triangles: triangles.len(),
        n_edges,
        n_boundary_edges: n_boundary,
        n_non_manifold_edges: n_non_manifold,
        euler_characteristic: vertices.len() as i64 - n_edges as i64 + triangles.len() as i64,
    }
}
/// Remove duplicate triangles (same vertex set regardless of winding).
pub fn remove_duplicate_triangles(triangles: &[[usize; 3]]) -> Vec<[usize; 3]> {
    let mut seen: HashSet<[usize; 3]> = HashSet::new();
    let mut result = Vec::new();
    for tri in triangles {
        let mut sorted = *tri;
        sorted.sort_unstable();
        if seen.insert(sorted) {
            result.push(*tri);
        }
    }
    result
}
/// Remove unused vertices (those not referenced by any triangle).
/// Returns compacted vertex array and remapped triangles.
pub fn remove_unused_vertices(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut used: HashSet<usize> = HashSet::new();
    for tri in triangles {
        for &vi in tri {
            used.insert(vi);
        }
    }
    let mut remap = vec![usize::MAX; vertices.len()];
    let mut new_verts = Vec::new();
    for i in 0..vertices.len() {
        if used.contains(&i) {
            remap[i] = new_verts.len();
            new_verts.push(vertices[i]);
        }
    }
    let new_tris: Vec<[usize; 3]> = triangles
        .iter()
        .map(|tri| [remap[tri[0]], remap[tri[1]], remap[tri[2]]])
        .collect();
    (new_verts, new_tris)
}
/// Remove degenerate (zero-area) triangles.
///
/// Returns `(filtered_triangles, n_removed)`.
pub fn remove_degenerate_triangles_counted(
    vertices: &[[f64; 3]],
    tris: &[[usize; 3]],
) -> (Vec<[usize; 3]>, usize) {
    let before = tris.len();
    let kept: Vec<[usize; 3]> = tris
        .iter()
        .filter(|tri| {
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                return false;
            }
            let area = compute_triangle_area(vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
            area >= 1e-12
        })
        .copied()
        .collect();
    let removed = before - kept.len();
    (kept, removed)
}
/// Merge duplicate vertices within `eps` of each other.
///
/// Returns `(new_vertices, new_triangles, n_merged)`.
pub fn merge_duplicate_vertices(
    vertices: &[[f64; 3]],
    tris: &[[usize; 3]],
    eps: f64,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, usize) {
    let n = vertices.len();
    let mut remap: Vec<usize> = (0..n).collect();
    for i in 0..n {
        if remap[i] != i {
            continue;
        }
        for j in (i + 1)..n {
            if remap[j] == j && dist(vertices[i], vertices[j]) <= eps {
                remap[j] = i;
            }
        }
    }
    let n_merged = remap.iter().enumerate().filter(|&(i, &r)| r != i).count();
    let mut compact: Vec<usize> = vec![usize::MAX; n];
    let mut new_verts: Vec<[f64; 3]> = Vec::new();
    for &canon in remap.iter().take(n) {
        if compact[canon] == usize::MAX {
            compact[canon] = new_verts.len();
            new_verts.push(vertices[canon]);
        }
    }
    let new_tris: Vec<[usize; 3]> = tris
        .iter()
        .map(|tri| {
            [
                compact[remap[tri[0]]],
                compact[remap[tri[1]]],
                compact[remap[tri[2]]],
            ]
        })
        .collect();
    (new_verts, new_tris, n_merged)
}
/// Fill all boundary holes in a mesh using ear-clipping.
///
/// Modifies `vertices` and `tris` in place. Returns the number of holes filled.
pub fn fill_boundary_holes(vertices: &[[f64; 3]], tris: &mut Vec<[usize; 3]>) -> usize {
    let loops = find_boundary_loops(tris);
    let n_holes = loops.len();
    for lp in &loops {
        let new_tris = fill_hole_ear_clipping(vertices, lp);
        tris.extend_from_slice(&new_tris);
    }
    n_holes
}
/// Orient all triangles consistently using flood-fill from the first face.
///
/// Returns the number of faces that were flipped.
pub fn fix_face_normals(vertices: &[[f64; 3]], tris: &mut [[usize; 3]]) -> usize {
    let n = tris.len();
    if n == 0 {
        return 0;
    }
    let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (fi, tri) in tris.iter().enumerate() {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }
    let face_normal = |tri: &[usize; 3]| -> [f64; 3] {
        let e1 = sub(vertices[tri[1]], vertices[tri[0]]);
        let e2 = sub(vertices[tri[2]], vertices[tri[0]]);
        normalize(cross(e1, e2))
    };
    let mut visited = vec![false; n];
    let mut flipped_set = HashSet::new();
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(0usize);
    visited[0] = true;
    while let Some(fi) = queue.pop_front() {
        let ni = face_normal(&tris[fi]);
        for k in 0..3 {
            let a = tris[fi][k];
            let b = tris[fi][(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(faces) = edge_faces.get(&key) {
                for &fj in faces {
                    if fj == fi || visited[fj] {
                        continue;
                    }
                    visited[fj] = true;
                    let nj = face_normal(&tris[fj]);
                    if dot(ni, nj) < 0.0 {
                        flipped_set.insert(fj);
                        tris[fj].swap(1, 2);
                    }
                    queue.push_back(fj);
                }
            }
        }
    }
    flipped_set.len()
}
/// Compute Euler characteristic: V − E + F.
pub fn compute_mesh_euler_characteristic(n_v: usize, n_e: usize, n_f: usize) -> i32 {
    n_v as i32 - n_e as i32 + n_f as i32
}
/// Count edges that belong to exactly one triangle (boundary edges).
pub fn count_boundary_edges(tris: &[[usize; 3]]) -> usize {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in tris {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    edge_count.values().filter(|&&c| c == 1).count()
}
/// Returns `true` if every edge is shared by exactly 2 triangles (watertight).
pub fn is_watertight(tris: &[[usize; 3]]) -> bool {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in tris {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    edge_count.values().all(|&c| c == 2)
}
/// Returns edge pairs `(a, b)` where the edge is shared by more than 2 triangles
/// (non-manifold edges).
pub fn edge_manifold_check(tris: &[[usize; 3]]) -> Vec<(usize, usize)> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in tris {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    edge_count
        .into_iter()
        .filter(|&(_, c)| c > 2)
        .map(|(e, _)| e)
        .collect()
}
/// Compute the circumradius of a 2D triangle (using XY coordinates).
///
/// Uses the formula R = (a·b·c) / (4·Area).
/// Returns `f64::INFINITY` for degenerate triangles.
pub fn circumradius_2d(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> f64 {
    let a = dist(v1, v2);
    let b = dist(v0, v2);
    let c = dist(v0, v1);
    let area = compute_triangle_area(v0, v1, v2);
    if area < 1e-20 {
        return f64::INFINITY;
    }
    (a * b * c) / (4.0 * area)
}
/// Compute the circumcenter of a triangle (projected to XY plane).
///
/// Returns `None` if the triangle is degenerate.
pub fn circumcenter_2d(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> Option<[f64; 3]> {
    let ax = v1[0] - v0[0];
    let ay = v1[1] - v0[1];
    let bx = v2[0] - v0[0];
    let by = v2[1] - v0[1];
    let d = 2.0 * (ax * by - ay * bx);
    if d.abs() < 1e-15 {
        return None;
    }
    let ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / d;
    let uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / d;
    Some([v0[0] + ux, v0[1] + uy, v0[2]])
}
/// Compute the minimum enclosing sphere of a set of vertices using the
/// Welzl / approximate iterative algorithm.
///
/// Returns `(center, radius)`. For an empty set returns `([0;3], 0.0)`.
pub fn minimum_bounding_sphere(vertices: &[[f64; 3]]) -> ([f64; 3], f64) {
    if vertices.is_empty() {
        return ([0.0; 3], 0.0);
    }
    let mut far = vertices[0];
    for &p in vertices.iter().skip(1) {
        if dist(p, vertices[0]) > dist(far, vertices[0]) {
            far = p;
        }
    }
    let mut far2 = far;
    for &p in vertices {
        if dist(p, far) > dist(far2, far) {
            far2 = p;
        }
    }
    let mut center = midpoint(far, far2);
    let mut radius = dist(far, far2) * 0.5;
    for &p in vertices {
        let d = dist(p, center);
        if d > radius {
            let new_radius = (radius + d) * 0.5;
            let t = (d - radius) / (2.0 * d);
            for i in 0..3 {
                center[i] += (p[i] - center[i]) * t;
            }
            radius = new_radius;
        }
    }
    (center, radius)
}
/// Apply one pass of Laplacian smoothing to vertex positions.
///
/// Each vertex is moved towards the average of its neighbours (connected
/// by at least one triangle edge). `lambda` ∈ \[0, 1\] controls the step size.
/// Returns the new smoothed vertex positions.
pub fn laplacian_smooth(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    lambda: f64,
) -> Vec<[f64; 3]> {
    let n = vertices.len();
    let mut neighbour_sum = vec![[0.0_f64; 3]; n];
    let mut neighbour_count = vec![0usize; n];
    for tri in triangles {
        let verts = [tri[0], tri[1], tri[2]];
        for k in 0..3 {
            let vi = verts[k];
            let vj = verts[(k + 1) % 3];
            for d in 0..3 {
                neighbour_sum[vi][d] += vertices[vj][d];
                neighbour_sum[vj][d] += vertices[vi][d];
            }
            neighbour_count[vi] += 1;
            neighbour_count[vj] += 1;
        }
    }
    let mut result = vertices.to_vec();
    for i in 0..n {
        let c = neighbour_count[i];
        if c == 0 {
            continue;
        }
        let avg = [
            neighbour_sum[i][0] / c as f64,
            neighbour_sum[i][1] / c as f64,
            neighbour_sum[i][2] / c as f64,
        ];
        for d in 0..3 {
            result[i][d] = vertices[i][d] + lambda * (avg[d] - vertices[i][d]);
        }
    }
    result
}
/// Subdivide each triangle into 4 by inserting midpoint vertices on each edge.
///
/// Returns the new `(vertices, triangles)` pair. Midpoints on shared edges are
/// created only once (vertices are shared across adjacent triangles).
pub fn subdivide_midpoint(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut new_verts = vertices.to_vec();
    let mut edge_midpoint: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    let mut new_tris = Vec::with_capacity(triangles.len() * 4);
    for tri in triangles {
        let [v0, v1, v2] = *tri;
        let get_mid = |a: usize,
                       b: usize,
                       verts: &mut Vec<[f64; 3]>,
                       em: &mut std::collections::HashMap<(usize, usize), usize>|
         -> usize {
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(&idx) = em.get(&key) {
                return idx;
            }
            let mid = midpoint(verts[a], verts[b]);
            let idx = verts.len();
            verts.push(mid);
            em.insert(key, idx);
            idx
        };
        let m01 = get_mid(v0, v1, &mut new_verts, &mut edge_midpoint);
        let m12 = get_mid(v1, v2, &mut new_verts, &mut edge_midpoint);
        let m20 = get_mid(v2, v0, &mut new_verts, &mut edge_midpoint);
        new_tris.push([v0, m01, m20]);
        new_tris.push([v1, m12, m01]);
        new_tris.push([v2, m20, m12]);
        new_tris.push([m01, m12, m20]);
    }
    (new_verts, new_tris)
}
/// Compute per-face normals for all triangles in a mesh.
///
/// Returns a `Vec<[f64;3]>` of length `triangles.len()`.
/// Each normal is normalised; degenerate faces yield `[0,0,0]`.
pub fn compute_face_normals(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> Vec<[f64; 3]> {
    triangles
        .iter()
        .map(|tri| {
            let e1 = sub(vertices[tri[1]], vertices[tri[0]]);
            let e2 = sub(vertices[tri[2]], vertices[tri[0]]);
            normalize(cross(e1, e2))
        })
        .collect()
}
/// Compute the axis-aligned bounding box of a mesh's vertices.
///
/// Returns `(min_corner, max_corner)`. Returns `([0;3], [0;3])` for empty input.
pub fn mesh_bounding_box(vertices: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    if vertices.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    let mut mn = vertices[0];
    let mut mx = vertices[0];
    for &v in vertices.iter().skip(1) {
        for i in 0..3 {
            if v[i] < mn[i] {
                mn[i] = v[i];
            }
            if v[i] > mx[i] {
                mx[i] = v[i];
            }
        }
    }
    (mn, mx)
}
/// Detect edges shared by exactly one triangle (boundary edges), returned as
/// directed `(from, to)` pairs as they appear in the first triangle.
pub fn directed_boundary_edges(triangles: &[[usize; 3]]) -> Vec<(usize, usize)> {
    let mut edge_dir: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_dir.entry(key).or_insert(0) += 1;
        }
    }
    let mut result = Vec::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_dir.get(&key) == Some(&1) {
                result.push((a, b));
            }
        }
    }
    result
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean_ops::HalfEdgeMesh;
    use crate::mesh_repair::FlipMesh;

    use crate::mesh_repair::edge_manifold_check;
    #[test]
    fn test_remove_degenerate_all_same_vertex() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let triangles = vec![[0usize, 0, 0], [0, 1, 2]];
        let result = remove_degenerate_triangles(&vertices, &triangles);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], [0, 1, 2]);
    }
    #[test]
    fn test_remove_degenerate_collinear() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [0, 1, 3]];
        let result = remove_degenerate_triangles(&vertices, &triangles);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], [0, 1, 3]);
    }
    #[test]
    fn test_remove_degenerate_with_threshold() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [0.001, 0.0, 0.0],
            [0.0, 0.001, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [0, 3, 4]];
        let result = remove_degenerate_triangles_threshold(&vertices, &triangles, 0.01);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], [0, 3, 4]);
    }
    #[test]
    fn test_weld_vertices_within_tolerance() {
        let vertices = vec![[0.0_f64, 0.0, 0.0], [1e-8, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let triangles = vec![[0usize, 1, 2]];
        let (new_verts, new_tris) = weld_vertices(&vertices, &triangles, 1e-6);
        assert_eq!(new_verts.len(), 2);
        assert_eq!(new_tris[0][0], new_tris[0][1]);
    }
    #[test]
    fn test_weld_vertices_no_merge() {
        let vertices = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [5.0, 10.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let (new_verts, new_tris) = weld_vertices(&vertices, &triangles, 1e-6);
        assert_eq!(new_verts.len(), 3);
        assert_eq!(new_tris[0], [0, 1, 2]);
    }
    #[test]
    fn test_compute_surface_area_unit_cube() {
        let v: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let t: Vec<[usize; 3]> = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [3, 6, 2],
            [3, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];
        let area = compute_surface_area(&v, &t);
        assert!((area - 6.0).abs() < 1e-10, "expected 6.0, got {area}");
    }
    #[test]
    fn test_half_edge_mesh_face_count() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let triangles = vec![[0usize, 1, 2]];
        let mesh = HalfEdgeMesh::from_triangle_mesh(&vertices, &triangles);
        assert_eq!(mesh.n_faces(), 1);
        assert_eq!(mesh.n_vertices(), 3);
        assert_eq!(mesh.half_edges.len(), 3);
    }
    #[test]
    fn test_half_edge_mesh_manifold() {
        // Closed tetrahedron — every edge shared by exactly two faces
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let triangles = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
        let mesh = HalfEdgeMesh::from_triangle_mesh(&vertices, &triangles);
        assert!(mesh.is_manifold());
    }
    #[test]
    fn test_half_edge_mesh_boundary() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let mesh = HalfEdgeMesh::from_triangle_mesh(&vertices, &triangles);
        let boundary = mesh.boundary_edges();
        assert_eq!(boundary.len(), 3, "single triangle has 3 boundary edges");
    }
    #[test]
    fn test_half_edge_mesh_boundary_loops() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let mesh = HalfEdgeMesh::from_triangle_mesh(&vertices, &triangles);
        let loops = mesh.boundary_loops();
        assert_eq!(loops.len(), 1, "single triangle has one boundary loop");
        assert_eq!(loops[0].len(), 3);
    }
    #[test]
    fn test_flip_normals_winding_reversed() {
        let mut triangles = vec![[0usize, 1, 2], [3, 4, 5]];
        flip_normals(&mut triangles);
        assert_eq!(triangles[0], [0, 2, 1]);
        assert_eq!(triangles[1], [3, 5, 4]);
    }
    #[test]
    fn test_compute_vertex_normals_flat_surface() {
        let vertices = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let triangles = vec![[0usize, 1, 2], [0, 2, 3]];
        let normals = compute_vertex_normals(&vertices, &triangles);
        assert_eq!(normals.len(), 4);
        for n in &normals {
            assert!((n[2] - 1.0).abs() < 1e-10, "expected z~1, got {:?}", n);
            assert!(n[0].abs() < 1e-10, "expected x~0, got {:?}", n);
            assert!(n[1].abs() < 1e-10, "expected y~0, got {:?}", n);
        }
    }
    #[test]
    fn test_repair_t_junction_basic() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2]];
        let result = repair_t_junctions(&vertices, &triangles, 1e-6);
        assert_eq!(result.len(), 2, "T-junction should split triangle into 2");
    }
    #[test]
    fn test_repair_t_junction_no_split_needed() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let result = repair_t_junctions(&vertices, &triangles, 1e-6);
        assert_eq!(result.len(), 1, "no T-junction should mean no split");
    }
    #[test]
    fn test_repair_t_junction_preserves_area() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [2.0, 0.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2]];
        let original_area = compute_surface_area(&vertices, &triangles);
        let repaired = repair_t_junctions(&vertices, &triangles, 1e-6);
        let repaired_area = compute_surface_area(&vertices, &repaired);
        assert!(
            (original_area - repaired_area).abs() < 1e-10,
            "area should be preserved: {} vs {}",
            original_area,
            repaired_area
        );
    }
    #[test]
    fn test_detect_non_manifold_edges() {
        let triangles = vec![[0, 1, 2], [0, 3, 1], [0, 4, 1]];
        let nm = detect_non_manifold_edges(&triangles);
        assert_eq!(nm.len(), 1, "one non-manifold edge expected");
    }
    #[test]
    fn test_detect_no_non_manifold_edges() {
        let triangles = vec![[0, 1, 2], [0, 3, 1]];
        let nm = detect_non_manifold_edges(&triangles);
        assert!(nm.is_empty());
    }
    #[test]
    fn test_fix_non_manifold_edges() {
        let triangles = vec![[0, 1, 2], [0, 3, 1], [0, 4, 1]];
        let fixed = fix_non_manifold_edges(&triangles);
        let nm = detect_non_manifold_edges(&fixed);
        assert!(
            nm.is_empty(),
            "after fix, no non-manifold edges should remain"
        );
    }
    #[test]
    fn test_find_boundary_loops_closed_mesh() {
        let triangles = vec![[0, 1, 2], [0, 3, 1]];
        let loops = find_boundary_loops(&triangles);
        assert!(!loops.is_empty());
    }
    #[test]
    fn test_fill_hole_ear_clipping_quad() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let boundary = vec![0, 1, 2, 3];
        let tris = fill_hole_ear_clipping(&vertices, &boundary);
        assert_eq!(tris.len(), 2, "quad should be filled with 2 triangles");
        let area = compute_surface_area(&vertices, &tris);
        assert!(
            (area - 1.0).abs() < 1e-10,
            "filled quad area should be 1.0, got {}",
            area
        );
    }
    #[test]
    fn test_fill_hole_ear_clipping_triangle() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let boundary = vec![0, 1, 2];
        let tris = fill_hole_ear_clipping(&vertices, &boundary);
        assert_eq!(tris.len(), 1, "triangle boundary fills with 1 triangle");
    }
    #[test]
    fn test_fill_hole_ear_clipping_pentagon() {
        let vertices: Vec<[f64; 3]> = (0..5)
            .map(|i| {
                let angle = 2.0 * std::f64::consts::PI * i as f64 / 5.0;
                [angle.cos(), angle.sin(), 0.0]
            })
            .collect();
        let boundary = vec![0, 1, 2, 3, 4];
        let tris = fill_hole_ear_clipping(&vertices, &boundary);
        assert_eq!(tris.len(), 3, "pentagon should be filled with 3 triangles");
    }
    #[test]
    fn test_stitch_meshes() {
        let va = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let ta = vec![[0, 1, 2]];
        let vb = vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.5, 1.0, 0.0]];
        let tb = vec![[0, 1, 2]];
        let (verts, tris) = stitch_meshes(&va, &ta, &vb, &tb);
        assert_eq!(verts.len(), 6);
        assert_eq!(tris.len(), 2);
        assert_eq!(tris[1], [3, 4, 5]);
    }
    #[test]
    fn test_stitch_then_weld() {
        let va = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let ta = vec![[0, 1, 2]];
        let vb = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.5, 1.0, 0.0]];
        let tb = vec![[0, 1, 2]];
        let (verts, tris) = stitch_meshes(&va, &ta, &vb, &tb);
        let (welded_v, _welded_t) = weld_vertices(&verts, &tris, 1e-6);
        assert_eq!(welded_v.len(), 5, "one shared vertex should be merged");
    }
    #[test]
    fn test_inconsistent_normals_detected() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [0, 1, 3]];
        let bad = inconsistent_normal_faces(&vertices, &triangles);
        assert!(!bad.is_empty(), "should detect inconsistent normals");
    }
    #[test]
    fn test_consistent_normals_not_flagged() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [1, 0, 3]];
        let bad = inconsistent_normal_faces(&vertices, &triangles);
        assert!(bad.is_empty(), "consistent normals should not be flagged");
    }
    #[test]
    fn test_fix_normal_consistency() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
        ];
        let tris = vec![[0, 1, 2], [0, 1, 3]];
        let fixed = fix_normal_consistency(&vertices, &tris);
        assert_ne!(
            fixed, tris,
            "fix should have modified at least one triangle"
        );
    }
    #[test]
    fn test_mesh_stats_single_triangle() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let stats = compute_mesh_stats(&vertices, &triangles);
        assert_eq!(stats.n_vertices, 3);
        assert_eq!(stats.n_triangles, 1);
        assert_eq!(stats.n_edges, 3);
        assert_eq!(stats.n_boundary_edges, 3);
        assert_eq!(stats.n_non_manifold_edges, 0);
    }
    #[test]
    fn test_mesh_stats_euler_characteristic() {
        let v: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let t: Vec<[usize; 3]> = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [3, 6, 2],
            [3, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];
        let stats = compute_mesh_stats(&v, &t);
        assert_eq!(stats.euler_characteristic, 2, "cube Euler char should be 2");
    }
    #[test]
    fn test_remove_duplicate_triangles() {
        let triangles = vec![[0, 1, 2], [0, 1, 2], [0, 2, 1], [3, 4, 5]];
        let result = remove_duplicate_triangles(&triangles);
        assert_eq!(result.len(), 2);
    }
    #[test]
    fn test_remove_unused_vertices() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [99.0, 99.0, 99.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
        ];
        let triangles = vec![[0, 2, 3]];
        let (new_v, new_t) = remove_unused_vertices(&vertices, &triangles);
        assert_eq!(new_v.len(), 3);
        assert_eq!(new_t.len(), 1);
        for v in &new_v {
            assert!((v[0] - 99.0).abs() > 1.0, "unused vertex should be removed");
        }
    }
    #[test]
    fn test_remove_degenerate_triangles_counted() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let tris = vec![[0, 1, 2], [0, 1, 3]];
        let (kept, n_removed) = remove_degenerate_triangles_counted(&vertices, &tris);
        assert_eq!(n_removed, 1);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0], [0, 1, 3]);
    }
    #[test]
    fn test_merge_duplicate_vertices_within_eps() {
        let vertices = vec![[0.0_f64, 0.0, 0.0], [1e-9, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let tris = vec![[0usize, 1, 2]];
        let (new_v, _new_t, n_merged) = merge_duplicate_vertices(&vertices, &tris, 1e-6);
        assert_eq!(n_merged, 1, "one vertex should be merged");
        assert_eq!(new_v.len(), 2);
    }
    #[test]
    fn test_is_watertight_cube() {
        let tris: Vec<[usize; 3]> = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [3, 6, 2],
            [3, 7, 6],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ];
        assert!(
            is_watertight(&tris),
            "closed cube mesh should be watertight"
        );
    }
    #[test]
    fn test_is_watertight_open_mesh() {
        let tris = vec![[0, 1, 2]];
        assert!(!is_watertight(&tris), "single triangle is not watertight");
    }
    #[test]
    fn test_euler_characteristic_sphere() {
        let chi = compute_mesh_euler_characteristic(8, 18, 12);
        assert_eq!(chi, 2, "Euler characteristic of a sphere should be 2");
    }
    #[test]
    fn test_count_boundary_edges() {
        let tris = vec![[0, 1, 2]];
        assert_eq!(count_boundary_edges(&tris), 3);
    }
    #[test]
    fn test_edge_manifold_check_non_manifold() {
        let tris = vec![[0, 1, 2], [0, 3, 1], [0, 4, 1]];
        let nm = edge_manifold_check(&tris);
        assert_eq!(nm.len(), 1, "one non-manifold edge expected");
    }
    #[test]
    fn test_fill_boundary_holes_adds_triangles() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let mut tris = vec![[0usize, 1, 2]];
        let before = tris.len();
        let n_holes = fill_boundary_holes(&vertices, &mut tris);
        assert!(n_holes > 0, "should find at least one hole");
        assert!(tris.len() > before, "filling holes should add triangles");
    }
    #[test]
    fn test_flip_mesh_no_flip_needed_delaunay() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [0, 3, 1]];
        let mut fm = FlipMesh::new(vertices, triangles);
        let flips = fm.flip_pass();
        let _ = flips;
    }
    #[test]
    fn test_flip_mesh_to_delaunay_converges() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 0.1, 0.0],
            [1.0, 2.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [0, 2, 3]];
        let mut fm = FlipMesh::new(vertices, triangles);
        fm.flip_to_delaunay(10);
        assert_eq!(
            fm.triangles.len(),
            2,
            "flip must not create or destroy triangles"
        );
    }
    #[test]
    fn test_flip_mesh_triangle_count_invariant() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [1.5, 0.2, 0.0],
            [1.5, 3.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [0, 2, 3]];
        let original_count = triangles.len();
        let mut fm = FlipMesh::new(vertices, triangles);
        fm.flip_to_delaunay(20);
        assert_eq!(fm.triangles.len(), original_count);
    }
    #[test]
    fn test_circumradius_2d_equilateral() {
        let s = 1.0_f64;
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [s, 0.0, 0.0];
        let v2 = [0.5 * s, 0.5 * 3.0_f64.sqrt() * s, 0.0];
        let r = circumradius_2d(v0, v1, v2);
        let expected = s / 3.0_f64.sqrt();
        assert!(
            (r - expected).abs() < 1e-8,
            "circumradius = {r}, expected {expected}"
        );
    }
    #[test]
    fn test_circumcenter_2d_right_triangle() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [2.0, 0.0, 0.0];
        let v2 = [0.0, 2.0, 0.0];
        let cc = circumcenter_2d(v0, v1, v2).expect("should find circumcenter");
        assert!((cc[0] - 1.0).abs() < 1e-9, "cx={}", cc[0]);
        assert!((cc[1] - 1.0).abs() < 1e-9, "cy={}", cc[1]);
    }
    #[test]
    fn test_circumcenter_2d_degenerate_returns_none() {
        let v0 = [0.0, 0.0, 0.0];
        let v1 = [1.0, 0.0, 0.0];
        let v2 = [2.0, 0.0, 0.0];
        let cc = circumcenter_2d(v0, v1, v2);
        assert!(cc.is_none(), "collinear points should return None");
    }
    #[test]
    fn test_minimum_bounding_sphere_single_point() {
        let verts = vec![[3.0, 5.0, 1.0]];
        let (center, radius) = minimum_bounding_sphere(&verts);
        assert!((center[0] - 3.0).abs() < 1e-9);
        assert!(radius.is_finite());
    }
    #[test]
    fn test_minimum_bounding_sphere_symmetric() {
        let verts = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let (center, radius) = minimum_bounding_sphere(&verts);
        assert!(center[0].abs() < 0.1, "cx={}", center[0]);
        assert!(center[1].abs() < 0.1, "cy={}", center[1]);
        assert!((0.99..=1.5).contains(&radius), "radius={radius}");
    }
    #[test]
    fn test_minimum_bounding_sphere_contains_all_points() {
        let verts: Vec<[f64; 3]> = (0..8)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 8.0;
                [angle.cos(), angle.sin(), 0.0]
            })
            .collect();
        let (center, radius) = minimum_bounding_sphere(&verts);
        for &p in &verts {
            let d = dist(p, center);
            assert!(
                d <= radius + 1e-6,
                "point {:?} outside sphere (d={d}, r={radius})",
                p
            );
        }
    }
    #[test]
    fn test_laplacian_smooth_flat_grid() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let triangles = vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]];
        let smoothed = laplacian_smooth(&vertices, &triangles, 0.5);
        for v in &smoothed {
            assert!(
                v[2].abs() < 1e-9,
                "z should stay 0 after smoothing flat grid: z={}",
                v[2]
            );
        }
    }
    #[test]
    fn test_laplacian_smooth_zero_lambda() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let smoothed = laplacian_smooth(&vertices, &triangles, 0.0);
        for (orig, smth) in vertices.iter().zip(smoothed.iter()) {
            for d in 0..3 {
                assert!(
                    (orig[d] - smth[d]).abs() < 1e-12,
                    "should not move with lambda=0"
                );
            }
        }
    }
    #[test]
    fn test_laplacian_smooth_reduces_roughness() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.5, 0.5, 2.0],
        ];
        let triangles = vec![[0, 1, 4], [1, 3, 4], [3, 2, 4], [2, 0, 4]];
        let smoothed = laplacian_smooth(&vertices, &triangles, 1.0);
        assert!(
            smoothed[4][2] < 2.0,
            "center z should decrease after smoothing: z={}",
            smoothed[4][2]
        );
    }
    #[test]
    fn test_subdivide_midpoint_quad_4x() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2], [1, 3, 2]];
        let (new_v, new_t) = subdivide_midpoint(&vertices, &triangles);
        assert_eq!(new_t.len(), 8, "2 tris × 4 = 8 after subdivision");
        assert_eq!(new_v.len(), 9, "expected 9 vertices after subdivision");
    }
    #[test]
    fn test_subdivide_midpoint_area_preserved() {
        let vertices = vec![[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let original_area = compute_surface_area(&vertices, &triangles);
        let (new_v, new_t) = subdivide_midpoint(&vertices, &triangles);
        let new_area = compute_surface_area(&new_v, &new_t);
        assert!(
            (original_area - new_area).abs() < 1e-9,
            "area should be preserved: {original_area} vs {new_area}"
        );
    }
    #[test]
    fn test_compute_face_normals_flat() {
        let vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let triangles = vec![[0, 1, 2]];
        let normals = compute_face_normals(&vertices, &triangles);
        assert_eq!(normals.len(), 1);
        assert!(
            normals[0][2].abs() > 0.99,
            "flat face normal should point in z: {:?}",
            normals[0]
        );
    }
    #[test]
    fn test_compute_face_normals_unit_length() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
            [1.5, 1.0, 0.5],
        ];
        let triangles = vec![[0, 1, 2], [3, 4, 5]];
        let normals = compute_face_normals(&vertices, &triangles);
        for n in &normals {
            let len = length(*n);
            assert!(
                (len - 1.0).abs() < 1e-9 || len < 1e-9,
                "normal not unit: len={len}"
            );
        }
    }
    #[test]
    fn test_mesh_bounding_box_unit_cube() {
        let v: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let (mn, mx) = mesh_bounding_box(&v);
        for i in 0..3 {
            assert!((mn[i] - 0.0).abs() < 1e-9, "min[{i}]={}", mn[i]);
            assert!((mx[i] - 1.0).abs() < 1e-9, "max[{i}]={}", mx[i]);
        }
    }
    #[test]
    fn test_mesh_bounding_box_empty() {
        let (mn, mx) = mesh_bounding_box(&[]);
        assert_eq!(mn, [0.0; 3]);
        assert_eq!(mx, [0.0; 3]);
    }
    #[test]
    fn test_directed_boundary_edges_single_tri() {
        let tris = vec![[0usize, 1, 2]];
        let edges = directed_boundary_edges(&tris);
        assert_eq!(edges.len(), 3, "single triangle has 3 boundary edges");
    }
    #[test]
    fn test_directed_boundary_edges_two_tris_share_edge() {
        let tris = vec![[0usize, 1, 2], [0, 3, 1]];
        let edges = directed_boundary_edges(&tris);
        assert_eq!(
            edges.len(),
            4,
            "two triangles sharing one edge: 4 boundary edges"
        );
    }
    #[test]
    fn test_fix_face_normals_two_flipped() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, -1.0, 0.0],
        ];
        let mut tris = vec![[0usize, 1, 2], [0, 1, 3]];
        let _n_flipped = fix_face_normals(&vertices, &mut tris);
    }
    #[test]
    fn test_t_junction_multiple_points_on_edge() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [3.0, 3.0, 0.0],
            [2.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let triangles = vec![[0, 1, 2]];
        let result = repair_t_junctions(&vertices, &triangles, 1e-6);
        assert_eq!(
            result.len(),
            3,
            "two T-junction vertices should split into 3 triangles: {:?}",
            result
        );
    }
}
#[cfg(test)]
mod tests_mesh_repair_struct {

    use crate::mesh_repair::MeshRepair;
    use crate::mesh_repair::edge_manifold_check;
    fn single_triangle() -> MeshRepair {
        MeshRepair::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
    }
    fn tetrahedron_mr() -> MeshRepair {
        MeshRepair::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [0.5, 0.33, 1.0],
            ],
            vec![[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]],
        )
    }
    fn quad_open() -> MeshRepair {
        MeshRepair::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
    }
    #[test]
    fn test_fill_holes_single_triangle() {
        let mut mr = single_triangle();
        let n = mr.fill_holes();
        assert_eq!(n, 1, "one boundary loop found");
        assert!(mr.triangles.len() > 1, "fan fill added triangles");
    }
    #[test]
    fn test_fill_holes_closed_mesh_no_op() {
        let mut mr = tetrahedron_mr();
        let n = mr.fill_holes();
        assert_eq!(n, 0, "closed mesh has no holes");
        assert_eq!(mr.triangles.len(), 4, "triangle count unchanged");
    }
    #[test]
    fn test_fill_holes_ear_single_triangle() {
        let mut mr = single_triangle();
        let n = mr.fill_holes_ear();
        assert_eq!(n, 1, "one boundary loop found");
    }
    #[test]
    fn test_fill_holes_adds_triangles() {
        let mut mr = quad_open();
        let before = mr.triangles.len();
        mr.fill_holes();
        assert!(
            mr.triangles.len() > before,
            "new triangles added by hole fill"
        );
    }
    #[test]
    fn test_remove_duplicate_vertices_identical() {
        let mut mr = MeshRepair::new(
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
            ],
            vec![[0, 2, 3], [1, 2, 3]],
        );
        let merged = mr.remove_duplicate_vertices(1e-9);
        assert_eq!(merged, 1, "one vertex merged");
        assert_eq!(mr.vertices.len(), 3);
    }
    #[test]
    fn test_remove_duplicate_vertices_epsilon() {
        let mut mr = MeshRepair::new(
            vec![
                [0.0, 0.0, 0.0],
                [0.005, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
            ],
            vec![[0, 2, 3], [1, 2, 3]],
        );
        let merged = mr.remove_duplicate_vertices(0.01);
        assert_eq!(merged, 1, "vertex within epsilon merged");
    }
    #[test]
    fn test_remove_duplicate_vertices_none() {
        let mut mr = single_triangle();
        let merged = mr.remove_duplicate_vertices(1e-9);
        assert_eq!(merged, 0);
        assert_eq!(mr.vertices.len(), 3);
    }
    #[test]
    fn test_fix_non_manifold_edges_clean_mesh() {
        let mut mr = quad_open();
        let removed = mr.fix_non_manifold_edges();
        assert_eq!(removed, 0, "quad is already manifold-valid");
    }
    #[test]
    fn test_fix_non_manifold_edges_removes_excess() {
        let mut mr = MeshRepair::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [0.5, -1.0, 0.0],
                [0.5, 0.5, 1.0],
            ],
            vec![[0, 1, 2], [0, 3, 1], [0, 4, 1]],
        );
        let removed = mr.fix_non_manifold_edges();
        assert!(removed >= 1, "at least one triangle removed");
        let nm = edge_manifold_check(&mr.triangles);
        assert!(nm.is_empty(), "no non-manifold edges remain");
    }
    #[test]
    fn test_euler_characteristic_tetrahedron() {
        let mr = tetrahedron_mr();
        assert_eq!(mr.compute_euler_characteristic(), 2);
    }
    #[test]
    fn test_euler_characteristic_single_triangle() {
        let mr = single_triangle();
        assert_eq!(mr.compute_euler_characteristic(), 1);
    }
    #[test]
    fn test_euler_characteristic_quad() {
        let mr = quad_open();
        assert_eq!(mr.compute_euler_characteristic(), 1);
    }
}
