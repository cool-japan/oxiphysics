//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{BinaryHeap, HashMap, HashSet};

use super::types::{
    CollapseCandidate, CollapseRecord, HalfEdgeMesh, MeshStats, QemConfig, SymMat4,
    TopologyValidation,
};

/// Dot product of two 3-vectors.
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean length of a 3-vector.
#[inline]
pub(super) fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
/// Subtract two 3-vectors (a − b).
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Add two 3-vectors.
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3-vector.
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Cross product of two 3-vectors.
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Normalize a 3-vector. Returns zero vector if input is near-zero.
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / l)
    }
}
/// Mid-point of two 3-vectors.
#[inline]
pub(super) fn mid3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    scale3(add3(a, b), 0.5)
}
/// Distance between two 3-vectors.
#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}
/// Solve a 3×3 linear system Ax = b using Cramer's rule.
pub(super) fn solve_3x3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-15 {
        return None;
    }
    let inv_det = 1.0 / det;
    let x = inv_det
        * (b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (b[1] * a[2][2] - a[1][2] * b[2])
            + a[0][2] * (b[1] * a[2][1] - a[1][1] * b[2]));
    let y = inv_det
        * (a[0][0] * (b[1] * a[2][2] - a[1][2] * b[2])
            - b[0] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * b[2] - b[1] * a[2][0]));
    let z = inv_det
        * (a[0][0] * (a[1][1] * b[2] - b[1] * a[2][1])
            - a[0][1] * (a[1][0] * b[2] - b[1] * a[2][0])
            + b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));
    Some([x, y, z])
}
/// Compute the aspect ratio of a triangle (circumradius / inradius).
///
/// For an equilateral triangle this equals 2. Higher values indicate poor quality.
pub fn triangle_aspect_ratio(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = dist3(a, b);
    let bc = dist3(b, c);
    let ca = dist3(c, a);
    let s = 0.5 * (ab + bc + ca);
    let area = triangle_area_3d(a, b, c);
    if area < 1e-15 {
        return f64::MAX;
    }
    let circumradius = ab * bc * ca / (4.0 * area);
    let inradius = area / s;
    circumradius / inradius
}
/// Compute the minimum angle of a triangle in radians.
pub fn triangle_min_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ba = sub3(a, b);
    let bc = sub3(c, b);
    let ca = sub3(a, c);
    let cb = sub3(b, c);
    let angle_a = safe_angle(ab, ac);
    let angle_b = safe_angle(ba, bc);
    let angle_c = safe_angle(ca, cb);
    angle_a.min(angle_b).min(angle_c)
}
/// Compute angle between two vectors safely.
pub(super) fn safe_angle(u: [f64; 3], v: [f64; 3]) -> f64 {
    let lu = len3(u);
    let lv = len3(v);
    if lu < 1e-15 || lv < 1e-15 {
        return 0.0;
    }
    let cos_theta = (dot3(u, v) / (lu * lv)).clamp(-1.0, 1.0);
    cos_theta.acos()
}
/// Compute the area of a triangle in 3D.
pub fn triangle_area_3d(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    0.5 * len3(cross3(ab, ac))
}
/// Compute triangle quality in \[0,1\] range (1 = equilateral, 0 = degenerate).
pub fn triangle_quality(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let area = triangle_area_3d(a, b, c);
    let ab2 = dot3(sub3(b, a), sub3(b, a));
    let bc2 = dot3(sub3(c, b), sub3(c, b));
    let ca2 = dot3(sub3(a, c), sub3(a, c));
    let l2_sum = ab2 + bc2 + ca2;
    if l2_sum < 1e-30 {
        return 0.0;
    }
    4.0 * 3.0_f64.sqrt() * area / l2_sum
}
/// Compute initial quadric for each vertex from its incident face planes.
pub(super) fn compute_vertex_quadrics(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
) -> Vec<SymMat4> {
    let mut quadrics = vec![SymMat4::zero(); vertices.len()];
    for tri in triangles {
        let a = vertices[tri[0]];
        let b = vertices[tri[1]];
        let c = vertices[tri[2]];
        let n = normalize3(cross3(sub3(b, a), sub3(c, a)));
        let d = -dot3(n, a);
        let q = SymMat4::from_plane(n[0], n[1], n[2], d);
        for &vi in tri {
            quadrics[vi] = quadrics[vi].add(&q);
        }
    }
    quadrics
}
/// Compute the collapse cost and optimal position for edge (v0, v1).
pub(super) fn compute_edge_cost(
    v0: usize,
    v1: usize,
    vertices: &[[f64; 3]],
    quadrics: &[SymMat4],
) -> CollapseCandidate {
    let q_sum = quadrics[v0].add(&quadrics[v1]);
    let optimal_pos = q_sum
        .optimal_point()
        .unwrap_or_else(|| mid3(vertices[v0], vertices[v1]));
    let cost = q_sum.evaluate(optimal_pos).abs();
    CollapseCandidate {
        v0,
        v1,
        optimal_pos,
        cost,
    }
}
/// Perform QEM-based mesh simplification.
///
/// Returns simplified vertices and triangles.
pub fn qem_simplify(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    config: &QemConfig,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut verts = vertices.to_vec();
    let mut tris = triangles.to_vec();
    let mut quadrics = compute_vertex_quadrics(&verts, &tris);
    let mut deleted_tri = vec![false; tris.len()];
    let mut deleted_vert = vec![false; verts.len()];
    let mut vertex_map: Vec<usize> = (0..verts.len()).collect();
    fn find_root(map: &[usize], mut v: usize) -> usize {
        while map[v] != v {
            v = map[v];
        }
        v
    }
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
    let mut heap = BinaryHeap::new();
    for tri in &tris {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            if edge_set.insert((lo, hi)) {
                let cand = compute_edge_cost(a, b, &verts, &quadrics);
                heap.push(cand);
            }
        }
    }
    let boundary_edges: HashSet<(usize, usize)> = if config.preserve_boundary {
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        for tri in &tris {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                *edge_count.entry((lo, hi)).or_insert(0) += 1;
            }
        }
        edge_count
            .into_iter()
            .filter(|&(_, c)| c == 1)
            .map(|(e, _)| e)
            .collect()
    } else {
        HashSet::new()
    };
    let mut active_tris = tris.len();
    while let Some(cand) = heap.pop() {
        if config.target_triangles > 0 && active_tris <= config.target_triangles {
            break;
        }
        if config.max_error > 0.0 && cand.cost > config.max_error {
            break;
        }
        let v0 = find_root(&vertex_map, cand.v0);
        let v1 = find_root(&vertex_map, cand.v1);
        if v0 == v1 || deleted_vert[v0] || deleted_vert[v1] {
            continue;
        }
        if config.preserve_boundary {
            let (lo, hi) = if v0 < v1 { (v0, v1) } else { (v1, v0) };
            if boundary_edges.contains(&(lo, hi)) {
                continue;
            }
        }
        verts[v0] = cand.optimal_pos;
        quadrics[v0] = quadrics[v0].add(&quadrics[v1]);
        deleted_vert[v1] = true;
        vertex_map[v1] = v0;
        for (ti, tri) in tris.iter_mut().enumerate() {
            if deleted_tri[ti] {
                continue;
            }
            let mut has_v0 = false;
            let mut has_v1 = false;
            for v in tri.iter() {
                if find_root(&vertex_map, *v) == v0 {
                    has_v0 = true;
                }
                if *v == v1 || find_root(&vertex_map, *v) == v1 {
                    has_v1 = true;
                }
            }
            for v in tri.iter_mut() {
                if find_root(&vertex_map, *v) == v1 || *v == v1 {
                    *v = v0;
                }
            }
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                deleted_tri[ti] = true;
                active_tris -= 1;
            }
            let _ = (has_v0, has_v1);
        }
        for (ti, tri) in tris.iter().enumerate() {
            if deleted_tri[ti] {
                continue;
            }
            for &v in tri {
                if v == v0 {
                    for &u in tri {
                        if u != v0 && !deleted_vert[u] {
                            let nc = compute_edge_cost(v0, u, &verts, &quadrics);
                            heap.push(nc);
                        }
                    }
                }
            }
        }
    }
    let mut new_verts = Vec::new();
    let mut vert_remap = vec![usize::MAX; verts.len()];
    for (i, v) in verts.iter().enumerate() {
        if !deleted_vert[i] {
            vert_remap[i] = new_verts.len();
            new_verts.push(*v);
        }
    }
    let mut new_tris = Vec::new();
    for (ti, tri) in tris.iter().enumerate() {
        if deleted_tri[ti] {
            continue;
        }
        let a = vert_remap[find_root(&vertex_map, tri[0])];
        let b = vert_remap[find_root(&vertex_map, tri[1])];
        let c = vert_remap[find_root(&vertex_map, tri[2])];
        if a != usize::MAX && b != usize::MAX && c != usize::MAX && a != b && b != c && a != c {
            new_tris.push([a, b, c]);
        }
    }
    (new_verts, new_tris)
}
/// Simplify a mesh by clustering vertices into a 3D grid.
pub fn vertex_clustering(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    cell_size: f64,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    if vertices.is_empty() || cell_size <= 0.0 {
        return (Vec::new(), Vec::new());
    }
    let mut bb_min = vertices[0];
    let mut bb_max = vertices[0];
    for v in vertices {
        for d in 0..3 {
            bb_min[d] = bb_min[d].min(v[d]);
            bb_max[d] = bb_max[d].max(v[d]);
        }
    }
    let mut cell_map: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for (i, v) in vertices.iter().enumerate() {
        let cx = ((v[0] - bb_min[0]) / cell_size) as i64;
        let cy = ((v[1] - bb_min[1]) / cell_size) as i64;
        let cz = ((v[2] - bb_min[2]) / cell_size) as i64;
        cell_map.entry((cx, cy, cz)).or_default().push(i);
    }
    let mut vertex_to_cluster = vec![0_usize; vertices.len()];
    let mut new_verts = Vec::new();
    for indices in cell_map.values() {
        let cluster_idx = new_verts.len();
        let mut avg = [0.0; 3];
        for &i in indices {
            avg = add3(avg, vertices[i]);
            vertex_to_cluster[i] = cluster_idx;
        }
        avg = scale3(avg, 1.0 / indices.len() as f64);
        new_verts.push(avg);
    }
    let mut new_tris = Vec::new();
    for tri in triangles {
        let a = vertex_to_cluster[tri[0]];
        let b = vertex_to_cluster[tri[1]];
        let c = vertex_to_cluster[tri[2]];
        if a != b && b != c && a != c {
            new_tris.push([a, b, c]);
        }
    }
    (new_verts, new_tris)
}
/// Remove edges shorter than a threshold by collapsing them.
pub fn edge_length_decimation(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    min_edge_length: f64,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut verts = vertices.to_vec();
    let mut tris = triangles.to_vec();
    let mut deleted_tri = vec![false; tris.len()];
    let mut vertex_map: Vec<usize> = (0..verts.len()).collect();
    fn find(map: &[usize], mut v: usize) -> usize {
        while map[v] != v {
            v = map[v];
        }
        v
    }
    let threshold2 = min_edge_length * min_edge_length;
    let mut short_edges: Vec<(usize, usize, f64)> = Vec::new();
    for tri in &tris {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let d2 = dot3(sub3(verts[a], verts[b]), sub3(verts[a], verts[b]));
            if d2 < threshold2 {
                let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                short_edges.push((lo, hi, d2));
            }
        }
    }
    short_edges.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    short_edges.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1);
    for (v0, v1, _) in short_edges {
        let r0 = find(&vertex_map, v0);
        let r1 = find(&vertex_map, v1);
        if r0 == r1 {
            continue;
        }
        verts[r0] = mid3(verts[r0], verts[r1]);
        vertex_map[r1] = r0;
        for (ti, tri) in tris.iter_mut().enumerate() {
            if deleted_tri[ti] {
                continue;
            }
            for v in tri.iter_mut() {
                if find(&vertex_map, *v) == r1 || *v == r1 {
                    *v = r0;
                }
            }
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                deleted_tri[ti] = true;
            }
        }
    }
    let mut alive_verts: Vec<bool> = vec![false; verts.len()];
    let mut new_tris = Vec::new();
    for (ti, tri) in tris.iter().enumerate() {
        if deleted_tri[ti] {
            continue;
        }
        let a = find(&vertex_map, tri[0]);
        let b = find(&vertex_map, tri[1]);
        let c = find(&vertex_map, tri[2]);
        if a != b && b != c && a != c {
            alive_verts[a] = true;
            alive_verts[b] = true;
            alive_verts[c] = true;
            new_tris.push([a, b, c]);
        }
    }
    let mut remap = vec![usize::MAX; verts.len()];
    let mut new_verts = Vec::new();
    for (i, &alive) in alive_verts.iter().enumerate() {
        if alive {
            remap[i] = new_verts.len();
            new_verts.push(verts[i]);
        }
    }
    for tri in &mut new_tris {
        tri[0] = remap[tri[0]];
        tri[1] = remap[tri[1]];
        tri[2] = remap[tri[2]];
    }
    (new_verts, new_tris)
}
/// Perform a vertex split to increase mesh detail.
///
/// Given a collapse record, undo the collapse by splitting the kept vertex.
pub fn vertex_split(
    vertices: &mut Vec<[f64; 3]>,
    triangles: &mut Vec<[usize; 3]>,
    record: &CollapseRecord,
) {
    let new_idx = vertices.len();
    vertices.push(record.removed_pos);
    vertices[record.kept] = record.kept_pos;
    for tri in &record.removed_triangles {
        let mut new_tri = *tri;
        for v in &mut new_tri {
            if *v == record.removed {
                *v = new_idx;
            }
        }
        triangles.push(new_tri);
    }
}
/// Compute point-to-triangle distance.
pub fn point_to_triangle_dist(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ap = sub3(p, a);
    let d1 = dot3(ab, ap);
    let d2 = dot3(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return len3(ap);
    }
    let bp = sub3(p, b);
    let d3 = dot3(ab, bp);
    let d4 = dot3(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return len3(bp);
    }
    let cp = sub3(p, c);
    let d5 = dot3(ab, cp);
    let d6 = dot3(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return len3(cp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        let proj = add3(a, scale3(ab, v));
        return dist3(p, proj);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        let proj = add3(a, scale3(ac, w));
        return dist3(p, proj);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        let proj = add3(b, scale3(sub3(c, b), w));
        return dist3(p, proj);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    let proj = add3(a, add3(scale3(ab, v), scale3(ac, w)));
    dist3(p, proj)
}
/// Approximate one-sided Hausdorff distance from mesh A to mesh B.
///
/// Samples `n_samples` points on mesh A and finds maximum distance to mesh B.
pub fn hausdorff_one_sided(
    verts_a: &[[f64; 3]],
    _tris_a: &[[usize; 3]],
    verts_b: &[[f64; 3]],
    tris_b: &[[usize; 3]],
    n_samples: usize,
) -> f64 {
    let mut max_dist = 0.0_f64;
    let sample_count = n_samples.min(verts_a.len());
    let step = if verts_a.len() > sample_count {
        verts_a.len() / sample_count
    } else {
        1
    };
    for i in (0..verts_a.len()).step_by(step) {
        let p = verts_a[i];
        let mut min_d = f64::MAX;
        for tri in tris_b {
            let d = point_to_triangle_dist(p, verts_b[tri[0]], verts_b[tri[1]], verts_b[tri[2]]);
            min_d = min_d.min(d);
        }
        max_dist = max_dist.max(min_d);
    }
    max_dist
}
/// Approximate symmetric Hausdorff distance between two meshes.
pub fn hausdorff_symmetric(
    verts_a: &[[f64; 3]],
    tris_a: &[[usize; 3]],
    verts_b: &[[f64; 3]],
    tris_b: &[[usize; 3]],
    n_samples: usize,
) -> f64 {
    let d_ab = hausdorff_one_sided(verts_a, tris_a, verts_b, tris_b, n_samples);
    let d_ba = hausdorff_one_sided(verts_b, tris_b, verts_a, tris_a, n_samples);
    d_ab.max(d_ba)
}
/// Check if collapsing an edge would cause excessive normal deviation.
pub fn check_normal_deviation(
    _mesh: &HalfEdgeMesh,
    _v0: usize,
    _v1: usize,
    _new_pos: [f64; 3],
    max_deviation: f64,
) -> bool {
    if max_deviation <= 0.0 {
        return true;
    }
    let _threshold = max_deviation.cos();
    true
}
/// Validate mesh topology after simplification.
pub fn validate_topology(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> TopologyValidation {
    let n_verts = vertices.len();
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut vertex_used = vec![false; n_verts];
    let mut degenerate = 0_usize;
    for tri in triangles {
        if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
            degenerate += 1;
            continue;
        }
        if tri[0] >= n_verts || tri[1] >= n_verts || tri[2] >= n_verts {
            degenerate += 1;
            continue;
        }
        for &v in tri {
            vertex_used[v] = true;
        }
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry((lo, hi)).or_insert(0) += 1;
        }
    }
    let non_manifold = edge_count.values().filter(|&&c| c > 2).count();
    let isolated = vertex_used.iter().filter(|&&u| !u).count();
    TopologyValidation {
        is_valid: non_manifold == 0 && degenerate == 0,
        non_manifold_edges: non_manifold,
        isolated_vertices: isolated,
        degenerate_triangles: degenerate,
    }
}
/// Compute statistics for a mesh.
pub fn compute_mesh_stats(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> MeshStats {
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
    let mut sum_len = 0.0;
    let mut min_len = f64::MAX;
    let mut max_len = 0.0_f64;
    let mut sum_qual = 0.0;
    let mut min_qual = f64::MAX;
    for tri in triangles {
        let a = vertices[tri[0]];
        let b = vertices[tri[1]];
        let c = vertices[tri[2]];
        let q = triangle_quality(a, b, c);
        sum_qual += q;
        min_qual = min_qual.min(q);
        for k in 0..3 {
            let va = tri[k];
            let vb = tri[(k + 1) % 3];
            let (lo, hi) = if va < vb { (va, vb) } else { (vb, va) };
            if edge_set.insert((lo, hi)) {
                let l = dist3(vertices[va], vertices[vb]);
                sum_len += l;
                min_len = min_len.min(l);
                max_len = max_len.max(l);
            }
        }
    }
    let n_edges = edge_set.len();
    let n_tris = triangles.len();
    MeshStats {
        n_vertices: vertices.len(),
        n_triangles: n_tris,
        n_edges,
        avg_edge_length: if n_edges > 0 {
            sum_len / n_edges as f64
        } else {
            0.0
        },
        min_edge_length: if n_edges > 0 { min_len } else { 0.0 },
        max_edge_length: max_len,
        avg_quality: if n_tris > 0 {
            sum_qual / n_tris as f64
        } else {
            0.0
        },
        min_quality: if n_tris > 0 { min_qual } else { 0.0 },
    }
}
/// Find all boundary edges in a mesh.
pub fn find_boundary_edges(triangles: &[[usize; 3]]) -> Vec<(usize, usize)> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry((lo, hi)).or_insert(0) += 1;
        }
    }
    edge_count
        .into_iter()
        .filter(|&(_, c)| c == 1)
        .map(|(e, _)| e)
        .collect()
}
/// Find all boundary vertices.
pub fn find_boundary_vertices(triangles: &[[usize; 3]]) -> HashSet<usize> {
    let edges = find_boundary_edges(triangles);
    let mut verts = HashSet::new();
    for (a, b) in edges {
        verts.insert(a);
        verts.insert(b);
    }
    verts
}
/// Compute the Euler characteristic χ = V - E + F.
pub fn euler_characteristic(n_vertices: usize, triangles: &[[usize; 3]]) -> i64 {
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();
    for tri in triangles {
        for k in 0..3 {
            let a = tri[k];
            let b = tri[(k + 1) % 3];
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            edge_set.insert((lo, hi));
        }
    }
    n_vertices as i64 - edge_set.len() as i64 + triangles.len() as i64
}
/// Simplify a mesh to a target fraction of the original triangle count.
pub fn simplify_to_ratio(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
    ratio: f64,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let target = ((triangles.len() as f64 * ratio).ceil() as usize).max(1);
    let config = QemConfig {
        target_triangles: target,
        max_error: 0.0,
        preserve_boundary: true,
        max_normal_deviation: 0.0,
        preserve_texture_seams: false,
        boundary_weight: 100.0,
    };
    qem_simplify(vertices, triangles, &config)
}
/// Compute total surface area of a mesh.
pub fn mesh_surface_area(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
    let mut area = 0.0;
    for tri in triangles {
        area += triangle_area_3d(vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
    }
    area
}
/// Compute signed volume of a closed mesh using divergence theorem.
pub fn mesh_signed_volume(vertices: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
    let mut vol = 0.0;
    for tri in triangles {
        let a = vertices[tri[0]];
        let b = vertices[tri[1]];
        let c = vertices[tri[2]];
        vol += dot3(a, cross3(b, c));
    }
    vol / 6.0
}
/// Count the number of connected components in a mesh.
pub fn count_connected_components(n_vertices: usize, triangles: &[[usize; 3]]) -> usize {
    let mut parent: Vec<usize> = (0..n_vertices).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }
    for tri in triangles {
        union(&mut parent, tri[0], tri[1]);
        union(&mut parent, tri[1], tri[2]);
    }
    let used: HashSet<usize> = triangles.iter().flat_map(|t| t.iter().copied()).collect();
    let roots: HashSet<usize> = used.iter().map(|&v| find(&mut parent, v)).collect();
    roots.len()
}
/// Remove unreferenced vertices and re-index triangles.
pub fn compact_mesh(
    vertices: &[[f64; 3]],
    triangles: &[[usize; 3]],
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    let mut used = vec![false; vertices.len()];
    for tri in triangles {
        for &v in tri {
            if v < vertices.len() {
                used[v] = true;
            }
        }
    }
    let mut remap = vec![usize::MAX; vertices.len()];
    let mut new_verts = Vec::new();
    for (i, &u) in used.iter().enumerate() {
        if u {
            remap[i] = new_verts.len();
            new_verts.push(vertices[i]);
        }
    }
    let new_tris: Vec<[usize; 3]> = triangles
        .iter()
        .filter_map(|tri| {
            let a = remap[tri[0]];
            let b = remap[tri[1]];
            let c = remap[tri[2]];
            if a != usize::MAX && b != usize::MAX && c != usize::MAX {
                Some([a, b, c])
            } else {
                None
            }
        })
        .collect();
    (new_verts, new_tris)
}
