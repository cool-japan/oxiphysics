//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, HashSet};

use super::types::{DecimationMetrics, EdgeFeature, QuadricMatrix, SimpleMesh};

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
pub(super) fn length(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}
pub(super) fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = length(v);
    if len < 1e-300 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}
/// Compute per-vertex quadric matrices by summing Kp over adjacent triangles.
pub fn compute_vertex_quadrics(mesh: &SimpleMesh) -> Vec<QuadricMatrix> {
    let n = mesh.vertex_count();
    let mut quadrics: Vec<QuadricMatrix> = (0..n).map(|_| QuadricMatrix::zero()).collect();
    for t_idx in 0..mesh.triangle_count() {
        let [a, b, c] = mesh.triangles[t_idx];
        let normal = mesh.triangle_normal(t_idx);
        let va = mesh.vertices[a];
        let d = -(normal[0] * va[0] + normal[1] * va[1] + normal[2] * va[2]);
        let kp = QuadricMatrix::from_plane(normal[0], normal[1], normal[2], d);
        quadrics[a] = quadrics[a].add(&kp);
        quadrics[b] = quadrics[b].add(&kp);
        quadrics[c] = quadrics[c].add(&kp);
    }
    quadrics
}
/// Compute the optimal collapse target and error for collapsing edge (v1→v2).
///
/// Returns `(error, optimal_position)`. The optimal position minimises the
/// combined quadric error; if the combined Q is not invertible the midpoint is
/// used instead.
pub fn edge_collapse_cost(
    v1: [f64; 3],
    v2: [f64; 3],
    q1: &QuadricMatrix,
    q2: &QuadricMatrix,
) -> (f64, [f64; 3]) {
    let combined = q1.add(q2);
    let mid = [
        (v1[0] + v2[0]) * 0.5,
        (v1[1] + v2[1]) * 0.5,
        (v1[2] + v2[2]) * 0.5,
    ];
    let a00 = combined.q[0][0];
    let a01 = combined.q[0][1];
    let a02 = combined.q[0][2];
    let a11 = combined.q[1][1];
    let a12 = combined.q[1][2];
    let a22 = combined.q[2][2];
    let b0 = -combined.q[0][3];
    let b1 = -combined.q[1][3];
    let b2 = -combined.q[2][3];
    let det = a00 * (a11 * a22 - a12 * a12) - a01 * (a01 * a22 - a12 * a02)
        + a02 * (a01 * a12 - a11 * a02);
    let opt = if det.abs() > 1e-10 {
        let inv_det = 1.0 / det;
        let x = inv_det
            * (b0 * (a11 * a22 - a12 * a12)
                + b1 * (a02 * a12 - a01 * a22)
                + b2 * (a01 * a12 - a11 * a02));
        let y = inv_det
            * (b0 * (a12 * a02 - a01 * a22)
                + b1 * (a00 * a22 - a02 * a02)
                + b2 * (a01 * a02 - a00 * a12));
        let z = inv_det
            * (b0 * (a01 * a12 - a02 * a11)
                + b1 * (a02 * a01 - a00 * a12)
                + b2 * (a00 * a11 - a01 * a01));
        [x, y, z]
    } else {
        mid
    };
    let error = combined.vertex_error(opt);
    (error, opt)
}
/// Compute the grid cell coordinates for vertex `v` with cell size `grid_size`.
pub fn cluster_id(v: [f64; 3], grid_size: f64) -> (i64, i64, i64) {
    (
        (v[0] / grid_size).floor() as i64,
        (v[1] / grid_size).floor() as i64,
        (v[2] / grid_size).floor() as i64,
    )
}
/// Decimate `mesh` by merging all vertices that fall inside the same voxel cell.
/// Each cluster is represented by the average position of its members.
/// Triangles with two or more vertices in the same cluster are discarded.
pub fn vertex_clustering(mesh: &SimpleMesh, grid_size: f64) -> SimpleMesh {
    let mut cell_to_cluster: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut cluster_sum: Vec<[f64; 3]> = Vec::new();
    let mut cluster_count: Vec<f64> = Vec::new();
    let mut vertex_to_cluster: Vec<usize> = Vec::with_capacity(mesh.vertex_count());
    for &v in &mesh.vertices {
        let cid = cluster_id(v, grid_size);
        let idx = cell_to_cluster.entry(cid).or_insert_with(|| {
            let new_idx = cluster_sum.len();
            cluster_sum.push([0.0; 3]);
            cluster_count.push(0.0);
            new_idx
        });
        cluster_sum[*idx][0] += v[0];
        cluster_sum[*idx][1] += v[1];
        cluster_sum[*idx][2] += v[2];
        cluster_count[*idx] += 1.0;
        vertex_to_cluster.push(*idx);
    }
    let mut out = SimpleMesh::new();
    for i in 0..cluster_sum.len() {
        let cnt = cluster_count[i];
        out.add_vertex([
            cluster_sum[i][0] / cnt,
            cluster_sum[i][1] / cnt,
            cluster_sum[i][2] / cnt,
        ]);
    }
    for &[a, b, c] in &mesh.triangles {
        let ca = vertex_to_cluster[a];
        let cb = vertex_to_cluster[b];
        let cc = vertex_to_cluster[c];
        if ca != cb && cb != cc && ca != cc {
            out.add_triangle(ca, cb, cc);
        }
    }
    out
}
/// Collapse edge (v1, v2): move v1 to `new_pos`, redirect all references to
/// v2 toward v1, then remove any resulting degenerate (zero-area) triangles.
/// Returns the index of the surviving vertex (v1, updated in-place).
pub fn collapse_edge(mesh: &mut SimpleMesh, v1: usize, v2: usize, new_pos: [f64; 3]) -> usize {
    mesh.vertices[v1] = new_pos;
    for tri in &mut mesh.triangles {
        for idx in tri.iter_mut() {
            if *idx == v2 {
                *idx = v1;
            }
        }
    }
    mesh.triangles
        .retain(|&[a, b, c]| a != b && b != c && a != c);
    v1
}
/// Return all edges that belong to exactly one triangle (boundary edges).
pub fn find_boundary_edges(mesh: &SimpleMesh) -> Vec<(usize, usize)> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for &[a, b, c] in &mesh.triangles {
        let edges = [
            (a.min(b), a.max(b)),
            (b.min(c), b.max(c)),
            (a.min(c), a.max(c)),
        ];
        for e in edges {
            *edge_count.entry(e).or_insert(0) += 1;
        }
    }
    edge_count
        .into_iter()
        .filter(|&(_, count)| count == 1)
        .map(|(e, _)| e)
        .collect()
}
/// Check whether collapsing edge (v1, v2) preserves manifold topology.
///
/// Uses the *link condition*: the collapse is safe if and only if the
/// intersection of the one-ring neighbourhoods of v1 and v2 equals the
/// set of vertices shared by triangles that contain both v1 and v2.
pub fn is_manifold_after_collapse(mesh: &SimpleMesh, v1: usize, v2: usize) -> bool {
    let mut ring1: HashSet<usize> = HashSet::new();
    let mut ring2: HashSet<usize> = HashSet::new();
    let mut shared: HashSet<usize> = HashSet::new();
    for &[a, b, c] in &mesh.triangles {
        let verts = [a, b, c];
        let has_v1 = verts.contains(&v1);
        let has_v2 = verts.contains(&v2);
        if has_v1 {
            for &v in &verts {
                if v != v1 {
                    ring1.insert(v);
                }
            }
        }
        if has_v2 {
            for &v in &verts {
                if v != v2 {
                    ring2.insert(v);
                }
            }
        }
        if has_v1 && has_v2 {
            for &v in &verts {
                if v != v1 && v != v2 {
                    shared.insert(v);
                }
            }
        }
    }
    ring1.remove(&v2);
    ring2.remove(&v1);
    let intersection: HashSet<usize> = ring1.intersection(&ring2).copied().collect();
    intersection == shared
}
/// Signed volume of the mesh via the divergence theorem.
/// Positive for outward-facing normals with CCW winding.
pub fn mesh_volume(mesh: &SimpleMesh) -> f64 {
    let mut vol = 0.0;
    for &[a, b, c] in &mesh.triangles {
        let va = mesh.vertices[a];
        let vb = mesh.vertices[b];
        let vc = mesh.vertices[c];
        vol += va[0] * (vb[1] * vc[2] - vc[1] * vb[2])
            + vb[0] * (vc[1] * va[2] - va[1] * vc[2])
            + vc[0] * (va[1] * vb[2] - vb[1] * va[2]);
    }
    vol / 6.0
}
/// Area-weighted centroid of the mesh surface.
pub fn mesh_centroid(mesh: &SimpleMesh) -> [f64; 3] {
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    let mut total_area = 0.0;
    for t_idx in 0..mesh.triangle_count() {
        let [a, b, c] = mesh.triangles[t_idx];
        let va = mesh.vertices[a];
        let vb = mesh.vertices[b];
        let vc = mesh.vertices[c];
        let area = mesh.triangle_area(t_idx);
        let tcx = (va[0] + vb[0] + vc[0]) / 3.0;
        let tcy = (va[1] + vb[1] + vc[1]) / 3.0;
        let tcz = (va[2] + vb[2] + vc[2]) / 3.0;
        cx += area * tcx;
        cy += area * tcy;
        cz += area * tcz;
        total_area += area;
    }
    if total_area < 1e-300 {
        return [0.0, 0.0, 0.0];
    }
    [cx / total_area, cy / total_area, cz / total_area]
}
/// Classify every vertex as *smooth* or *feature* based on the maximum dihedral
/// angle between its adjacent triangle faces.
///
/// A vertex is marked as a feature vertex when any pair of its adjacent faces
/// has a dihedral angle exceeding `threshold_rad`.
pub fn classify_feature_vertices(mesh: &SimpleMesh, threshold_rad: f64) -> Vec<bool> {
    let n = mesh.vertex_count();
    let mut is_feature = vec![false; n];
    let mut vtris: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (t_idx, &[a, b, c]) in mesh.triangles.iter().enumerate() {
        vtris[a].push(t_idx);
        vtris[b].push(t_idx);
        vtris[c].push(t_idx);
    }
    for v in 0..n {
        'outer: for &t1 in &vtris[v] {
            for &t2 in &vtris[v] {
                if t1 >= t2 {
                    continue;
                }
                let n1 = mesh.triangle_normal(t1);
                let n2 = mesh.triangle_normal(t2);
                let cos_ang = dot(n1, n2).clamp(-1.0, 1.0);
                let angle = cos_ang.acos();
                if angle > threshold_rad {
                    is_feature[v] = true;
                    break 'outer;
                }
            }
        }
    }
    is_feature
}
/// Simplify `mesh` using QEM edge collapses, preserving feature vertices
/// (those for which `is_feature[v] == true` will never be the removed vertex).
///
/// Reduces until at most `target_triangles` remain or no more safe collapses
/// are available.
pub fn simplify_with_features(mesh: &mut SimpleMesh, is_feature: &[bool], target_triangles: usize) {
    while mesh.triangle_count() > target_triangles {
        let quadrics = compute_vertex_quadrics(mesh);
        let mut edge_set: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        for &[a, b, c] in &mesh.triangles {
            edge_set.insert((a.min(b), a.max(b)));
            edge_set.insert((b.min(c), b.max(c)));
            edge_set.insert((a.min(c), a.max(c)));
        }
        let mut best_cost = f64::MAX;
        let mut best_v1 = usize::MAX;
        let mut best_v2 = usize::MAX;
        let mut best_pos = [0.0f64; 3];
        for (v1, v2) in &edge_set {
            let can_remove_v2 = !is_feature.get(*v2).copied().unwrap_or(false);
            let can_remove_v1 = !is_feature.get(*v1).copied().unwrap_or(false);
            if !can_remove_v2 && !can_remove_v1 {
                continue;
            }
            let (v_kept, v_rem) = if can_remove_v2 {
                (*v1, *v2)
            } else {
                (*v2, *v1)
            };
            let (cost, pos) = edge_collapse_cost(
                mesh.vertices[v_kept],
                mesh.vertices[v_rem],
                &quadrics[v_kept],
                &quadrics[v_rem],
            );
            if cost < best_cost {
                best_cost = cost;
                best_v1 = v_kept;
                best_v2 = v_rem;
                best_pos = pos;
            }
        }
        if best_v1 == usize::MAX {
            break;
        }
        collapse_edge(mesh, best_v1, best_v2, best_pos);
    }
}
/// Compute the Hausdorff-like error between two meshes (approximate).
///
/// For each vertex in `mesh_a`, finds the nearest vertex in `mesh_b` and
/// records the maximum distance.  This is a one-sided approximation of the
/// true Hausdorff distance.
pub fn one_sided_hausdorff(mesh_a: &SimpleMesh, mesh_b: &SimpleMesh) -> f64 {
    if mesh_a.vertex_count() == 0 || mesh_b.vertex_count() == 0 {
        return 0.0;
    }
    let mut max_dist = 0.0f64;
    for &va in &mesh_a.vertices {
        let mut min_d = f64::MAX;
        for &vb in &mesh_b.vertices {
            let d = {
                let dx = va[0] - vb[0];
                let dy = va[1] - vb[1];
                let dz = va[2] - vb[2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            };
            if d < min_d {
                min_d = d;
            }
        }
        if min_d > max_dist {
            max_dist = min_d;
        }
    }
    max_dist
}
/// Compute the average edge length of the mesh.
pub fn average_edge_length(mesh: &SimpleMesh) -> f64 {
    let mut total = 0.0;
    let mut count = 0usize;
    for &[a, b, c] in &mesh.triangles {
        for (i, j) in [(a, b), (b, c), (a, c)] {
            let vi = mesh.vertices[i];
            let vj = mesh.vertices[j];
            let d = length([vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]]);
            total += d;
            count += 1;
        }
    }
    if count == 0 {
        0.0
    } else {
        total / count as f64
    }
}
/// Compactness ratio: 36π V² / A³  (= 1 for a sphere, < 1 otherwise).
///
/// Returns `0.0` if area is effectively zero.
pub fn mesh_compactness(mesh: &SimpleMesh) -> f64 {
    let v = mesh_volume(mesh).abs();
    let a = mesh.surface_area();
    if a < 1e-300 {
        return 0.0;
    }
    36.0 * std::f64::consts::PI * v * v / (a * a * a)
}
/// Vertex clustering decimation: reduces mesh by merging vertices within cells
/// of a uniform grid, then rebuilding triangles.
///
/// The grid cell size `cell_size` controls the level of simplification.
pub fn vertex_clustering_decimate(mesh: &SimpleMesh, cell_size: f64) -> SimpleMesh {
    if mesh.vertices.is_empty() || cell_size <= 0.0 {
        return SimpleMesh::new();
    }
    let mut mn = mesh.vertices[0];
    let mut mx = mesh.vertices[0];
    for v in &mesh.vertices {
        for d in 0..3 {
            if v[d] < mn[d] {
                mn[d] = v[d];
            }
            if v[d] > mx[d] {
                mx[d] = v[d];
            }
        }
    }
    let mut cell_to_new_idx: HashMap<[i64; 3], usize> = HashMap::new();
    let mut vertex_remap: Vec<usize> = Vec::with_capacity(mesh.vertices.len());
    let mut new_vertices: Vec<[f64; 3]> = Vec::new();
    for v in &mesh.vertices {
        let key = [
            ((v[0] - mn[0]) / cell_size).floor() as i64,
            ((v[1] - mn[1]) / cell_size).floor() as i64,
            ((v[2] - mn[2]) / cell_size).floor() as i64,
        ];
        let idx = *cell_to_new_idx.entry(key).or_insert_with(|| {
            let new_idx = new_vertices.len();
            new_vertices.push(*v);
            new_idx
        });
        vertex_remap.push(idx);
    }
    let mut new_tris: Vec<[usize; 3]> = Vec::new();
    let mut seen: HashSet<[usize; 3]> = HashSet::new();
    for tri in &mesh.triangles {
        let a = vertex_remap[tri[0]];
        let b = vertex_remap[tri[1]];
        let c = vertex_remap[tri[2]];
        if a == b || b == c || a == c {
            continue;
        }
        let mut key = [a, b, c];
        key.sort_unstable();
        if seen.insert(key) {
            new_tris.push([a, b, c]);
        }
    }
    SimpleMesh {
        vertices: new_vertices,
        triangles: new_tris,
    }
}
/// Classify an edge between two triangles by their dihedral angle.
///
/// `n0`, `n1`: outward normals of the two triangles.
/// `crease_angle_deg`: threshold for classifying as a crease.
pub fn classify_edge(n0: [f64; 3], n1: [f64; 3], crease_angle_deg: f64) -> EdgeFeature {
    let cos_angle = dot3(normalize3(n0), normalize3(n1)).clamp(-1.0, 1.0);
    let angle_deg = cos_angle.acos().to_degrees();
    if angle_deg > crease_angle_deg {
        EdgeFeature::Crease
    } else {
        EdgeFeature::Smooth
    }
}
/// Feature-aware QEM cost: penalizes collapsing crease or boundary edges.
///
/// Returns the base QEM cost multiplied by a feature weight:
/// - Smooth: weight = 1.0
/// - Crease: weight = `crease_weight`
/// - Boundary: weight = `boundary_weight`
pub fn feature_aware_cost(
    base_cost: f64,
    feature: EdgeFeature,
    crease_weight: f64,
    boundary_weight: f64,
) -> f64 {
    let w = match feature {
        EdgeFeature::Smooth => 1.0,
        EdgeFeature::Crease => crease_weight,
        EdgeFeature::Boundary => boundary_weight,
    };
    base_cost * w
}
/// Collect decimation statistics before and after simplification.
pub fn collect_decimation_metrics(
    original: &SimpleMesh,
    reduced: &SimpleMesh,
    n_collapses: usize,
    total_cost: f64,
    max_cost: f64,
) -> DecimationMetrics {
    DecimationMetrics {
        original_vertices: original.vertex_count(),
        original_triangles: original.triangle_count(),
        reduced_vertices: reduced.vertex_count(),
        reduced_triangles: reduced.triangle_count(),
        n_collapses,
        total_qem_cost: total_cost,
        max_qem_cost: max_cost,
    }
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = dot3(v, v).sqrt();
    if len < 1e-15 {
        return [0.0, 0.0, 1.0];
    }
    [v[0] / len, v[1] / len, v[2] / len]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::decimation::DecimationStats;
    use crate::decimation::ProgressiveMesh;
    fn unit_tetrahedron() -> SimpleMesh {
        let mut m = SimpleMesh::new();
        let v0 = m.add_vertex([0.0, 0.0, 0.0]);
        let v1 = m.add_vertex([1.0, 0.0, 0.0]);
        let v2 = m.add_vertex([0.0, 1.0, 0.0]);
        let v3 = m.add_vertex([0.0, 0.0, 1.0]);
        m.add_triangle(v0, v2, v1);
        m.add_triangle(v0, v1, v3);
        m.add_triangle(v1, v2, v3);
        m.add_triangle(v0, v3, v2);
        m
    }
    fn single_triangle() -> SimpleMesh {
        let mut m = SimpleMesh::new();
        let v0 = m.add_vertex([0.0, 0.0, 0.0]);
        let v1 = m.add_vertex([1.0, 0.0, 0.0]);
        let v2 = m.add_vertex([0.0, 1.0, 0.0]);
        m.add_triangle(v0, v1, v2);
        m
    }
    fn unit_cube() -> SimpleMesh {
        let mut m = SimpleMesh::new();
        let verts: [[f64; 3]; 8] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        for v in &verts {
            m.add_vertex(*v);
        }
        m.add_triangle(0, 2, 1);
        m.add_triangle(0, 3, 2);
        m.add_triangle(4, 5, 6);
        m.add_triangle(4, 6, 7);
        m.add_triangle(0, 1, 5);
        m.add_triangle(0, 5, 4);
        m.add_triangle(2, 3, 7);
        m.add_triangle(2, 7, 6);
        m.add_triangle(0, 4, 7);
        m.add_triangle(0, 7, 3);
        m.add_triangle(1, 2, 6);
        m.add_triangle(1, 6, 5);
        m
    }
    #[test]
    fn test_mesh_new_empty() {
        let m = SimpleMesh::new();
        assert_eq!(m.vertex_count(), 0);
        assert_eq!(m.triangle_count(), 0);
    }
    #[test]
    fn test_add_vertex_returns_index() {
        let mut m = SimpleMesh::new();
        let i0 = m.add_vertex([1.0, 2.0, 3.0]);
        let i1 = m.add_vertex([4.0, 5.0, 6.0]);
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
        assert_eq!(m.vertex_count(), 2);
    }
    #[test]
    fn test_add_triangle_increments_count() {
        let mut m = single_triangle();
        assert_eq!(m.triangle_count(), 1);
        m.add_triangle(0, 1, 2);
        assert_eq!(m.triangle_count(), 2);
    }
    #[test]
    fn test_triangle_area_unit_right_triangle() {
        let m = single_triangle();
        let area = m.triangle_area(0);
        assert!((area - 0.5).abs() < 1e-12, "area={area}");
    }
    #[test]
    fn test_surface_area_tetrahedron() {
        let m = unit_tetrahedron();
        let expected = 3.0 * 0.5 + 3.0_f64.sqrt() / 2.0;
        let got = m.surface_area();
        assert!(
            (got - expected).abs() < 1e-10,
            "got={got}, expected={expected}"
        );
    }
    #[test]
    fn test_triangle_normal_z_up() {
        let m = single_triangle();
        let n = m.triangle_normal(0);
        assert!(n[2] > 0.9, "expected z-up normal, got {:?}", n);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "normal not unit length");
    }
    #[test]
    fn test_triangle_normal_unit_length() {
        let m = unit_tetrahedron();
        for i in 0..m.triangle_count() {
            let n = m.triangle_normal(i);
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-12, "face {i} normal not unit: {len}");
        }
    }
    #[test]
    fn test_quadric_zero() {
        let q = QuadricMatrix::zero();
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(q.q[i][j], 0.0);
            }
        }
    }
    #[test]
    fn test_quadric_from_plane_vertex_error_zero_on_plane() {
        let q = QuadricMatrix::from_plane(0.0, 0.0, 1.0, 0.0);
        let err = q.vertex_error([3.0, -2.0, 0.0]);
        assert!(err.abs() < 1e-12, "err={err}");
    }
    #[test]
    fn test_quadric_from_plane_vertex_error_nonzero_off_plane() {
        let q = QuadricMatrix::from_plane(0.0, 0.0, 1.0, 0.0);
        let err = q.vertex_error([0.0, 0.0, 1.0]);
        assert!((err - 1.0).abs() < 1e-12, "err={err}");
    }
    #[test]
    fn test_quadric_add() {
        let q1 = QuadricMatrix::from_plane(1.0, 0.0, 0.0, 0.0);
        let q2 = QuadricMatrix::from_plane(0.0, 1.0, 0.0, 0.0);
        let sum = q1.add(&q2);
        let err = sum.vertex_error([1.0, 1.0, 0.0]);
        assert!((err - 2.0).abs() < 1e-12, "err={err}");
    }
    #[test]
    fn test_compute_vertex_quadrics_count() {
        let m = unit_tetrahedron();
        let qs = compute_vertex_quadrics(&m);
        assert_eq!(qs.len(), m.vertex_count());
    }
    #[test]
    fn test_compute_vertex_quadrics_nonzero() {
        let m = single_triangle();
        let qs = compute_vertex_quadrics(&m);
        let total: f64 = qs[0].q.iter().flatten().map(|x| x.abs()).sum();
        assert!(total > 0.0);
    }
    #[test]
    fn test_edge_collapse_cost_midpoint_fallback() {
        let q1 = QuadricMatrix::zero();
        let q2 = QuadricMatrix::zero();
        let v1 = [0.0, 0.0, 0.0];
        let v2 = [2.0, 0.0, 0.0];
        let (error, pos) = edge_collapse_cost(v1, v2, &q1, &q2);
        assert!(error.abs() < 1e-12);
        assert!((pos[0] - 1.0).abs() < 1e-12, "pos={:?}", pos);
    }
    #[test]
    fn test_edge_collapse_cost_returns_finite() {
        let m = unit_tetrahedron();
        let qs = compute_vertex_quadrics(&m);
        let (err, pos) = edge_collapse_cost(m.vertices[0], m.vertices[1], &qs[0], &qs[1]);
        assert!(err.is_finite(), "error not finite");
        assert!(pos.iter().all(|x| x.is_finite()), "pos not finite");
    }
    #[test]
    fn test_cluster_id_basic() {
        let id = cluster_id([0.5, 1.5, 2.5], 1.0);
        assert_eq!(id, (0, 1, 2));
    }
    #[test]
    fn test_cluster_id_negative() {
        let id = cluster_id([-0.5, -1.5, -2.5], 1.0);
        assert_eq!(id, (-1, -2, -3));
    }
    #[test]
    fn test_vertex_clustering_reduces_vertices() {
        let mut m = SimpleMesh::new();
        let eps = 0.01;
        let v0 = m.add_vertex([0.0, 0.0, 0.0]);
        let v1 = m.add_vertex([1.0, 0.0, 0.0]);
        let v2 = m.add_vertex([0.0, 1.0, 0.0]);
        let v3 = m.add_vertex([0.0 + eps, 0.0, 0.0]);
        let v4 = m.add_vertex([1.0 + eps, 0.0, 0.0]);
        let v5 = m.add_vertex([0.0 + eps, 1.0, 0.0]);
        m.add_triangle(v0, v1, v2);
        m.add_triangle(v3, v4, v5);
        let simplified = vertex_clustering(&m, 0.1);
        assert!(simplified.vertex_count() < m.vertex_count());
    }
    #[test]
    fn test_vertex_clustering_no_degenerate_triangles() {
        let mut m = SimpleMesh::new();
        let v0 = m.add_vertex([0.1, 0.1, 0.1]);
        let v1 = m.add_vertex([0.2, 0.1, 0.1]);
        let v2 = m.add_vertex([0.1, 0.2, 0.1]);
        m.add_triangle(v0, v1, v2);
        let simplified = vertex_clustering(&m, 1.0);
        assert_eq!(simplified.triangle_count(), 0);
    }
    #[test]
    fn test_collapse_edge_reduces_triangles() {
        let mut m = unit_tetrahedron();
        let before_tris = m.triangle_count();
        collapse_edge(&mut m, 0, 1, [0.5, 0.0, 0.0]);
        assert!(m.triangle_count() < before_tris);
    }
    #[test]
    fn test_collapse_edge_no_degenerate_triangles_remain() {
        let mut m = single_triangle();
        collapse_edge(&mut m, 0, 1, [0.5, 0.0, 0.0]);
        for &[a, b, c] in &m.triangles {
            assert!(a != b && b != c && a != c, "degenerate triangle found");
        }
    }
    #[test]
    fn test_find_boundary_edges_single_triangle() {
        let m = single_triangle();
        let boundary = find_boundary_edges(&m);
        assert_eq!(boundary.len(), 3);
    }
    #[test]
    fn test_find_boundary_edges_closed_tetrahedron() {
        let m = unit_tetrahedron();
        let boundary = find_boundary_edges(&m);
        assert_eq!(boundary.len(), 0);
    }
    #[test]
    fn test_is_manifold_after_collapse_tetrahedron_edge() {
        let m = unit_tetrahedron();
        let result = is_manifold_after_collapse(&m, 0, 1);
        assert!(result, "manifold check failed on tetrahedron edge");
    }
    #[test]
    fn test_mesh_volume_unit_cube() {
        let m = unit_cube();
        let vol = mesh_volume(&m).abs();
        assert!((vol - 1.0).abs() < 1e-10, "cube volume={vol}");
    }
    #[test]
    fn test_mesh_volume_single_triangle_zero() {
        let m = single_triangle();
        let vol = mesh_volume(&m);
        assert!(vol.abs() < 1e-12, "single triangle volume={vol}");
    }
    #[test]
    fn test_mesh_centroid_single_triangle() {
        let m = single_triangle();
        let c = mesh_centroid(&m);
        assert!((c[0] - 1.0 / 3.0).abs() < 1e-12, "cx={}", c[0]);
        assert!((c[1] - 1.0 / 3.0).abs() < 1e-12, "cy={}", c[1]);
        assert!(c[2].abs() < 1e-12, "cz={}", c[2]);
    }
    #[test]
    fn test_mesh_centroid_unit_cube() {
        let m = unit_cube();
        let c = mesh_centroid(&m);
        assert!((c[0] - 0.5).abs() < 1e-10, "cx={}", c[0]);
        assert!((c[1] - 0.5).abs() < 1e-10, "cy={}", c[1]);
        assert!((c[2] - 0.5).abs() < 1e-10, "cz={}", c[2]);
    }
    #[test]
    fn test_progressive_mesh_collapse_one_reduces_triangles() {
        let m = unit_tetrahedron();
        let before = m.triangle_count();
        let mut pm = ProgressiveMesh::new(m);
        let ok = pm.collapse_one();
        assert!(ok, "collapse_one should succeed on a tetrahedron");
        assert!(
            pm.mesh.triangle_count() < before,
            "triangle count should decrease"
        );
    }
    #[test]
    fn test_progressive_mesh_decimate_to() {
        let m = unit_cube();
        let mut pm = ProgressiveMesh::new(m);
        pm.decimate_to(4);
        assert!(
            pm.mesh.triangle_count() <= 4,
            "got {}",
            pm.mesh.triangle_count()
        );
    }
    #[test]
    fn test_progressive_mesh_history_recorded() {
        let m = unit_tetrahedron();
        let mut pm = ProgressiveMesh::new(m);
        pm.collapse_one();
        assert_eq!(pm.history.len(), 1);
        let rec = &pm.history[0];
        assert!(rec.old_pos.iter().all(|x| x.is_finite()));
        assert!(rec.new_pos.iter().all(|x| x.is_finite()));
    }
    #[test]
    fn test_progressive_mesh_stops_at_minimum() {
        let m = single_triangle();
        let mut pm = ProgressiveMesh::new(m);
        let ok = pm.collapse_one();
        assert!(!ok, "should not collapse a 1-triangle mesh");
    }
    #[test]
    fn test_classify_feature_vertices_tetrahedron() {
        let m = unit_tetrahedron();
        let feat = classify_feature_vertices(&m, 0.01);
        assert!(
            feat.iter().any(|&f| f),
            "at least one feature vertex expected"
        );
    }
    #[test]
    fn test_classify_feature_vertices_flat_mesh() {
        let mut m = SimpleMesh::new();
        let v0 = m.add_vertex([0.0, 0.0, 0.0]);
        let v1 = m.add_vertex([1.0, 0.0, 0.0]);
        let v2 = m.add_vertex([0.0, 1.0, 0.0]);
        let v3 = m.add_vertex([1.0, 1.0, 0.0]);
        m.add_triangle(v0, v1, v2);
        m.add_triangle(v1, v3, v2);
        let feat = classify_feature_vertices(&m, 0.1);
        assert!(
            feat.iter().all(|&f| !f),
            "flat mesh should have no feature vertices"
        );
    }
    #[test]
    fn test_simplify_with_features_preserves_feature_vertex() {
        let m = unit_cube();
        let n = m.vertex_count();
        let mut is_feature = vec![false; n];
        is_feature[0] = true;
        let mut mesh = m;
        simplify_with_features(&mut mesh, &is_feature, 2);
        let v0_used = mesh.triangles.iter().any(|t| t.contains(&0));
        assert!(v0_used, "feature vertex 0 should survive decimation");
    }
    #[test]
    fn test_simplify_with_features_reduces_triangles() {
        let m = unit_cube();
        let is_feature = vec![false; m.vertex_count()];
        let target = 4;
        let mut mesh = m;
        simplify_with_features(&mut mesh, &is_feature, target);
        assert!(mesh.triangle_count() <= target || mesh.triangle_count() < 12);
    }
    #[test]
    fn test_decimation_stats_reduction_ratio() {
        let m = unit_cube();
        let orig_v = m.vertex_count();
        let orig_t = m.triangle_count();
        let mut pm = ProgressiveMesh::new(m);
        pm.decimate_to(6);
        let stats = DecimationStats::compute(orig_v, orig_t, &pm.mesh, 0.0);
        assert!(stats.reduction_ratio >= 0.0 && stats.reduction_ratio <= 1.0);
        assert_eq!(stats.original_triangles, orig_t);
        assert_eq!(stats.result_triangles, pm.mesh.triangle_count());
    }
    #[test]
    fn test_decimation_stats_no_change() {
        let m = single_triangle();
        let stats = DecimationStats::compute(m.vertex_count(), m.triangle_count(), &m, 0.0);
        assert!((stats.reduction_ratio - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_one_sided_hausdorff_same_mesh_zero() {
        let m = single_triangle();
        let d = one_sided_hausdorff(&m, &m);
        assert!(d < 1e-12, "same mesh: d={d}");
    }
    #[test]
    fn test_one_sided_hausdorff_offset_mesh() {
        let mut m1 = SimpleMesh::new();
        m1.add_vertex([0.0, 0.0, 0.0]);
        let mut m2 = SimpleMesh::new();
        m2.add_vertex([3.0, 0.0, 0.0]);
        let d = one_sided_hausdorff(&m1, &m2);
        assert!((d - 3.0).abs() < 1e-12, "d={d}");
    }
    #[test]
    fn test_one_sided_hausdorff_empty_returns_zero() {
        let m = single_triangle();
        let empty = SimpleMesh::new();
        let d = one_sided_hausdorff(&m, &empty);
        assert_eq!(d, 0.0);
    }
    #[test]
    fn test_average_edge_length_unit_right_triangle() {
        let m = single_triangle();
        let avg = average_edge_length(&m);
        let expected = (2.0 + 2.0_f64.sqrt()) / 3.0;
        assert!((avg - expected).abs() < 1e-10, "avg={avg}");
    }
    #[test]
    fn test_average_edge_length_empty_mesh_zero() {
        let m = SimpleMesh::new();
        assert_eq!(average_edge_length(&m), 0.0);
    }
    #[test]
    fn test_mesh_compactness_single_triangle_zero() {
        let m = single_triangle();
        let c = mesh_compactness(&m);
        assert!(c >= 0.0);
    }
    #[test]
    fn test_mesh_compactness_cube_positive() {
        let m = unit_cube();
        let c = mesh_compactness(&m);
        assert!(c > 0.0 && c < 1.0, "compactness={c}");
    }
}
