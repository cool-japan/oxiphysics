// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{HashMap, HashSet, VecDeque};

use super::types::{
    CurveSkeleton, HalfEdgeMesh, MedialAxisPoint, MorseCriticalType, MorseVertex, NonManifoldEdge,
    NonManifoldVertex, PersistenceDiagram, SimplicialComplex,
};

/// Add two 3-vectors.
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors.
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector.
#[inline]
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Dot product.
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Length of a 3-vector.
#[inline]
pub(super) fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalize a 3-vector.
#[inline]
pub(super) fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}
/// Cross product.
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Compute a 0-dimensional persistence diagram from a scalar function on vertices.
///
/// Uses union-find to track connected component birth/death.
pub fn compute_0d_persistence(
    vertex_values: &[f64],
    edges: &[(usize, usize)],
) -> PersistenceDiagram {
    let n = vertex_values.len();
    let mut diagram = PersistenceDiagram::new();
    let mut sorted_verts: Vec<(usize, f64)> = vertex_values.iter().copied().enumerate().collect();
    sorted_verts.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank = vec![0usize; n];
    let mut birth = vertex_values.to_vec();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    fn union(parent: &mut Vec<usize>, rank: &mut [usize], x: usize, y: usize) -> (usize, usize) {
        let rx = find(parent, x);
        let ry = find(parent, y);
        if rx == ry {
            return (rx, rx);
        }
        if rank[rx] < rank[ry] {
            parent[rx] = ry;
            (rx, ry)
        } else if rank[rx] > rank[ry] {
            parent[ry] = rx;
            (ry, rx)
        } else {
            parent[ry] = rx;
            rank[rx] += 1;
            (ry, rx)
        }
    }
    let mut sorted_edges: Vec<(usize, usize, f64)> = edges
        .iter()
        .map(|&(u, v)| {
            let val = vertex_values[u].max(vertex_values[v]);
            (u, v, val)
        })
        .collect();
    sorted_edges.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    for &(vi, val) in &sorted_verts {
        birth[vi] = val;
    }
    for &(u, v, death_val) in &sorted_edges {
        let ru = find(&mut parent, u);
        let rv = find(&mut parent, v);
        if ru != rv {
            let birth_val = birth[ru].min(birth[rv]);
            let dying = if birth[ru] > birth[rv] { ru } else { rv };
            let (killed, _survived) = union(&mut parent, &mut rank, u, v);
            let _ = (killed, birth_val, dying);
            diagram.add_pair(0, birth[dying], death_val);
        }
    }
    let min_val = vertex_values.iter().copied().fold(f64::MAX, f64::min);
    diagram.add_pair(0, min_val, f64::MAX);
    diagram
}
/// Compute the Euler characteristic from vertex, edge, and face counts.
///
/// χ = V - E + F
pub fn euler_characteristic_vef(v: usize, e: usize, f: usize) -> i64 {
    v as i64 - e as i64 + f as i64
}
/// Compute the genus of a closed orientable surface.
///
/// g = (2 - χ) / 2
pub fn genus_from_euler(chi: i64) -> i64 {
    (2 - chi) / 2
}
/// Compute all Betti numbers for a simplicial complex (b0, b1, b2).
///
/// Uses the relation: χ = b0 - b1 + b2 and the counts of simplices.
pub fn betti_numbers_approx(sc: &SimplicialComplex) -> (usize, usize, usize) {
    let chi = sc.euler_characteristic();
    let c0 = sc.count(0) as i64;
    let c1 = sc.count(1) as i64;
    let _c2 = sc.count(2) as i64;
    let b0 = 1_i64.max(c0 - c1 + chi);
    let b2 = if chi < 0 { 0_i64 } else { 1_i64 };
    let b1 = (b0 + b2 - chi).max(0);
    (b0.max(0) as usize, b1 as usize, b2 as usize)
}
/// Classify critical points of a Morse function on a triangle mesh.
pub fn classify_morse_critical_points(mesh: &HalfEdgeMesh, values: &[f64]) -> Vec<MorseVertex> {
    let n = mesh.num_vertices();
    let mut result = Vec::with_capacity(n);
    for vi in 0..n {
        let neighbors = mesh.vertex_neighbors(vi);
        if neighbors.is_empty() {
            result.push(MorseVertex {
                idx: vi,
                value: values[vi],
                critical_type: MorseCriticalType::Regular,
            });
            continue;
        }
        let v_val = values[vi];
        let n_above = neighbors.iter().filter(|&&nj| values[nj] > v_val).count();
        let n_below = neighbors.iter().filter(|&&nj| values[nj] < v_val).count();
        let ctype = if n_above == 0 {
            MorseCriticalType::Maximum
        } else if n_below == 0 {
            MorseCriticalType::Minimum
        } else {
            let sign_changes = neighbors
                .windows(2)
                .filter(|w| (values[w[0]] > v_val) != (values[w[1]] > v_val))
                .count();
            if sign_changes > 2 {
                MorseCriticalType::Saddle
            } else {
                MorseCriticalType::Regular
            }
        };
        result.push(MorseVertex {
            idx: vi,
            value: v_val,
            critical_type: ctype,
        });
    }
    result
}
/// Trace a descending gradient flow path from `start_vertex`.
///
/// Returns the sequence of vertex indices along the descending path.
pub fn trace_descending_manifold(
    mesh: &HalfEdgeMesh,
    values: &[f64],
    start_vertex: usize,
    max_steps: usize,
) -> Vec<usize> {
    let mut path = vec![start_vertex];
    let mut cur = start_vertex;
    for _ in 0..max_steps {
        let neighbors = mesh.vertex_neighbors(cur);
        if neighbors.is_empty() {
            break;
        }
        let next = neighbors.iter().copied().min_by(|&a, &b| {
            values[a]
                .partial_cmp(&values[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if let Some(nxt) = next {
            if values[nxt] >= values[cur] {
                break;
            }
            path.push(nxt);
            cur = nxt;
        } else {
            break;
        }
    }
    path
}
/// Extract a simplified curve skeleton from a mesh using the mean curvature flow approach.
///
/// This is a very simplified version: contracts the mesh toward its medial axis
/// by iteratively moving vertices toward the centroid of their neighbors.
pub fn extract_curve_skeleton(
    mesh: &HalfEdgeMesh,
    iterations: usize,
    contraction_rate: f64,
) -> CurveSkeleton {
    let mut positions: Vec<[f64; 3]> = mesh.vertices.clone();
    let n = positions.len();
    for _ in 0..iterations {
        let old = positions.clone();
        for vi in 0..n {
            let neighbors = mesh.vertex_neighbors(vi);
            if neighbors.is_empty() {
                continue;
            }
            let centroid = neighbors
                .iter()
                .fold([0.0f64; 3], |acc, &nj| add3(acc, old[nj]));
            let centroid = scale3(centroid, 1.0 / neighbors.len() as f64);
            positions[vi] = add3(
                scale3(old[vi], 1.0 - contraction_rate),
                scale3(centroid, contraction_rate),
            );
        }
    }
    let mut skeleton = CurveSkeleton::new();
    for &p in &positions {
        skeleton.add_node(p);
    }
    for he in &mesh.half_edges {
        if he.face == usize::MAX {
            continue;
        }
        let a = he.origin;
        let b = mesh.half_edges[he.next].origin;
        if a < b {
            skeleton.add_edge(a, b);
        }
    }
    skeleton
}
/// Compute an approximate medial axis from a point cloud.
///
/// Uses the Voronoi-based approach: poles of Voronoi cells are medial axis candidates.
/// This simplified version uses a brute-force approach.
pub fn compute_medial_axis_approx(
    surface_vertices: &[[f64; 3]],
    num_candidates: usize,
) -> Vec<MedialAxisPoint> {
    let n = surface_vertices.len();
    if n < 4 {
        return Vec::new();
    }
    let mut candidates = Vec::new();
    let step = (n / num_candidates.max(1)).max(1);
    let mut i = 0;
    while i < n && candidates.len() < num_candidates {
        let j = (i + step) % n;
        let k = (i + 2 * step) % n;
        if i == j || j == k || i == k {
            i += 1;
            continue;
        }
        let a = surface_vertices[i];
        let b = surface_vertices[j];
        let c = surface_vertices[k];
        let ab = sub3(b, a);
        let ac = sub3(c, a);
        let ab_len2 = dot3(ab, ab);
        let ac_len2 = dot3(ac, ac);
        let ab_ac = dot3(ab, ac);
        let denom = 2.0 * (ab_len2 * ac_len2 - ab_ac * ab_ac);
        if denom.abs() < 1e-12 {
            i += 1;
            continue;
        }
        let s = (ab_len2 * ac_len2 - ab_ac * dot3(ab, ac)) / denom;
        let t = (ac_len2 * ab_len2 - ab_ac * dot3(ac, ab)) / denom;
        let center = add3(a, add3(scale3(ab, s), scale3(ac, t)));
        let r = len3(sub3(center, a));
        candidates.push(MedialAxisPoint {
            position: center,
            radius: r,
            closest_verts: [i, j],
        });
        i += 1;
    }
    candidates
}
/// Detect non-manifold edges in a half-edge mesh.
pub fn detect_non_manifold_edges(mesh: &HalfEdgeMesh) -> Vec<NonManifoldEdge> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for he in &mesh.half_edges {
        if he.face == usize::MAX {
            continue;
        }
        let a = he.origin;
        let b = mesh.half_edges[he.next].origin;
        let key = if a < b { (a, b) } else { (b, a) };
        *edge_count.entry(key).or_insert(0) += 1;
    }
    edge_count
        .into_iter()
        .filter(|(_, c)| *c > 2)
        .map(|((v0, v1), c)| NonManifoldEdge {
            v0,
            v1,
            face_count: c,
        })
        .collect()
}
/// Detect non-manifold vertices using face connectivity.
pub fn detect_non_manifold_vertices(mesh: &HalfEdgeMesh) -> Vec<NonManifoldVertex> {
    let nv = mesh.num_vertices();
    let mut result = Vec::new();
    let mut vf: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (fi, _) in mesh.faces.iter().enumerate() {
        for he_idx in mesh.face_half_edges(fi) {
            let v = mesh.half_edges[he_idx].origin;
            if v < nv {
                vf[v].push(fi);
            }
        }
    }
    for (v, faces) in vf.iter().enumerate().take(nv) {
        if faces.len() <= 1 {
            continue;
        }
        let mut visited: HashSet<usize> = HashSet::new();
        let mut fan_count = 0;
        for &start_fi in faces {
            if visited.contains(&start_fi) {
                continue;
            }
            fan_count += 1;
            let mut queue = VecDeque::new();
            queue.push_back(start_fi);
            while let Some(fi) = queue.pop_front() {
                if !visited.insert(fi) {
                    continue;
                }
                for he_idx in mesh.face_half_edges(fi) {
                    if mesh.half_edges[he_idx].twin == usize::MAX {
                        continue;
                    }
                    let twin_he = mesh.half_edges[he_idx].twin;
                    let adj_fi = mesh.half_edges[twin_he].face;
                    if adj_fi != usize::MAX && faces.contains(&adj_fi) && !visited.contains(&adj_fi)
                    {
                        queue.push_back(adj_fi);
                    }
                }
            }
        }
        if fan_count > 1 {
            result.push(NonManifoldVertex { v, fan_count });
        }
    }
    result
}
/// Extract boundary loops from a half-edge mesh.
///
/// Returns each boundary loop as an ordered list of vertex indices.
pub fn extract_boundary_loops(mesh: &HalfEdgeMesh) -> Vec<Vec<usize>> {
    let boundary_hes = mesh.boundary_half_edges();
    let mut remaining: HashSet<usize> = boundary_hes.iter().copied().collect();
    let mut loops = Vec::new();
    while !remaining.is_empty() {
        let &start_he = remaining
            .iter()
            .next()
            .expect("iterator should have elements");
        let mut loop_verts = Vec::new();
        let mut cur = start_he;
        loop {
            loop_verts.push(mesh.half_edges[cur].origin);
            remaining.remove(&cur);
            let dst = mesh.half_edges[mesh.half_edges[cur].next].origin;
            let next_candidate = remaining
                .iter()
                .find(|&&he| mesh.half_edges[he].origin == dst)
                .copied();
            match next_candidate {
                Some(next) => cur = next,
                None => break,
            }
            if loop_verts.len() > mesh.num_vertices() {
                break;
            }
        }
        if !loop_verts.is_empty() {
            loops.push(loop_verts);
        }
    }
    loops
}
/// Compute connected components of the vertex graph using BFS.
///
/// Returns a vector of length `num_vertices` with component IDs.
pub fn connected_components_vertices(num_vertices: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); num_vertices];
    for &(a, b) in edges {
        if a < num_vertices && b < num_vertices {
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    let mut comp = vec![usize::MAX; num_vertices];
    let mut comp_id = 0;
    for start in 0..num_vertices {
        if comp[start] != usize::MAX {
            continue;
        }
        let mut queue = VecDeque::new();
        queue.push_back(start);
        comp[start] = comp_id;
        while let Some(v) = queue.pop_front() {
            for &nv in &adj[v] {
                if comp[nv] == usize::MAX {
                    comp[nv] = comp_id;
                    queue.push_back(nv);
                }
            }
        }
        comp_id += 1;
    }
    comp
}
/// Count the number of connected components.
pub fn count_connected_components(num_vertices: usize, edges: &[(usize, usize)]) -> usize {
    let comps = connected_components_vertices(num_vertices, edges);
    let max_comp = comps.iter().copied().filter(|&c| c != usize::MAX).max();
    max_comp.map_or(0, |m| m + 1)
}
/// Count connected components (β₀) from an adjacency edge list.
///
/// Uses union-find (DSU). Returns the number of connected components.
pub fn betti_0(adjacency: &[(usize, usize)], n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    for &(a, b) in adjacency {
        let ra = find(&mut parent, a.min(n - 1));
        let rb = find(&mut parent, b.min(n - 1));
        if ra != rb {
            parent[ra] = rb;
        }
    }
    let mut roots = std::collections::HashSet::new();
    for i in 0..n {
        roots.insert(find(&mut parent, i));
    }
    roots.len()
}
/// Estimate β₁ (first Betti number) from Euler characteristic:
///
/// β₁ = β₀ + E - V (for a surface, β₁ = E - V + 1 when β₀ = 1).
pub fn betti_1_estimate(n_v: usize, n_e: usize, n_f: usize) -> i32 {
    let chi = n_v as i32 - n_e as i32 + n_f as i32;
    2 - chi
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::BirthDeathPair;
    use crate::topology::IsoCurve;
    use crate::topology::MeshSimplification;
    use crate::topology::MeshTopology;
    use crate::topology::PolygonOffset;
    use crate::topology::QuadEdgeMesh;
    use crate::topology::Simplex;
    use crate::topology::TopologicalSurgery;
    use crate::topology::WingedEdgeMesh;
    /// Build a simple tetrahedron mesh for testing.
    fn make_tet_mesh() -> HalfEdgeMesh {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let faces = vec![[0, 1, 2], [0, 3, 1], [1, 3, 2], [0, 2, 3]];
        HalfEdgeMesh::from_triangles(verts, &faces)
    }
    /// Build a single triangle.
    fn make_triangle_mesh() -> HalfEdgeMesh {
        let verts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let faces = vec![[0, 1, 2]];
        HalfEdgeMesh::from_triangles(verts, &faces)
    }
    #[test]
    fn test_half_edge_mesh_tet_vertex_count() {
        let m = make_tet_mesh();
        assert_eq!(m.num_vertices(), 4);
    }
    #[test]
    fn test_half_edge_mesh_tet_face_count() {
        let m = make_tet_mesh();
        assert_eq!(m.num_faces(), 4);
    }
    #[test]
    fn test_half_edge_mesh_triangle_is_open() {
        let m = make_triangle_mesh();
        assert!(!m.is_closed(), "single triangle should have boundary");
    }
    #[test]
    fn test_half_edge_mesh_tet_is_closed() {
        let m = make_tet_mesh();
        assert!(m.is_closed(), "tetrahedron should be closed");
    }
    #[test]
    fn test_half_edge_mesh_tet_euler() {
        let m = make_tet_mesh();
        let chi = m.euler_characteristic();
        assert_eq!(chi, 2, "tetrahedron Euler char should be 2, got {chi}");
    }
    #[test]
    fn test_half_edge_mesh_tet_genus() {
        let m = make_tet_mesh();
        assert_eq!(m.genus(), 0, "tetrahedron genus should be 0");
    }
    #[test]
    fn test_half_edge_mesh_is_manifold() {
        let m = make_tet_mesh();
        assert!(m.is_manifold(), "tetrahedron should be manifold");
    }
    #[test]
    fn test_half_edge_face_vertices() {
        let m = make_triangle_mesh();
        let verts = m.face_vertices(0);
        assert_eq!(verts.len(), 3);
    }
    #[test]
    fn test_half_edge_face_normal_upward() {
        let m = make_triangle_mesh();
        let n = m.face_normal(0);
        assert!(
            n[2].abs() > 0.5,
            "normal z component should dominate: n={n:?}"
        );
    }
    #[test]
    fn test_half_edge_face_area() {
        let m = make_triangle_mesh();
        let area = m.face_area(0);
        assert!((area - 0.5).abs() < 1e-8, "area={area}");
    }
    #[test]
    fn test_half_edge_total_area() {
        let m = make_triangle_mesh();
        assert!((m.total_area() - 0.5).abs() < 1e-8);
    }
    #[test]
    fn test_half_edge_vertex_normals() {
        let m = make_tet_mesh();
        let norms = m.vertex_normals();
        assert_eq!(norms.len(), 4);
        for n in &norms {
            let l = len3(*n);
            assert!(
                (l - 1.0).abs() < 1e-6 || l < 1e-10,
                "normal should be unit or zero: l={l}"
            );
        }
    }
    #[test]
    fn test_winged_edge_add_vertex() {
        let mut we = WingedEdgeMesh::new();
        let v = we.add_vertex([1.0, 2.0, 3.0]);
        assert_eq!(v, 0);
        assert_eq!(we.vertices.len(), 1);
    }
    #[test]
    fn test_winged_edge_add_edge() {
        let mut we = WingedEdgeMesh::new();
        we.add_vertex([0.0; 3]);
        we.add_vertex([1.0, 0.0, 0.0]);
        let e = we.add_edge(0, 1);
        assert_eq!(e, 0);
        assert_eq!(we.num_edges(), 1);
    }
    #[test]
    fn test_quad_edge_make_edge() {
        let mut qe = QuadEdgeMesh::new();
        let v0 = qe.add_vertex([0.0; 3]);
        let v1 = qe.add_vertex([1.0, 0.0, 0.0]);
        let e = qe.make_edge(v0, v1);
        assert_eq!(e, 0);
        assert_eq!(qe.num_edges(), 1);
    }
    #[test]
    fn test_simplicial_complex_add_triangle() {
        let mut sc = SimplicialComplex::new(2);
        sc.add_simplex(Simplex::new(vec![0, 1, 2]));
        assert_eq!(sc.count(0), 3);
        assert_eq!(sc.count(1), 3);
        assert_eq!(sc.count(2), 1);
    }
    #[test]
    fn test_simplicial_complex_euler_triangle() {
        let mut sc = SimplicialComplex::new(2);
        sc.add_simplex(Simplex::new(vec![0, 1, 2]));
        let chi = sc.euler_characteristic();
        assert_eq!(chi, 1, "triangle chi = {chi}");
    }
    #[test]
    fn test_simplicial_complex_is_valid() {
        let mut sc = SimplicialComplex::new(2);
        sc.add_simplex(Simplex::new(vec![0, 1, 2]));
        assert!(sc.is_valid());
    }
    #[test]
    fn test_simplex_boundary_faces() {
        let s = Simplex::new(vec![0, 1, 2]);
        let faces = s.boundary_faces();
        assert_eq!(faces.len(), 3);
    }
    #[test]
    fn test_simplex_is_face_of() {
        let edge = Simplex::new(vec![0, 1]);
        let tri = Simplex::new(vec![0, 1, 2]);
        assert!(edge.is_face_of(&tri));
    }
    #[test]
    fn test_simplex_dim() {
        assert_eq!(Simplex::new(vec![0]).dim(), 0);
        assert_eq!(Simplex::new(vec![0, 1]).dim(), 1);
        assert_eq!(Simplex::new(vec![0, 1, 2]).dim(), 2);
        assert_eq!(Simplex::new(vec![0, 1, 2, 3]).dim(), 3);
    }
    #[test]
    fn test_persistence_diagram_add_pair() {
        let mut pd = PersistenceDiagram::new();
        pd.add_pair(0, 0.0, 1.0);
        assert_eq!(pd.pairs.len(), 1);
    }
    #[test]
    fn test_persistence_betti_at() {
        let mut pd = PersistenceDiagram::new();
        pd.add_pair(0, 0.0, f64::MAX);
        pd.add_pair(0, 0.5, 2.0);
        assert_eq!(pd.betti_at(0, 1.0), 2);
        assert_eq!(pd.betti_at(0, 3.0), 1);
    }
    #[test]
    fn test_persistence_max_persistence() {
        let mut pd = PersistenceDiagram::new();
        pd.add_pair(0, 0.0, 3.0);
        pd.add_pair(0, 1.0, 4.0);
        assert!((pd.max_persistence() - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_birth_death_pair_persistence() {
        let p = BirthDeathPair {
            dim: 0,
            birth: 1.0,
            death: 3.0,
        };
        assert!((p.persistence() - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_essential_class_infinite_persistence() {
        let p = BirthDeathPair {
            dim: 0,
            birth: 0.0,
            death: f64::MAX,
        };
        assert!(p.is_essential());
    }
    #[test]
    fn test_0d_persistence_line_graph() {
        let values = [0.0, 1.0, 2.0, 3.0];
        let edges = [(0, 1), (1, 2), (2, 3)];
        let pd = compute_0d_persistence(&values, &edges);
        assert!(!pd.pairs.is_empty());
    }
    #[test]
    fn test_euler_characteristic_sphere() {
        let chi = euler_characteristic_vef(4, 6, 4);
        assert_eq!(chi, 2);
    }
    #[test]
    fn test_genus_from_euler_sphere() {
        assert_eq!(genus_from_euler(2), 0);
    }
    #[test]
    fn test_genus_from_euler_torus() {
        assert_eq!(genus_from_euler(0), 1);
    }
    #[test]
    fn test_edge_collapse() {
        let mut m = make_tet_mesh();
        let mut surgery = TopologicalSurgery::new(&mut m);
        assert!(surgery.edge_collapse(0));
    }
    #[test]
    fn test_edge_split() {
        let mut m = make_triangle_mesh();
        let orig_nv = m.num_vertices();
        let new_v = {
            let mut surgery = TopologicalSurgery::new(&mut m);
            surgery.edge_split(0)
        };
        assert_ne!(new_v, usize::MAX);
        assert_eq!(m.num_vertices(), orig_nv + 1);
    }
    #[test]
    fn test_attach_handle_needs_boundary() {
        let mut m = make_triangle_mesh();
        let mut surgery = TopologicalSurgery::new(&mut m);
        let result = surgery.attach_handle(0, 1);
        assert!(result);
    }
    #[test]
    fn test_morse_critical_classification() {
        let m = make_tet_mesh();
        let values = [0.0, 1.0, 2.0, 3.0];
        let crits = classify_morse_critical_points(&m, &values);
        assert_eq!(crits.len(), 4);
        assert_eq!(crits[0].critical_type, MorseCriticalType::Minimum);
        assert_eq!(crits[3].critical_type, MorseCriticalType::Maximum);
    }
    #[test]
    fn test_descending_manifold_goes_downhill() {
        let m = make_tet_mesh();
        let values = [0.0, 1.0, 2.0, 3.0];
        let path = trace_descending_manifold(&m, &values, 3, 10);
        assert!(!path.is_empty());
        assert_eq!(path[0], 3);
    }
    #[test]
    fn test_curve_skeleton_creation() {
        let m = make_tet_mesh();
        let skel = extract_curve_skeleton(&m, 5, 0.3);
        assert_eq!(skel.num_nodes(), 4);
    }
    #[test]
    fn test_curve_skeleton_total_length_positive() {
        let mut skel = CurveSkeleton::new();
        let a = skel.add_node([0.0; 3]);
        let b = skel.add_node([1.0, 0.0, 0.0]);
        skel.add_edge(a, b);
        assert!((skel.total_length() - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_medial_axis_approx_returns_candidates() {
        let verts: Vec<[f64; 3]> = (0..16)
            .map(|i| {
                let t = i as f64 / 16.0 * std::f64::consts::TAU;
                [t.cos(), t.sin(), 0.0]
            })
            .collect();
        let ma = compute_medial_axis_approx(&verts, 4);
        assert!(!ma.is_empty());
    }
    #[test]
    fn test_no_non_manifold_edges_in_tet() {
        let m = make_tet_mesh();
        let nme = detect_non_manifold_edges(&m);
        assert!(
            nme.is_empty(),
            "tetrahedron should have no non-manifold edges"
        );
    }
    #[test]
    fn test_no_non_manifold_vertices_in_tet() {
        let m = make_tet_mesh();
        let nmv = detect_non_manifold_vertices(&m);
        assert!(
            nmv.is_empty(),
            "tetrahedron should have no non-manifold vertices"
        );
    }
    #[test]
    fn test_boundary_loops_triangle_has_one_loop() {
        let m = make_triangle_mesh();
        let loops = extract_boundary_loops(&m);
        assert_eq!(
            loops.len(),
            1,
            "single triangle should have 1 boundary loop"
        );
        assert_eq!(loops[0].len(), 3, "boundary loop should have 3 vertices");
    }
    #[test]
    fn test_boundary_loops_tet_is_empty() {
        let m = make_tet_mesh();
        let loops = extract_boundary_loops(&m);
        assert!(loops.is_empty(), "closed tet should have no boundary loops");
    }
    #[test]
    fn test_connected_components_single() {
        let edges = [(0, 1), (1, 2), (2, 3)];
        let n = count_connected_components(4, &edges);
        assert_eq!(n, 1);
    }
    #[test]
    fn test_connected_components_two_islands() {
        let edges = [(0, 1), (2, 3)];
        let n = count_connected_components(4, &edges);
        assert_eq!(n, 2);
    }
    #[test]
    fn test_connected_components_isolated() {
        let edges: [(usize, usize); 0] = [];
        let n = count_connected_components(4, &edges);
        assert_eq!(n, 4);
    }
    #[test]
    fn test_mesh_topology_euler_characteristic_sphere() {
        let t = MeshTopology::new(4, 6, 4);
        assert_eq!(t.euler_characteristic(), 2);
    }
    #[test]
    fn test_mesh_topology_genus_sphere() {
        let t = MeshTopology::new(4, 6, 4);
        assert_eq!(t.genus(), 0);
    }
    #[test]
    fn test_mesh_topology_genus_torus() {
        let t = MeshTopology::new(9, 18, 9);
        assert_eq!(t.euler_characteristic(), 0);
        assert_eq!(t.genus(), 1);
    }
    #[test]
    fn test_mesh_topology_default() {
        let t = MeshTopology::default();
        assert_eq!(t.euler_characteristic(), 0);
    }
    #[test]
    fn test_mesh_topology_cycle_rank_tree() {
        let t = MeshTopology::new(5, 4, 0);
        assert_eq!(t.cycle_rank(), 0);
    }
    #[test]
    fn test_mesh_topology_cycle_rank_triangle() {
        let t = MeshTopology::new(3, 3, 1);
        assert_eq!(t.cycle_rank(), 1);
    }
    #[test]
    fn test_mesh_topology_betti_0_estimate() {
        let t = MeshTopology::new(4, 6, 4);
        assert_eq!(t.betti_0_estimate(), 1);
    }
    #[test]
    fn test_betti_0_connected_graph() {
        let edges = [(0, 1), (1, 2), (2, 3)];
        assert_eq!(betti_0(&edges, 4), 1);
    }
    #[test]
    fn test_betti_0_two_components() {
        let edges = [(0, 1), (2, 3)];
        assert_eq!(betti_0(&edges, 4), 2);
    }
    #[test]
    fn test_betti_0_isolated_vertices() {
        let edges: [(usize, usize); 0] = [];
        assert_eq!(betti_0(&edges, 5), 5);
    }
    #[test]
    fn test_betti_0_single_vertex() {
        let edges: [(usize, usize); 0] = [];
        assert_eq!(betti_0(&edges, 1), 1);
    }
    #[test]
    fn test_betti_0_zero_vertices() {
        let edges: [(usize, usize); 0] = [];
        assert_eq!(betti_0(&edges, 0), 0);
    }
    #[test]
    fn test_betti_1_estimate_sphere() {
        assert_eq!(betti_1_estimate(4, 6, 4), 0);
    }
    #[test]
    fn test_betti_1_estimate_torus() {
        assert_eq!(betti_1_estimate(9, 18, 9), 2);
    }
    #[test]
    fn test_betti_1_estimate_double_torus() {
        assert_eq!(betti_1_estimate(4, 10, 4), 4);
    }
    #[test]
    fn test_quadric_error_matrix_zero_faces() {
        let q = MeshSimplification::quadric_error_matrix([0.0; 3], &[]);
        for row in q.iter() {
            for &v in row.iter() {
                assert_eq!(v, 0.0);
            }
        }
    }
    #[test]
    fn test_quadric_error_matrix_one_face() {
        let face = [[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let q = MeshSimplification::quadric_error_matrix([0.0; 3], &[face]);
        let sum: f64 = q.iter().flat_map(|r| r.iter()).sum();
        assert!(sum.abs() > 0.0);
    }
    #[test]
    fn test_edge_collapse_cost_identical_vertices() {
        let q = [[0.0; 4]; 4];
        let cost = MeshSimplification::edge_collapse_cost([0.0; 3], [0.0; 3], q, q);
        assert_eq!(cost, 0.0);
    }
    #[test]
    fn test_simplify_qem_reduces_faces() {
        let mut positions = vec![
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ];
        let mut triangles = vec![
            [0usize, 1, 4],
            [1, 3, 4],
            [3, 2, 4],
            [2, 0, 4],
            [0, 1, 3],
            [0, 3, 2],
        ];
        let initial_count = triangles.len();
        MeshSimplification::simplify_qem(&mut positions, &mut triangles, 2);
        assert!(triangles.len() < initial_count);
        assert!(triangles.len() <= 2 || triangles.len() < initial_count);
    }
    #[test]
    fn test_simplify_qem_target_already_met() {
        let mut positions = vec![[0.0f64; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut triangles = vec![[0usize, 1, 2]];
        MeshSimplification::simplify_qem(&mut positions, &mut triangles, 5);
        assert_eq!(triangles.len(), 1);
    }
    #[test]
    fn test_simplify_qem_empty_mesh() {
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut triangles: Vec<[usize; 3]> = Vec::new();
        MeshSimplification::simplify_qem(&mut positions, &mut triangles, 0);
        assert!(triangles.is_empty());
    }
    #[test]
    fn test_marching_squares_empty_field() {
        let field: Vec<Vec<f64>> = Vec::new();
        let segs = IsoCurve::marching_squares(&field, 0, 0, 1.0, 1.0, 0.5);
        assert!(segs.is_empty());
    }
    #[test]
    fn test_marching_squares_all_below() {
        let field = vec![vec![0.0f64; 4]; 4];
        let segs = IsoCurve::marching_squares(&field, 4, 4, 1.0, 1.0, 1.0);
        assert!(segs.is_empty());
    }
    #[test]
    fn test_marching_squares_all_above() {
        let field = vec![vec![2.0f64; 4]; 4];
        let segs = IsoCurve::marching_squares(&field, 4, 4, 1.0, 1.0, 1.0);
        assert!(segs.is_empty());
    }
    #[test]
    fn test_marching_squares_single_corner() {
        let field = vec![vec![2.0, 0.0], vec![0.0, 0.0]];
        let segs = IsoCurve::marching_squares(&field, 2, 2, 1.0, 1.0, 1.0);
        assert_eq!(segs.len(), 1);
    }
    #[test]
    fn test_marching_squares_horizontal_contour() {
        let field = vec![vec![0.0f64, 0.0, 0.0, 0.0], vec![2.0, 2.0, 2.0, 2.0]];
        let segs = IsoCurve::marching_squares(&field, 4, 2, 1.0, 1.0, 1.0);
        assert!(!segs.is_empty());
    }
    #[test]
    fn test_marching_squares_returns_line_segments() {
        let field = vec![
            vec![0.0f64, 0.0, 2.0],
            vec![0.0, 2.0, 2.0],
            vec![2.0, 2.0, 2.0],
        ];
        let segs = IsoCurve::marching_squares(&field, 3, 3, 1.0, 1.0, 1.0);
        for seg in &segs {
            assert_eq!(seg.len(), 2);
        }
    }
    #[test]
    fn test_marching_squares_circle_like() {
        let n = 10usize;
        let mut field = vec![vec![0.0f64; n]; n];
        for (y, row) in field.iter_mut().enumerate() {
            for (x, cell) in row.iter_mut().enumerate() {
                let cx = 4.5f64;
                let cy = 4.5f64;
                let r2 = (x as f64 - cx).powi(2) + (y as f64 - cy).powi(2);
                *cell = (-r2 / 4.0).exp();
            }
        }
        let segs = IsoCurve::marching_squares(&field, n, n, 1.0, 1.0, 0.3);
        assert!(!segs.is_empty());
    }
    #[test]
    fn test_polygon_offset_preserves_vertex_count() {
        let square = vec![[0.0f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let offset = PolygonOffset::offset_polygon_2d(&square, 0.1);
        assert_eq!(offset.len(), 4);
    }
    #[test]
    fn test_polygon_offset_too_small_input() {
        let pts = vec![[0.0f64, 0.0], [1.0, 0.0]];
        let offset = PolygonOffset::offset_polygon_2d(&pts, 0.1);
        assert_eq!(offset.len(), 2);
    }
    #[test]
    fn test_polygon_offset_empty_input() {
        let pts: Vec<[f64; 2]> = Vec::new();
        let offset = PolygonOffset::offset_polygon_2d(&pts, 0.1);
        assert!(offset.is_empty());
    }
    #[test]
    fn test_polygon_offset_outward_increases_size() {
        let square = vec![[0.0f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let area_orig = PolygonOffset::signed_area(&square).abs();
        let offset = PolygonOffset::offset_polygon_2d(&square, 0.5);
        let area_offset = PolygonOffset::signed_area(&offset).abs();
        assert!(
            area_offset > area_orig,
            "area_orig={}, area_offset={}",
            area_orig,
            area_offset
        );
    }
    #[test]
    fn test_signed_area_ccw_square() {
        let square = vec![[0.0f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let area = PolygonOffset::signed_area(&square);
        assert!((area - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_signed_area_cw_square_negative() {
        let square = vec![[0.0f64, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]];
        let area = PolygonOffset::signed_area(&square);
        assert!(area < 0.0);
    }
    #[test]
    fn test_signed_area_triangle() {
        let tri = vec![[0.0f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let area = PolygonOffset::signed_area(&tri);
        assert!((area - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_is_ccw_square() {
        let square = vec![[0.0f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        assert!(PolygonOffset::is_ccw(&square));
    }
    #[test]
    fn test_is_ccw_cw_polygon() {
        let cw = vec![[0.0f64, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]];
        assert!(!PolygonOffset::is_ccw(&cw));
    }
    #[test]
    fn test_polygon_offset_zero_dist_unchanged() {
        let square = vec![[0.0f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let offset = PolygonOffset::offset_polygon_2d(&square, 0.0);
        for (a, b) in square.iter().zip(offset.iter()) {
            assert!((a[0] - b[0]).abs() < 1e-10);
            assert!((a[1] - b[1]).abs() < 1e-10);
        }
    }
}
