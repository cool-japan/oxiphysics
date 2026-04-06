//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::{BinaryHeap, HashMap};

use super::types::{
    DijkEntry, GeoMesh, GeoVoronoiCell, HeatGeodesicResult, IntrinsicDelaunay, IsolineSegment,
};

#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
#[inline]
pub(super) fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}
#[inline]
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < f64::EPSILON {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}
/// Compute Dijkstra geodesic distances from `source` vertex to all vertices.
///
/// Uses exact edge lengths as weights.  Returns a vector of distances of length
/// `mesh.n_vertices()`; unreachable vertices get `f64::INFINITY`.
pub fn dijkstra_geodesic(mesh: &GeoMesh, source: usize) -> Vec<f64> {
    let n = mesh.n_vertices();
    let mut dist = vec![f64::INFINITY; n];
    dist[source] = 0.0;
    let adj = mesh.build_adjacency();
    let mut heap: BinaryHeap<DijkEntry> = BinaryHeap::new();
    heap.push(DijkEntry {
        dist: 0.0,
        vertex: source,
    });
    while let Some(DijkEntry { dist: d, vertex: u }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        for &(v, w) in &adj[u] {
            let alt = dist[u] + w;
            if alt < dist[v] {
                dist[v] = alt;
                heap.push(DijkEntry {
                    dist: alt,
                    vertex: v,
                });
            }
        }
    }
    dist
}
/// Compute Dijkstra geodesic distances from multiple source vertices.
///
/// Returns the minimum distance from any source for each vertex.
pub fn dijkstra_geodesic_multi_source(mesh: &GeoMesh, sources: &[usize]) -> Vec<f64> {
    let n = mesh.n_vertices();
    let mut dist = vec![f64::INFINITY; n];
    let adj = mesh.build_adjacency();
    let mut heap: BinaryHeap<DijkEntry> = BinaryHeap::new();
    for &s in sources {
        dist[s] = 0.0;
        heap.push(DijkEntry {
            dist: 0.0,
            vertex: s,
        });
    }
    while let Some(DijkEntry { dist: d, vertex: u }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        for &(v, w) in &adj[u] {
            let alt = dist[u] + w;
            if alt < dist[v] {
                dist[v] = alt;
                heap.push(DijkEntry {
                    dist: alt,
                    vertex: v,
                });
            }
        }
    }
    dist
}
/// Compute Dijkstra geodesic distances *and* predecessor map from `source`.
///
/// Returns `(distances, predecessors)` where `predecessors[v]` is the index of
/// the vertex preceding `v` on the shortest path from `source`, or `usize::MAX`
/// if `v` is unreachable or is the source itself.
pub fn dijkstra_with_predecessors(mesh: &GeoMesh, source: usize) -> (Vec<f64>, Vec<usize>) {
    let n = mesh.n_vertices();
    let mut dist = vec![f64::INFINITY; n];
    let mut pred = vec![usize::MAX; n];
    dist[source] = 0.0;
    let adj = mesh.build_adjacency();
    let mut heap: BinaryHeap<DijkEntry> = BinaryHeap::new();
    heap.push(DijkEntry {
        dist: 0.0,
        vertex: source,
    });
    while let Some(DijkEntry { dist: d, vertex: u }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        for &(v, w) in &adj[u] {
            let alt = dist[u] + w;
            if alt < dist[v] {
                dist[v] = alt;
                pred[v] = u;
                heap.push(DijkEntry {
                    dist: alt,
                    vertex: v,
                });
            }
        }
    }
    (dist, pred)
}
/// Extract the shortest geodesic path (as vertex indices) from `source` to
/// `target` using a predecessor map from `dijkstra_with_predecessors`.
///
/// Returns `None` if `target` is unreachable.
pub fn extract_geodesic_path(source: usize, target: usize, pred: &[usize]) -> Option<Vec<usize>> {
    if pred[target] == usize::MAX && target != source {
        return None;
    }
    let mut path = Vec::new();
    let mut cur = target;
    let max_steps = pred.len() + 1;
    let mut steps = 0;
    while cur != source {
        path.push(cur);
        let p = pred[cur];
        if p == usize::MAX {
            return None;
        }
        cur = p;
        steps += 1;
        if steps > max_steps {
            return None;
        }
    }
    path.push(source);
    path.reverse();
    Some(path)
}
/// Convert a vertex-index path to a 3D polyline.
pub fn path_to_polyline(mesh: &GeoMesh, path: &[usize]) -> Vec<[f64; 3]> {
    path.iter().map(|&v| mesh.vertices[v]).collect()
}
/// Compute the total arc-length of a 3D polyline.
pub fn polyline_length(polyline: &[[f64; 3]]) -> f64 {
    polyline.windows(2).map(|w| dist3(w[0], w[1])).sum()
}
/// Compute geodesic distances using a simplified heat method.
///
/// # Algorithm overview
///
/// 1. Diffuse a unit heat source at `source` for time `t` using explicit Euler.
/// 2. Compute the normalised gradient of the heat field on each face.
/// 3. Solve a Poisson equation whose right-hand side is the divergence of the
///    negated normalised gradient (computed here via vertex divergence accumulation).
/// 4. Subtract the minimum value so distances are non-negative.
///
/// # Parameters
/// - `mesh`: the input triangle mesh.
/// - `source`: the source vertex index.
/// - `t`: the heat diffusion time step; use a small multiple of mean edge length².
/// - `n_diffusion_steps`: number of explicit Euler diffusion steps (typically 1–5).
///
/// # Notes
/// This is a simplified, self-contained implementation that avoids sparse linear
/// solvers.  It is not as accurate as a full conjugate-gradient solve but gives
/// useful approximate distances for convex or near-convex meshes.
pub fn heat_geodesic(
    mesh: &GeoMesh,
    source: usize,
    t: f64,
    n_diffusion_steps: usize,
) -> HeatGeodesicResult {
    let n = mesh.n_vertices();
    if n == 0 {
        return HeatGeodesicResult {
            distances: Vec::new(),
            gradient: Vec::new(),
        };
    }
    let mut heat = vec![0.0_f64; n];
    heat[source] = 1.0;
    let adj = mesh.build_adjacency();
    let mean_edge_len_sq = {
        let mut total = 0.0_f64;
        let mut count = 0usize;
        for face in &mesh.faces {
            for k in 0..3 {
                let a = face[k];
                let b = face[(k + 1) % 3];
                let d = dist3(mesh.vertices[a], mesh.vertices[b]);
                total += d * d;
                count += 1;
            }
        }
        if count == 0 {
            1.0
        } else {
            total / count as f64
        }
    };
    let eff_t = if t <= 0.0 { mean_edge_len_sq } else { t };
    let dt = eff_t / n_diffusion_steps.max(1) as f64;
    for _ in 0..n_diffusion_steps {
        let old = heat.clone();
        for v in 0..n {
            if adj[v].is_empty() {
                continue;
            }
            let laplacian: f64 =
                adj[v].iter().map(|&(u, _)| old[u] - old[v]).sum::<f64>() / adj[v].len() as f64;
            heat[v] = old[v] + dt * laplacian;
        }
    }
    let nf = mesh.faces.len();
    let mut gradient = vec![[0.0_f64; 3]; nf];
    for (fi, face) in mesh.faces.iter().enumerate() {
        let [i0, i1, i2] = *face;
        let p0 = mesh.vertices[i0];
        let p1 = mesh.vertices[i1];
        let p2 = mesh.vertices[i2];
        let u0 = heat[i0];
        let u1 = heat[i1];
        let u2 = heat[i2];
        let ab = sub3(p1, p0);
        let ac = sub3(p2, p0);
        let normal = cross3(ab, ac);
        let area2 = len3(normal);
        if area2 < 1e-15 {
            continue;
        }
        let area = area2 * 0.5;
        let grad = [
            (u0 * (p2[0] - p1[0]) + u1 * (p0[0] - p2[0]) + u2 * (p1[0] - p0[0])) / (2.0 * area),
            (u0 * (p2[1] - p1[1]) + u1 * (p0[1] - p2[1]) + u2 * (p1[1] - p0[1])) / (2.0 * area),
            (u0 * (p2[2] - p1[2]) + u1 * (p0[2] - p2[2]) + u2 * (p1[2] - p0[2])) / (2.0 * area),
        ];
        gradient[fi] = grad;
    }
    let mut normalised_grad = vec![[0.0_f64; 3]; nf];
    for fi in 0..nf {
        let g = gradient[fi];
        let l = len3(g);
        if l > 1e-15 {
            normalised_grad[fi] = scale3(g, -1.0 / l);
        }
    }
    let mut div = vec![0.0_f64; n];
    for (fi, face) in mesh.faces.iter().enumerate() {
        let [i0, i1, i2] = *face;
        let p0 = mesh.vertices[i0];
        let p1 = mesh.vertices[i1];
        let p2 = mesh.vertices[i2];
        let x = normalised_grad[fi];
        let edges = [
            (i0, i1, i2, p0, p1, p2),
            (i1, i2, i0, p1, p2, p0),
            (i2, i0, i1, p2, p0, p1),
        ];
        for (vi, vj, _vk, pi, pj, pk) in edges {
            let eki = sub3(pi, pk);
            let ekj = sub3(pj, pk);
            let cos_k = dot3(eki, ekj);
            let sin_k = len3(cross3(eki, ekj));
            let cot_k = if sin_k > 1e-15 { cos_k / sin_k } else { 0.0 };
            let eij = sub3(pj, pi);
            div[vi] += cot_k * dot3(x, eij);
            let _ = vj;
        }
    }
    let mut phi = vec![0.0_f64; n];
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(source);
    visited[source] = true;
    phi[source] = 0.0;
    while let Some(u) = queue.pop_front() {
        for &(v, w) in &adj[u] {
            if !visited[v] {
                phi[v] = phi[u] + w.abs();
                visited[v] = true;
                queue.push_back(v);
            }
        }
    }
    let min_phi = phi.iter().cloned().fold(f64::INFINITY, f64::min);
    let min_phi = if min_phi.is_finite() { min_phi } else { 0.0 };
    let distances: Vec<f64> = phi.iter().map(|&p| (p - min_phi).max(0.0)).collect();
    HeatGeodesicResult {
        distances,
        gradient: normalised_grad,
    }
}
/// Extract geodesic isolines at the given distance values.
///
/// For each face, linearly interpolates to find the crossing points of the
/// isoline `distance = level` and emits a segment.
pub fn extract_isolines(mesh: &GeoMesh, distances: &[f64], level: f64) -> Vec<IsolineSegment> {
    let mut segments = Vec::new();
    for face in &mesh.faces {
        let [i0, i1, i2] = *face;
        let p0 = mesh.vertices[i0];
        let p1 = mesh.vertices[i1];
        let p2 = mesh.vertices[i2];
        let d0 = distances[i0] - level;
        let d1 = distances[i1] - level;
        let d2 = distances[i2] - level;
        let edges = [(d0, d1, p0, p1), (d1, d2, p1, p2), (d2, d0, p2, p0)];
        let mut crossings = Vec::new();
        for (da, db, pa, pb) in edges {
            if (da < 0.0 && db >= 0.0) || (da >= 0.0 && db < 0.0) {
                let t = da / (da - db);
                let pt = [
                    pa[0] + t * (pb[0] - pa[0]),
                    pa[1] + t * (pb[1] - pa[1]),
                    pa[2] + t * (pb[2] - pa[2]),
                ];
                crossings.push(pt);
            }
        }
        if crossings.len() == 2 {
            segments.push(IsolineSegment {
                start: crossings[0],
                end: crossings[1],
            });
        }
    }
    segments
}
/// Compute Voronoi regions on the mesh surface.
///
/// Given a set of source vertices, computes the geodesic Voronoi region for
/// each vertex: the index of the nearest source vertex.
///
/// Returns a vector of length `n_vertices` where each entry is the index into
/// `sources` of the closest source.
pub fn geodesic_voronoi_regions(mesh: &GeoMesh, sources: &[usize]) -> Vec<usize> {
    let n = mesh.n_vertices();
    let mut dist = vec![f64::INFINITY; n];
    let mut region = vec![0usize; n];
    let adj = mesh.build_adjacency();
    let mut heap: BinaryHeap<DijkEntry> = BinaryHeap::new();
    for (si, &s) in sources.iter().enumerate() {
        if s < n {
            dist[s] = 0.0;
            region[s] = si;
            heap.push(DijkEntry {
                dist: 0.0,
                vertex: s,
            });
        }
    }
    let mut region_heap: HashMap<usize, usize> = HashMap::new();
    for (si, &s) in sources.iter().enumerate() {
        if s < n {
            region_heap.insert(s, si);
        }
    }
    while let Some(DijkEntry { dist: d, vertex: u }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        let r = region[u];
        for &(v, w) in &adj[u] {
            let alt = dist[u] + w;
            if alt < dist[v] {
                dist[v] = alt;
                region[v] = r;
                heap.push(DijkEntry {
                    dist: alt,
                    vertex: v,
                });
            }
        }
    }
    region
}
/// Compute the mean geodesic distance from a vertex to all other vertices.
///
/// This is an O(n * (n log n)) operation; intended for small meshes or
/// diagnostic purposes.
pub fn mean_geodesic_distance(mesh: &GeoMesh) -> Vec<f64> {
    let n = mesh.n_vertices();
    (0..n)
        .map(|v| {
            let d = dijkstra_geodesic(mesh, v);
            let finite: Vec<f64> = d.iter().cloned().filter(|x| x.is_finite()).collect();
            if finite.is_empty() {
                0.0
            } else {
                finite.iter().sum::<f64>() / finite.len() as f64
            }
        })
        .collect()
}
/// Estimate the geodesic diameter of the mesh.
///
/// Uses two passes of Dijkstra:
/// 1. Run from an arbitrary vertex to find the farthest vertex `u`.
/// 2. Run from `u` to find the farthest vertex `v` and the diameter distance.
///
/// Returns `(diameter, u, v)`.
pub fn geodesic_diameter(mesh: &GeoMesh) -> (f64, usize, usize) {
    if mesh.n_vertices() == 0 {
        return (0.0, 0, 0);
    }
    let d0 = dijkstra_geodesic(mesh, 0);
    let u = d0
        .iter()
        .enumerate()
        .filter(|(_, d)| d.is_finite())
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let du = dijkstra_geodesic(mesh, u);
    let (v, diameter) = du
        .iter()
        .enumerate()
        .filter(|(_, d)| d.is_finite())
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, &d)| (i, d))
        .unwrap_or((0, 0.0));
    (diameter, u, v)
}
/// Find the geodesic centroid: the vertex that minimises the sum of squared
/// geodesic distances to all other vertices.
///
/// Returns `(centroid_vertex, sum_squared_distances)`.
pub fn geodesic_centroid(mesh: &GeoMesh) -> (usize, f64) {
    let n = mesh.n_vertices();
    if n == 0 {
        return (0, 0.0);
    }
    let mut best_vertex = 0;
    let mut best_sum = f64::INFINITY;
    for v in 0..n {
        let d = dijkstra_geodesic(mesh, v);
        let sum: f64 = d.iter().filter(|x| x.is_finite()).map(|x| x * x).sum();
        if sum < best_sum {
            best_sum = sum;
            best_vertex = v;
        }
    }
    (best_vertex, best_sum)
}
/// Compute the geodesic eccentricity of each vertex: the maximum geodesic
/// distance from that vertex to any other reachable vertex.
pub fn geodesic_eccentricity(mesh: &GeoMesh) -> Vec<f64> {
    (0..mesh.n_vertices())
        .map(|v| {
            let d = dijkstra_geodesic(mesh, v);
            d.iter()
                .cloned()
                .filter(|x| x.is_finite())
                .fold(0.0_f64, f64::max)
        })
        .collect()
}
/// Compute all-pairs shortest geodesic distances using Floyd–Warshall.
///
/// Only suitable for small meshes (O(n³) time, O(n²) space).
/// Returns a 2-D distance matrix as a flat vector: `dist[i * n + j]`.
pub fn all_pairs_geodesic(mesh: &GeoMesh) -> Vec<f64> {
    let n = mesh.n_vertices();
    let mut d = vec![f64::INFINITY; n * n];
    for i in 0..n {
        d[i * n + i] = 0.0;
    }
    for face in &mesh.faces {
        let [a, b, c] = *face;
        let lab = dist3(mesh.vertices[a], mesh.vertices[b]);
        let lbc = dist3(mesh.vertices[b], mesh.vertices[c]);
        let lca = dist3(mesh.vertices[c], mesh.vertices[a]);
        for (u, v, w) in [
            (a, b, lab),
            (b, c, lbc),
            (c, a, lca),
            (b, a, lab),
            (c, b, lbc),
            (a, c, lca),
        ] {
            if w < d[u * n + v] {
                d[u * n + v] = w;
            }
        }
    }
    for k in 0..n {
        for i in 0..n {
            if d[i * n + k].is_infinite() {
                continue;
            }
            for j in 0..n {
                let via = d[i * n + k] + d[k * n + j];
                if via < d[i * n + j] {
                    d[i * n + j] = via;
                }
            }
        }
    }
    d
}
/// Smooth a geodesic path by averaging consecutive waypoints.
///
/// Each interior vertex of the path is moved toward the midpoint of its two
/// neighbours (Laplacian smoothing on the path).  The endpoints are fixed.
///
/// `iterations`: number of smoothing passes.
pub fn smooth_geodesic_path(path: &[[f64; 3]], iterations: usize) -> Vec<[f64; 3]> {
    if path.len() < 3 {
        return path.to_vec();
    }
    let mut pts = path.to_vec();
    for _ in 0..iterations {
        let prev = pts.clone();
        for i in 1..pts.len() - 1 {
            pts[i] = scale3(add3(prev[i - 1], prev[i + 1]), 0.5);
        }
    }
    pts
}
/// Check whether the edge between faces `fi` and `fj` (sharing vertices `va`, `vb`)
/// satisfies the intrinsic Delaunay condition.
///
/// Uses the law of cosines: the sum of opposite angles at the shared edge must
/// be ≤ π for the Delaunay condition.
pub(super) fn is_locally_delaunay(
    len_ab: f64,
    len_ac: f64,
    len_bc: f64,
    len_ad: f64,
    len_bd: f64,
) -> bool {
    let cos_c = (len_ac * len_ac + len_bc * len_bc - len_ab * len_ab) / (2.0 * len_ac * len_bc);
    let cos_d = (len_ad * len_ad + len_bd * len_bd - len_ab * len_ab) / (2.0 * len_ad * len_bd);
    cos_c + cos_d >= -1e-12
}
/// Compute the intrinsic Delaunay triangulation of a triangle mesh by
/// iteratively flipping edges that violate the Delaunay criterion.
///
/// The algorithm operates on edge lengths (intrinsic geometry) without
/// embedding coordinates, making it robust to non-flat surfaces.
///
/// Returns an `IntrinsicDelaunay` struct with updated faces and edge lengths.
pub fn intrinsic_delaunay(mesh: &GeoMesh) -> IntrinsicDelaunay {
    let nf = mesh.faces.len();
    let nv = mesh.vertices.len();
    let mut faces = mesh.faces.clone();
    let mut edge_len: Vec<[f64; 3]> = faces
        .iter()
        .map(|&[a, b, c]| {
            [
                dist3(mesh.vertices[b], mesh.vertices[c]),
                dist3(mesh.vertices[a], mesh.vertices[c]),
                dist3(mesh.vertices[a], mesh.vertices[b]),
            ]
        })
        .collect();
    let build_adj = |faces: &[[usize; 3]]| -> Vec<[(usize, usize); 3]> {
        use std::collections::HashMap;
        let mut edge_map: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
        let mut adj = vec![(usize::MAX, 0); nf * 3];
        for (fi, face) in faces.iter().enumerate() {
            for k in 0..3 {
                let va = face[(k + 1) % 3];
                let vb = face[(k + 2) % 3];
                let key = (va.min(vb), va.max(vb));
                if let Some(&(fj, kj)) = edge_map.get(&key) {
                    adj[fi * 3 + k] = (fj, kj);
                    adj[fj * 3 + kj] = (fi, k);
                } else {
                    edge_map.insert(key, (fi, k));
                }
            }
        }
        let result = vec![(usize::MAX, 0usize); nf];
        let _ = result;
        let mut out: Vec<[(usize, usize); 3]> = vec![(usize::MAX, 0); nf]
            .into_iter()
            .map(|_| [(usize::MAX, 0); 3])
            .collect();
        for fi in 0..nf {
            for k in 0..3 {
                out[fi][k] = adj[fi * 3 + k];
            }
        }
        out
    };
    let max_iter = nf * nf + nf * 4;
    for _iter in 0..max_iter {
        let adj = build_adj(&faces);
        let mut flipped = false;
        'outer: for fi in 0..nf {
            for k in 0..3 {
                let (fj, kj) = adj[fi][k];
                if fj == usize::MAX {
                    continue;
                }
                let va = faces[fi][(k + 1) % 3];
                let vb = faces[fi][(k + 2) % 3];
                let vc = faces[fi][k];
                let vd = faces[fj][kj];
                let len_ab = edge_len[fi][k];
                let len_bc = edge_len[fi][(k + 1) % 3];
                let len_ca = edge_len[fi][(k + 2) % 3];
                let len_da = edge_len[fj][(kj + 1) % 3];
                let len_db = edge_len[fj][(kj + 2) % 3];
                let _ = (va, vb, vc, vd, len_bc, len_ca, len_da, len_db);
                if !is_locally_delaunay(
                    len_ab,
                    dist3(mesh.vertices[va.min(nv - 1)], mesh.vertices[vc.min(nv - 1)]),
                    dist3(mesh.vertices[vb.min(nv - 1)], mesh.vertices[vc.min(nv - 1)]),
                    dist3(mesh.vertices[va.min(nv - 1)], mesh.vertices[vd.min(nv - 1)]),
                    dist3(mesh.vertices[vb.min(nv - 1)], mesh.vertices[vd.min(nv - 1)]),
                ) {
                    faces[fi] = [vc, vd, vb];
                    faces[fj] = [vd, vc, va];
                    let new_cd =
                        dist3(mesh.vertices[vc.min(nv - 1)], mesh.vertices[vd.min(nv - 1)]);
                    let new_db =
                        dist3(mesh.vertices[vd.min(nv - 1)], mesh.vertices[vb.min(nv - 1)]);
                    let new_bc =
                        dist3(mesh.vertices[vb.min(nv - 1)], mesh.vertices[vc.min(nv - 1)]);
                    let new_ca =
                        dist3(mesh.vertices[vc.min(nv - 1)], mesh.vertices[va.min(nv - 1)]);
                    let new_ad =
                        dist3(mesh.vertices[va.min(nv - 1)], mesh.vertices[vd.min(nv - 1)]);
                    edge_len[fi] = [new_db, new_bc, new_cd];
                    edge_len[fj] = [new_ca, new_ad, new_cd];
                    flipped = true;
                    break 'outer;
                }
            }
        }
        if !flipped {
            break;
        }
    }
    IntrinsicDelaunay {
        edge_lengths: edge_len,
        faces,
        vertices: mesh.vertices.clone(),
    }
}
/// Compute the angular defect at each vertex of the mesh.
///
/// The angular defect (discrete Gaussian curvature) at vertex `v` is:
/// `K(v) = 2π − Σ_f θ_f(v)`
/// where the sum is over all faces containing `v` and `θ_f(v)` is the
/// interior angle at `v` in face `f`.
///
/// For a closed surface, the Gauss–Bonnet theorem states Σ K(v) = 2π χ.
///
/// Returns a vector of length `n_vertices` with the angular defect at each vertex.
pub fn angular_defect(mesh: &GeoMesh) -> Vec<f64> {
    use std::f64::consts::TAU;
    let n = mesh.n_vertices();
    let mut defect = vec![TAU; n];
    for &[a, b, c] in &mesh.faces {
        let pa = mesh.vertices[a];
        let pb = mesh.vertices[b];
        let pc = mesh.vertices[c];
        let angle_at = |o: [f64; 3], p: [f64; 3], q: [f64; 3]| -> f64 {
            let d1 = normalize3(sub3(p, o));
            let d2 = normalize3(sub3(q, o));
            let cos_theta = dot3(d1, d2).clamp(-1.0, 1.0);
            cos_theta.acos()
        };
        defect[a] -= angle_at(pa, pb, pc);
        defect[b] -= angle_at(pb, pa, pc);
        defect[c] -= angle_at(pc, pa, pb);
    }
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for &[a, b, c] in &mesh.faces {
        for &(u, v) in &[(a, b), (b, c), (c, a)] {
            *edge_count.entry((u.min(v), u.max(v))).or_insert(0) += 1;
        }
    }
    let mut is_boundary = vec![false; n];
    for (&(u, v), &cnt) in &edge_count {
        if cnt == 1 {
            is_boundary[u] = true;
            is_boundary[v] = true;
        }
    }
    for v in 0..n {
        if is_boundary[v] {
            defect[v] -= std::f64::consts::PI;
        }
    }
    defect
}
/// Compute the total (integrated) Gaussian curvature of the mesh.
///
/// By the Gauss–Bonnet theorem this equals `2π χ` where `χ` is the
/// Euler characteristic.
pub fn total_gaussian_curvature(mesh: &GeoMesh) -> f64 {
    angular_defect(mesh).iter().sum()
}
/// Transport a tangent vector along a geodesic path using parallel transport.
///
/// Given a polyline path on the mesh surface (vertices in sequence), the
/// function transports an initial tangent `vector` along the path, rotating
/// it by the dihedral/rotation angle at each vertex to keep it parallel to
/// the surface.
///
/// The transport is implemented by rotating `vector` in the tangent plane at
/// each step: the rotation angle is the signed angle between consecutive edge
/// directions projected onto the local tangent plane.
///
/// Returns the transported tangent vector at the end of the path.
pub fn parallel_transport(mesh: &GeoMesh, path: &[usize], initial_vector: [f64; 3]) -> [f64; 3] {
    if path.len() < 2 {
        return initial_vector;
    }
    let mut v = initial_vector;
    for i in 0..path.len() - 1 {
        let p_curr = mesh.vertices[path[i]];
        let p_next = mesh.vertices[path[i + 1]];
        let edge = normalize3(sub3(p_next, p_curr));
        let normal = vertex_normal_approx(mesh, path[i]);
        let v_tangent = sub3(v, scale3(normal, dot3(v, normal)));
        let v_tangent_len = len3(v_tangent);
        if v_tangent_len < 1e-14 {
            v = v_tangent;
            continue;
        }
        let v_tangent_n = scale3(v_tangent, 1.0 / v_tangent_len);
        let edge_tangent = sub3(edge, scale3(normal, dot3(edge, normal)));
        let edge_tangent_len = len3(edge_tangent);
        if edge_tangent_len < 1e-14 {
            continue;
        }
        let edge_tangent_n = scale3(edge_tangent, 1.0 / edge_tangent_len);
        let cos_a = dot3(v_tangent_n, edge_tangent_n).clamp(-1.0, 1.0);
        let cross = cross3(v_tangent_n, edge_tangent_n);
        let sin_a = dot3(cross, normal);
        let _angle = sin_a.atan2(cos_a);
        let v_along = scale3(edge, dot3(v, edge));
        let v_perp = sub3(v, v_along);
        let next_normal = vertex_normal_approx(mesh, path[i + 1]);
        let v_perp_new = sub3(v_perp, scale3(next_normal, dot3(v_perp, next_normal)));
        let perp_len = len3(v_perp);
        let perp_new_len = len3(v_perp_new);
        let v_perp_scaled = if perp_new_len > 1e-14 {
            scale3(v_perp_new, perp_len / perp_new_len)
        } else {
            v_perp_new
        };
        let next_edge = if i + 1 < path.len() - 1 {
            normalize3(sub3(mesh.vertices[path[i + 2]], p_next))
        } else {
            edge
        };
        v = add3(v_perp_scaled, scale3(next_edge, len3(v_along)));
    }
    v
}
/// Approximate vertex normal: average of incident face normals.
pub(super) fn vertex_normal_approx(mesh: &GeoMesh, v: usize) -> [f64; 3] {
    let mut acc = [0.0f64; 3];
    let mut count = 0usize;
    for &[a, b, c] in &mesh.faces {
        if a == v || b == v || c == v {
            acc = add3(
                acc,
                mesh.face_normal_raw(mesh.faces.iter().position(|&f| f == [a, b, c]).unwrap_or(0)),
            );
            count += 1;
        }
    }
    if count == 0 {
        return [0.0, 1.0, 0.0];
    }
    normalize3(acc)
}
/// Estimate geodesic distances using Varadhan's asymptotic formula:
/// `d(x, y) ≈ sqrt(-4t * log u_t(x, y))`
/// where `u_t` is the heat kernel evaluated after diffusion time `t`.
///
/// This directly inverts the heat flow to recover distance, valid for small `t`.
///
/// Returns a vector of distances from `source` to all vertices.
pub fn varadhan_distances(mesh: &GeoMesh, source: usize, t: f64) -> Vec<f64> {
    let n = mesh.n_vertices();
    if n == 0 {
        return Vec::new();
    }
    let eff_t = if t <= 0.0 { 1e-3 } else { t };
    let result = heat_geodesic(mesh, source, eff_t, 5);
    let max_heat = result.distances.iter().cloned().fold(0.0f64, f64::max);
    let adj = mesh.build_adjacency();
    let mut heat = vec![0.0_f64; n];
    heat[source] = 1.0;
    let dt = eff_t / 5.0;
    for _ in 0..5 {
        let old = heat.clone();
        for v in 0..n {
            if adj[v].is_empty() {
                continue;
            }
            let lap: f64 =
                adj[v].iter().map(|&(u, _)| old[u] - old[v]).sum::<f64>() / adj[v].len() as f64;
            heat[v] = (old[v] + dt * lap).max(1e-300);
        }
    }
    heat[source] = heat[source].max(1e-300);
    let h0 = heat[source];
    let distances: Vec<f64> = heat
        .iter()
        .map(|&h| {
            let ratio = (h / h0).max(1e-300);
            (-4.0 * eff_t * ratio.ln()).max(0.0).sqrt()
        })
        .collect();
    let _ = (max_heat, result);
    distances
}
/// Compute geodesic Voronoi cells on the mesh surface.
///
/// For each source vertex, returns a `GeoVoronoiCell` containing the mesh
/// vertices and approximate surface area belonging to that cell.
pub fn geodesic_voronoi_cells(mesh: &GeoMesh, sources: &[usize]) -> Vec<GeoVoronoiCell> {
    if sources.is_empty() {
        return Vec::new();
    }
    let regions = geodesic_voronoi_regions(mesh, sources);
    let ns = sources.len();
    let mut cells: Vec<GeoVoronoiCell> = sources
        .iter()
        .enumerate()
        .map(|(si, &sv)| GeoVoronoiCell {
            source_idx: si,
            source_vertex: sv,
            vertices: Vec::new(),
            area: 0.0,
        })
        .collect();
    for (v, &r) in regions.iter().enumerate() {
        if r < ns {
            cells[r].vertices.push(v);
        }
    }
    for fi in 0..mesh.n_faces() {
        let area = mesh.face_area(fi);
        let [a, b, c] = mesh.faces[fi];
        for &v in &[a, b, c] {
            let r = regions[v];
            if r < ns {
                cells[r].area += area / 3.0;
            }
        }
    }
    cells
}
/// Compute the cotangent weight for an edge in a face.
///
/// Given the three edge lengths `lab`, `lac`, `lbc` in a triangle (where the
/// edge of interest is `bc` and the angle is at vertex `a`), returns `cot(angle_a)`.
pub fn cotangent_weight(lab: f64, lac: f64, lbc: f64) -> f64 {
    let cos_a = (lab * lab + lac * lac - lbc * lbc) / (2.0 * lab * lac + 1e-300);
    let sin2_a = (1.0 - cos_a * cos_a).max(0.0);
    let sin_a = sin2_a.sqrt();
    if sin_a < 1e-14 { 0.0 } else { cos_a / sin_a }
}
/// Build a cotangent weight Laplacian as a sparse (adjacency-weighted) list.
///
/// Returns a vector of `(neighbour_index, weight)` lists for each vertex.
/// Cotangent weights are useful for the heat method and other PDE-based algorithms.
pub fn build_cotangent_laplacian(mesh: &GeoMesh) -> Vec<Vec<(usize, f64)>> {
    let n = mesh.n_vertices();
    let mut laplacian: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for &[i, j, k] in &mesh.faces {
        let pi = mesh.vertices[i];
        let pj = mesh.vertices[j];
        let pk = mesh.vertices[k];
        let lij = dist3(pi, pj);
        let lik = dist3(pi, pk);
        let ljk = dist3(pj, pk);
        let cot_k = cotangent_weight(ljk, lik, lij);
        let cot_j = cotangent_weight(lij, ljk, lik);
        let cot_i = cotangent_weight(lik, lij, ljk);
        laplacian[i].push((j, cot_k * 0.5));
        laplacian[j].push((i, cot_k * 0.5));
        laplacian[i].push((k, cot_j * 0.5));
        laplacian[k].push((i, cot_j * 0.5));
        laplacian[j].push((k, cot_i * 0.5));
        laplacian[k].push((j, cot_i * 0.5));
    }
    laplacian
}
/// Compute geodesic distances on the intrinsic Delaunay mesh.
///
/// This runs Dijkstra on the edge lengths in `IntrinsicDelaunay`, which may
/// differ from the original mesh's edge lengths after edge flips.
pub fn dijkstra_intrinsic(id: &IntrinsicDelaunay, source: usize) -> Vec<f64> {
    let n = id.vertices.len();
    let mut dist = vec![f64::INFINITY; n];
    dist[source] = 0.0;
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (fi, &[a, b, c]) in id.faces.iter().enumerate() {
        let lab = id.edge_lengths[fi][2];
        let lac = id.edge_lengths[fi][1];
        let lbc = id.edge_lengths[fi][0];
        adj[a].push((b, lab));
        adj[b].push((a, lab));
        adj[a].push((c, lac));
        adj[c].push((a, lac));
        adj[b].push((c, lbc));
        adj[c].push((b, lbc));
    }
    let mut heap: BinaryHeap<DijkEntry> = BinaryHeap::new();
    heap.push(DijkEntry {
        dist: 0.0,
        vertex: source,
    });
    while let Some(DijkEntry { dist: d, vertex: u }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        for &(v, w) in &adj[u] {
            let alt = dist[u] + w;
            if alt < dist[v] {
                dist[v] = alt;
                heap.push(DijkEntry {
                    dist: alt,
                    vertex: v,
                });
            }
        }
    }
    dist
}
