//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::collections::BinaryHeap;

use super::types::{
    DistNode, ExpMapResult, FmmStatus, GeodesicMesh, GeodesicVoronoiCell, HeatMethodParams,
    LogMapResult,
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
#[inline]
pub(super) fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}
/// Compute the angle (in radians) at vertex `b` in triangle `a-b-c`.
#[inline]
pub(super) fn angle_at_vertex(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = normalize3(sub3(a, b));
    let bc = normalize3(sub3(c, b));
    let d = dot3(ba, bc).clamp(-1.0, 1.0);
    d.acos()
}
/// Compute geodesic distances from a set of source vertices using Dijkstra's
/// algorithm on mesh edges.
///
/// Returns a vector of distances for each vertex. Unreachable vertices get `f64::INFINITY`.
pub fn dijkstra_geodesic(mesh: &GeodesicMesh, sources: &[usize]) -> Vec<f64> {
    let n = mesh.num_vertices();
    let adj = mesh.vertex_adjacency();
    let mut dist = vec![f64::INFINITY; n];
    let mut heap = BinaryHeap::new();
    for &s in sources {
        dist[s] = 0.0;
        heap.push(DistNode {
            vertex: s,
            dist: 0.0,
        });
    }
    while let Some(DistNode { vertex: u, dist: d }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        for &v in &adj[u] {
            let w = dist3(mesh.vertices[u], mesh.vertices[v]);
            let new_dist = dist[u] + w;
            if new_dist < dist[v] {
                dist[v] = new_dist;
                heap.push(DistNode {
                    vertex: v,
                    dist: new_dist,
                });
            }
        }
    }
    dist
}
/// Solve the 2D Eikonal equation on a triangle for vertex `c` given known
/// distances at vertices `a` and `b`. Returns the candidate distance at `c`.
///
/// Uses the acute-triangle update from the fast marching method on surfaces.
pub(super) fn eikonal_update_triangle(
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    da: f64,
    db: f64,
) -> f64 {
    let ab = dist3(pa, pb);
    let ac = dist3(pa, pc);
    let bc = dist3(pb, pc);
    if ab < 1e-15 {
        return da.min(db) + ac;
    }
    let u = db - da;
    let cos_a = ((ac * ac + ab * ab - bc * bc) / (2.0 * ac * ab)).clamp(-1.0, 1.0);
    let sin_a = (1.0 - cos_a * cos_a).max(0.0).sqrt();
    if sin_a < 1e-15 {
        return (da + ac).min(db + bc);
    }
    let t_candidate = (u + ab * cos_a) / (ab * sin_a);
    if t_candidate.abs() < 1e-15 {
        return da + ac;
    }
    let _a_coeff = ab * ab;
    let _b_coeff = -2.0 * ab * u;
    let _c_coeff = u * u - ac * ac * sin_a * sin_a;
    let discriminant = _b_coeff * _b_coeff - 4.0 * _a_coeff * _c_coeff;
    if discriminant < 0.0 {
        return (da + ac).min(db + bc);
    }
    let sqrt_disc = discriminant.sqrt();
    let t1 = (-_b_coeff + sqrt_disc) / (2.0 * _a_coeff);
    let t2 = (-_b_coeff - sqrt_disc) / (2.0 * _a_coeff);
    let mut best = f64::INFINITY;
    for t in [t1, t2] {
        if (-1e-10..=1.0 + 1e-10).contains(&t) {
            let t_clamped = t.clamp(0.0, 1.0);
            let pt = lerp3(pa, pb, t_clamped);
            let d_pt = da * (1.0 - t_clamped) + db * t_clamped;
            let candidate = d_pt + dist3(pt, pc);
            best = best.min(candidate);
        }
    }
    best = best.min(da + ac);
    best = best.min(db + bc);
    best
}
/// Compute geodesic distances from source vertices using the Fast Marching
/// Method (FMM) on triangle meshes.
///
/// This produces more accurate distances than pure Dijkstra by solving the
/// Eikonal equation on triangle faces rather than restricting to mesh edges.
///
/// Returns a distance vector for all vertices.
pub fn fast_marching_geodesic(mesh: &GeodesicMesh, sources: &[usize]) -> Vec<f64> {
    let n = mesh.num_vertices();
    let adj = mesh.vertex_adjacency();
    let v2f = mesh.vertex_to_faces();
    let mut dist = vec![f64::INFINITY; n];
    let mut status = vec![FmmStatus::Far; n];
    let mut heap = BinaryHeap::new();
    for &s in sources {
        dist[s] = 0.0;
        status[s] = FmmStatus::Known;
    }
    for &s in sources {
        for &nb in &adj[s] {
            if status[nb] == FmmStatus::Known {
                continue;
            }
            let d = dist3(mesh.vertices[s], mesh.vertices[nb]);
            if d < dist[nb] {
                dist[nb] = d;
            }
            if status[nb] == FmmStatus::Far {
                status[nb] = FmmStatus::Trial;
                heap.push(DistNode {
                    vertex: nb,
                    dist: dist[nb],
                });
            }
        }
    }
    while let Some(DistNode { vertex: u, dist: d }) = heap.pop() {
        if status[u] == FmmStatus::Known {
            continue;
        }
        if d > dist[u] + 1e-15 {
            continue;
        }
        status[u] = FmmStatus::Known;
        for &nb in &adj[u] {
            if status[nb] == FmmStatus::Known {
                continue;
            }
            let mut best_d = dist[nb];
            for &fi in &v2f[u] {
                let f = mesh.faces[fi];
                if !f.contains(&nb) {
                    continue;
                }
                let third = f
                    .iter()
                    .copied()
                    .find(|&v| v != u && v != nb)
                    .expect("element must exist");
                if status[third] == FmmStatus::Known {
                    let candidate = eikonal_update_triangle(
                        mesh.vertices[u],
                        mesh.vertices[third],
                        mesh.vertices[nb],
                        dist[u],
                        dist[third],
                    );
                    best_d = best_d.min(candidate);
                }
                let edge_d = dist[u] + dist3(mesh.vertices[u], mesh.vertices[nb]);
                best_d = best_d.min(edge_d);
            }
            if best_d < dist[nb] {
                dist[nb] = best_d;
                status[nb] = FmmStatus::Trial;
                heap.push(DistNode {
                    vertex: nb,
                    dist: best_d,
                });
            }
        }
    }
    dist
}
/// Compute mean edge length of the mesh.
pub fn mean_edge_length(mesh: &GeodesicMesh) -> f64 {
    let mut total = 0.0;
    let mut count = 0usize;
    for f in &mesh.faces {
        for i in 0..3 {
            let a = mesh.vertices[f[i]];
            let b = mesh.vertices[f[(i + 1) % 3]];
            total += dist3(a, b);
            count += 1;
        }
    }
    if count == 0 {
        1.0
    } else {
        total / count as f64
    }
}
/// Simplified Laplacian matrix (cotangent weights) for the mesh.
/// Returns a sparse representation as `Vec<Vec<(usize, f64)>>` (row -> list of (col, value)).
pub(super) fn build_cotan_laplacian(mesh: &GeodesicMesh) -> Vec<Vec<(usize, f64)>> {
    let n = mesh.num_vertices();
    let mut lap: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for f in &mesh.faces {
        for k in 0..3 {
            let i = f[k];
            let j = f[(k + 1) % 3];
            let o = f[(k + 2) % 3];
            let ei = sub3(mesh.vertices[i], mesh.vertices[o]);
            let ej = sub3(mesh.vertices[j], mesh.vertices[o]);
            let cos_o = dot3(ei, ej) / (len3(ei) * len3(ej) + 1e-30);
            let sin_o = len3(cross3(ei, ej)) / (len3(ei) * len3(ej) + 1e-30);
            let cot_o = if sin_o.abs() > 1e-15 {
                cos_o / sin_o
            } else {
                0.0
            };
            add_sparse_entry(&mut lap[i], j, 0.5 * cot_o);
            add_sparse_entry(&mut lap[j], i, 0.5 * cot_o);
            add_sparse_entry(&mut lap[i], i, -0.5 * cot_o);
            add_sparse_entry(&mut lap[j], j, -0.5 * cot_o);
        }
    }
    lap
}
/// Add a value to a sparse row entry.
pub(super) fn add_sparse_entry(row: &mut Vec<(usize, f64)>, col: usize, val: f64) {
    if let Some(entry) = row.iter_mut().find(|(c, _)| *c == col) {
        entry.1 += val;
    } else {
        row.push((col, val));
    }
}
/// Solve a sparse linear system `A * x = b` using Jacobi iteration.
/// `a` is the sparse matrix, `b` is the right-hand side, returns `x`.
pub(super) fn jacobi_solve(
    a: &[Vec<(usize, f64)>],
    b: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = b.len();
    let mut x = b.to_vec();
    let mut x_new = vec![0.0; n];
    let diag: Vec<f64> = (0..n)
        .map(|i| {
            a[i].iter()
                .find(|(c, _)| *c == i)
                .map(|(_, v)| *v)
                .unwrap_or(1.0)
        })
        .collect();
    for _iter in 0..max_iter {
        let mut max_diff = 0.0_f64;
        for i in 0..n {
            let mut sum = b[i];
            for &(j, v) in &a[i] {
                if j != i {
                    sum -= v * x[j];
                }
            }
            let d = diag[i];
            x_new[i] = if d.abs() > 1e-30 { sum / d } else { x[i] };
            max_diff = max_diff.max((x_new[i] - x[i]).abs());
        }
        std::mem::swap(&mut x, &mut x_new);
        if max_diff < tol {
            break;
        }
    }
    x
}
/// Compute geodesic distances from source vertices using the heat method.
///
/// Steps:
/// 1. Integrate heat from source (solve `(M - t*L) u = delta`).
/// 2. Compute normalized negative gradient of `u`.
/// 3. Solve Poisson equation `L phi = div(X)` for the distance.
///
/// This is a simplified implementation suitable for moderate-size meshes.
pub fn heat_method_geodesic(
    mesh: &GeodesicMesh,
    sources: &[usize],
    params: &HeatMethodParams,
) -> Vec<f64> {
    let n = mesh.num_vertices();
    if n == 0 || sources.is_empty() {
        return vec![0.0; n];
    }
    let h = mean_edge_length(mesh);
    let t = params.time_factor * h * h;
    let lap = build_cotan_laplacian(mesh);
    let mass: Vec<f64> = (0..n).map(|i| mesh.voronoi_area(i)).collect();
    let mut sys: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for i in 0..n {
        add_sparse_entry(&mut sys[i], i, mass[i]);
        for &(j, v) in &lap[i] {
            add_sparse_entry(&mut sys[i], j, -t * v);
        }
    }
    let mut rhs = vec![0.0; n];
    for &s in sources {
        rhs[s] = 1.0;
    }
    let u = jacobi_solve(&sys, &rhs, 500, 1e-10);
    let mut grad_x: Vec<[f64; 3]> = vec![[0.0; 3]; mesh.num_faces()];
    for (fi, f) in mesh.faces.iter().enumerate() {
        let p0 = mesh.vertices[f[0]];
        let p1 = mesh.vertices[f[1]];
        let p2 = mesh.vertices[f[2]];
        let fn_ = cross3(sub3(p1, p0), sub3(p2, p0));
        let area2 = len3(fn_);
        if area2 < 1e-30 {
            continue;
        }
        let n_hat = scale3(fn_, 1.0 / area2);
        let u0 = u[f[0]];
        let u1 = u[f[1]];
        let u2 = u[f[2]];
        let e0 = sub3(p2, p1);
        let e1 = sub3(p0, p2);
        let e2 = sub3(p1, p0);
        let g = add3(
            add3(scale3(cross3(n_hat, e0), u0), scale3(cross3(n_hat, e1), u1)),
            scale3(cross3(n_hat, e2), u2),
        );
        let gl = len3(g);
        grad_x[fi] = if gl > 1e-30 {
            scale3(g, -1.0 / gl)
        } else {
            [0.0; 3]
        };
    }
    let mut div = vec![0.0; n];
    for (fi, f) in mesh.faces.iter().enumerate() {
        let x = grad_x[fi];
        for k in 0..3 {
            let i = f[k];
            let j = f[(k + 1) % 3];
            let o = f[(k + 2) % 3];
            let ei = sub3(mesh.vertices[i], mesh.vertices[o]);
            let ej = sub3(mesh.vertices[j], mesh.vertices[o]);
            let cos_o = dot3(ei, ej) / (len3(ei) * len3(ej) + 1e-30);
            let sin_o = len3(cross3(ei, ej)) / (len3(ei) * len3(ej) + 1e-30);
            let cot_o = if sin_o.abs() > 1e-15 {
                cos_o / sin_o
            } else {
                0.0
            };
            let edge_ij = sub3(mesh.vertices[j], mesh.vertices[i]);
            div[i] += 0.5 * cot_o * dot3(x, edge_ij);
        }
    }
    let phi = jacobi_solve(&lap, &div, 500, 1e-10);
    let min_phi = phi.iter().copied().fold(f64::INFINITY, f64::min);
    phi.iter().map(|&p| (p - min_phi).max(0.0)).collect()
}
/// Compute a geodesic Voronoi diagram on a mesh from a set of site vertices.
///
/// Each mesh vertex is assigned to the nearest site (by geodesic distance).
/// Returns `(assignment, cells)` where `assignment[i]` is the index into `sites`
/// of the nearest site, and `cells` contains the Voronoi cell info.
pub fn geodesic_voronoi(
    mesh: &GeodesicMesh,
    sites: &[usize],
) -> (Vec<usize>, Vec<GeodesicVoronoiCell>) {
    let n = mesh.num_vertices();
    let adj = mesh.vertex_adjacency();
    let num_sites = sites.len();
    let mut dist = vec![f64::INFINITY; n];
    let mut label = vec![0usize; n];
    let mut heap = BinaryHeap::new();
    for (si, &s) in sites.iter().enumerate() {
        dist[s] = 0.0;
        label[s] = si;
        heap.push(DistNode {
            vertex: s,
            dist: 0.0,
        });
    }
    while let Some(DistNode { vertex: u, dist: d }) = heap.pop() {
        if d > dist[u] {
            continue;
        }
        for &v in &adj[u] {
            let w = dist3(mesh.vertices[u], mesh.vertices[v]);
            let new_dist = dist[u] + w;
            if new_dist < dist[v] {
                dist[v] = new_dist;
                label[v] = label[u];
                heap.push(DistNode {
                    vertex: v,
                    dist: new_dist,
                });
            }
        }
    }
    let mut cells: Vec<GeodesicVoronoiCell> = (0..num_sites)
        .map(|si| GeodesicVoronoiCell {
            site: sites[si],
            vertices: Vec::new(),
            area: 0.0,
        })
        .collect();
    for i in 0..n {
        let si = label[i];
        if si < num_sites {
            cells[si].vertices.push(i);
            cells[si].area += mesh.voronoi_area(i);
        }
    }
    (label, cells)
}
/// Compute the discrete exponential map from a source vertex.
///
/// Given a tangent vector `v` in the tangent plane at `source`, trace a
/// geodesic-like path on the surface and return where it lands.
///
/// This uses an iterative unfolding approach: the tangent vector is projected
/// onto successive triangles along the traced direction.
#[allow(clippy::too_many_arguments)]
pub fn discrete_exp_map(
    mesh: &GeodesicMesh,
    source: usize,
    tangent_vec: [f64; 2],
    max_steps: usize,
) -> ExpMapResult {
    let (t_axis, b_axis, _n) = mesh.tangent_frame(source);
    let v2f = mesh.vertex_to_faces();
    let dir3d = add3(
        scale3(t_axis, tangent_vec[0]),
        scale3(b_axis, tangent_vec[1]),
    );
    let total_len = len3(dir3d);
    if total_len < 1e-15 || v2f[source].is_empty() {
        return ExpMapResult {
            vertex: source,
            position: mesh.vertices[source],
            face: v2f[source].first().copied().unwrap_or(0),
            bary: [1.0, 0.0, 0.0],
        };
    }
    let dir = normalize3(dir3d);
    let mut remaining = total_len;
    let mut current_pos = mesh.vertices[source];
    let mut current_face = v2f[source][0];
    let mut current_vertex = source;
    let mut best_dot = f64::NEG_INFINITY;
    for &fi in &v2f[source] {
        let f = mesh.faces[fi];
        let centroid = scale3(
            add3(
                add3(mesh.vertices[f[0]], mesh.vertices[f[1]]),
                mesh.vertices[f[2]],
            ),
            1.0 / 3.0,
        );
        let to_centroid = normalize3(sub3(centroid, mesh.vertices[source]));
        let d = dot3(to_centroid, dir);
        if d > best_dot {
            best_dot = d;
            current_face = fi;
        }
    }
    for _step in 0..max_steps {
        if remaining < 1e-15 {
            break;
        }
        let f = mesh.faces[current_face];
        let p0 = mesh.vertices[f[0]];
        let p1 = mesh.vertices[f[1]];
        let p2 = mesh.vertices[f[2]];
        let fn_ = normalize3(cross3(sub3(p1, p0), sub3(p2, p0)));
        let proj_dir = normalize3(sub3(dir, scale3(fn_, dot3(dir, fn_))));
        let trial_pos = add3(current_pos, scale3(proj_dir, remaining));
        let bary = barycentric_coords(p0, p1, p2, trial_pos);
        if bary[0] >= -1e-10 && bary[1] >= -1e-10 && bary[2] >= -1e-10 {
            let clamped = clamp_bary(bary);
            let final_pos = add3(
                add3(scale3(p0, clamped[0]), scale3(p1, clamped[1])),
                scale3(p2, clamped[2]),
            );
            let mut best_v = f[0];
            let mut best_d2 = f64::INFINITY;
            for &vi in &f {
                let d2 = dot3(
                    sub3(mesh.vertices[vi], final_pos),
                    sub3(mesh.vertices[vi], final_pos),
                );
                if d2 < best_d2 {
                    best_d2 = d2;
                    best_v = vi;
                }
            }
            return ExpMapResult {
                vertex: best_v,
                position: final_pos,
                face: current_face,
                bary: clamped,
            };
        }
        let edges = [(1, 2, 0), (2, 0, 1), (0, 1, 2)];
        let mut crossed_edge = None;
        for &(ev0, ev1, _bi) in &edges {
            let edge_start = mesh.vertices[f[ev0]];
            let edge_end = mesh.vertices[f[ev1]];
            if let Some(t) = ray_edge_intersect(current_pos, proj_dir, edge_start, edge_end)
                && t > 1e-10
                && t <= remaining + 1e-10
            {
                crossed_edge = Some((ev0, ev1, t));
                break;
            }
        }
        if let Some((ev0, ev1, t)) = crossed_edge {
            let step_dist = t.min(remaining);
            current_pos = add3(current_pos, scale3(proj_dir, step_dist));
            remaining -= step_dist;
            let edge_a = f[ev0];
            let edge_b = f[ev1];
            let mut found_next = false;
            for fj in 0..mesh.num_faces() {
                if fj == current_face {
                    continue;
                }
                let fv = mesh.faces[fj];
                if fv.contains(&edge_a) && fv.contains(&edge_b) {
                    current_face = fj;
                    found_next = true;
                    break;
                }
            }
            if !found_next {
                break;
            }
        } else {
            break;
        }
        let mut best_d2 = f64::INFINITY;
        for &vi in &mesh.faces[current_face] {
            let d2 = dot3(
                sub3(mesh.vertices[vi], current_pos),
                sub3(mesh.vertices[vi], current_pos),
            );
            if d2 < best_d2 {
                best_d2 = d2;
                current_vertex = vi;
            }
        }
    }
    let f = mesh.faces[current_face];
    let bary = barycentric_coords(
        mesh.vertices[f[0]],
        mesh.vertices[f[1]],
        mesh.vertices[f[2]],
        current_pos,
    );
    ExpMapResult {
        vertex: current_vertex,
        position: current_pos,
        face: current_face,
        bary: clamp_bary(bary),
    }
}
/// Compute barycentric coordinates of point `p` with respect to triangle `(a, b, c)`.
pub(super) fn barycentric_coords(a: [f64; 3], b: [f64; 3], c: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let v0 = sub3(b, a);
    let v1 = sub3(c, a);
    let v2 = sub3(p, a);
    let d00 = dot3(v0, v0);
    let d01 = dot3(v0, v1);
    let d11 = dot3(v1, v1);
    let d20 = dot3(v2, v0);
    let d21 = dot3(v2, v1);
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-30 {
        return [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;
    [u, v, w]
}
/// Clamp barycentric coordinates to be non-negative and sum to 1.
pub(super) fn clamp_bary(b: [f64; 3]) -> [f64; 3] {
    let mut r = [b[0].max(0.0), b[1].max(0.0), b[2].max(0.0)];
    let s = r[0] + r[1] + r[2];
    if s > 1e-30 {
        r[0] /= s;
        r[1] /= s;
        r[2] /= s;
    } else {
        r = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
    }
    r
}
/// Ray-edge intersection in 3D (projected onto the common plane).
/// Returns `Some(t)` where `t` is the distance along the ray to the intersection.
pub(super) fn ray_edge_intersect(
    origin: [f64; 3],
    dir: [f64; 3],
    edge_a: [f64; 3],
    edge_b: [f64; 3],
) -> Option<f64> {
    let edge = sub3(edge_b, edge_a);
    let oa = sub3(origin, edge_a);
    let dd = dot3(dir, dir);
    let de = dot3(dir, edge);
    let ee = dot3(edge, edge);
    let do_ = dot3(dir, oa);
    let eo = dot3(edge, oa);
    let denom = dd * ee - de * de;
    if denom.abs() < 1e-30 {
        return None;
    }
    let t = (de * eo - ee * do_) / denom;
    let s = (dd * eo - de * do_) / denom;
    if (-1e-10..=1.0 + 1e-10).contains(&s) && t > 0.0 {
        Some(t)
    } else {
        None
    }
}
/// Compute the discrete logarithmic map: given source and target vertices,
/// find the 2D tangent-plane coordinates at `source` that correspond to `target`.
///
/// This is the inverse of the exponential map: it unfolds the geodesic path
/// from `source` to `target` into the tangent plane.
pub fn discrete_log_map(mesh: &GeodesicMesh, source: usize, target: usize) -> LogMapResult {
    if source == target {
        return LogMapResult {
            coords: [0.0, 0.0],
            distance: 0.0,
            angle: 0.0,
        };
    }
    let dist = dijkstra_geodesic(mesh, &[source]);
    let geo_dist = dist[target];
    if geo_dist.is_infinite() {
        return LogMapResult {
            coords: [0.0, 0.0],
            distance: f64::INFINITY,
            angle: 0.0,
        };
    }
    let path = trace_geodesic_path(mesh, &dist, source, target);
    if path.len() < 2 {
        return LogMapResult {
            coords: [0.0, 0.0],
            distance: geo_dist,
            angle: 0.0,
        };
    }
    let (t_axis, b_axis, _n) = mesh.tangent_frame(source);
    let next_vertex = path[1];
    let direction = normalize3(sub3(mesh.vertices[next_vertex], mesh.vertices[source]));
    let tx = dot3(direction, t_axis);
    let ty = dot3(direction, b_axis);
    let angle = ty.atan2(tx);
    LogMapResult {
        coords: [geo_dist * angle.cos(), geo_dist * angle.sin()],
        distance: geo_dist,
        angle,
    }
}
/// Trace the shortest geodesic path from `source` to `target` given
/// precomputed distances from `source`.
///
/// Returns a sequence of vertex indices from `source` to `target`.
pub fn trace_geodesic_path(
    mesh: &GeodesicMesh,
    dist_from_source: &[f64],
    source: usize,
    target: usize,
) -> Vec<usize> {
    if source == target {
        return vec![source];
    }
    let adj = mesh.vertex_adjacency();
    let mut path = vec![target];
    let mut current = target;
    let max_iter = mesh.num_vertices();
    for _i in 0..max_iter {
        if current == source {
            break;
        }
        let mut best_nb = current;
        let mut best_d = dist_from_source[current];
        for &nb in &adj[current] {
            if dist_from_source[nb] < best_d {
                best_d = dist_from_source[nb];
                best_nb = nb;
            }
        }
        if best_nb == current {
            break;
        }
        current = best_nb;
        path.push(current);
    }
    path.reverse();
    path
}
/// Compute the length of a geodesic path (sum of edge lengths along the path).
pub fn geodesic_path_length(mesh: &GeodesicMesh, path: &[usize]) -> f64 {
    if path.len() < 2 {
        return 0.0;
    }
    let mut length = 0.0;
    for i in 0..path.len() - 1 {
        length += dist3(mesh.vertices[path[i]], mesh.vertices[path[i + 1]]);
    }
    length
}
/// Smooth a geodesic path by iteratively projecting interior vertices towards
/// the geodesic between their neighbors.
///
/// Returns a new smoothed path (same start and end vertices).
pub fn smooth_geodesic_path(
    mesh: &GeodesicMesh,
    path: &[usize],
    iterations: usize,
) -> Vec<[f64; 3]> {
    if path.len() < 3 {
        return path.iter().map(|&v| mesh.vertices[v]).collect();
    }
    let mut positions: Vec<[f64; 3]> = path.iter().map(|&v| mesh.vertices[v]).collect();
    for _iter in 0..iterations {
        let old = positions.clone();
        for i in 1..positions.len() - 1 {
            let mid = lerp3(old[i - 1], old[i + 1], 0.5);
            positions[i] = lerp3(old[i], mid, 0.3);
        }
    }
    positions
}
/// Transport a tangent vector along a path of vertices on the mesh surface.
///
/// Uses discrete parallel transport via Schild's ladder approximation:
/// at each step, the vector is rotated to align with the change in normals
/// between consecutive vertices.
///
/// - `mesh`: the triangle mesh.
/// - `path`: sequence of vertex indices defining the transport path.
/// - `initial_vector`: the initial tangent vector at `path[0]`.
///
/// Returns a vector of tangent vectors, one at each vertex along the path.
pub fn parallel_transport(
    mesh: &GeodesicMesh,
    path: &[usize],
    initial_vector: [f64; 3],
) -> Vec<[f64; 3]> {
    if path.is_empty() {
        return Vec::new();
    }
    let mut transported = vec![[0.0; 3]; path.len()];
    let n0 = mesh.vertex_normal(path[0]);
    let v0 = sub3(initial_vector, scale3(n0, dot3(initial_vector, n0)));
    transported[0] = v0;
    for i in 1..path.len() {
        let n_prev = mesh.vertex_normal(path[i - 1]);
        let n_curr = mesh.vertex_normal(path[i]);
        let v_prev = transported[i - 1];
        let v_rotated = rotate_vector_between_normals(v_prev, n_prev, n_curr);
        let projected = sub3(v_rotated, scale3(n_curr, dot3(v_rotated, n_curr)));
        transported[i] = projected;
    }
    transported
}
/// Rotate a vector by the rotation that maps normal `n1` to normal `n2`.
///
/// Uses Rodrigues' rotation formula around the axis `n1 x n2`.
pub(super) fn rotate_vector_between_normals(v: [f64; 3], n1: [f64; 3], n2: [f64; 3]) -> [f64; 3] {
    let axis = cross3(n1, n2);
    let sin_a = len3(axis);
    let cos_a = dot3(n1, n2).clamp(-1.0, 1.0);
    if sin_a < 1e-12 {
        if cos_a > 0.0 {
            return v;
        } else {
            return scale3(v, -1.0);
        }
    }
    let k = scale3(axis, 1.0 / sin_a);
    let kxv = cross3(k, v);
    let kdv = dot3(k, v);
    add3(
        add3(scale3(v, cos_a), scale3(kxv, sin_a)),
        scale3(k, kdv * (1.0 - cos_a)),
    )
}
/// Compute the total rotation angle accumulated during parallel transport
/// around a closed loop.
///
/// This is the holonomy, which equals the integral of Gaussian curvature
/// enclosed by the loop.
pub fn parallel_transport_holonomy(mesh: &GeodesicMesh, loop_path: &[usize]) -> f64 {
    if loop_path.len() < 3 {
        return 0.0;
    }
    let mut path = loop_path.to_vec();
    if path.first() != path.last() {
        path.push(path[0]);
    }
    let (t_axis, _b_axis, _n) = mesh.tangent_frame(path[0]);
    let transported = parallel_transport(mesh, &path, t_axis);
    if transported.len() < 2 {
        return 0.0;
    }
    let v_initial = transported[0];
    let v_final = transported[transported.len() - 1];
    let cos_angle = dot3(normalize3(v_initial), normalize3(v_final)).clamp(-1.0, 1.0);
    cos_angle.acos()
}
/// Perform Farthest Point Sampling on a mesh.
///
/// Starting from `initial_vertex`, iteratively select the vertex that is
/// farthest from all previously selected vertices (using geodesic distances).
///
/// Returns a vector of `num_samples` vertex indices.
pub fn farthest_point_sampling(
    mesh: &GeodesicMesh,
    initial_vertex: usize,
    num_samples: usize,
) -> Vec<usize> {
    let n = mesh.num_vertices();
    if num_samples == 0 || n == 0 {
        return Vec::new();
    }
    let mut samples = Vec::with_capacity(num_samples);
    samples.push(initial_vertex);
    let mut min_dist = dijkstra_geodesic(mesh, &[initial_vertex]);
    for _ in 1..num_samples {
        let mut best_v = 0;
        let mut best_d = f64::NEG_INFINITY;
        for (v, &d) in min_dist.iter().enumerate() {
            if d > best_d && d.is_finite() {
                best_d = d;
                best_v = v;
            }
        }
        if best_d <= 0.0 || !best_d.is_finite() {
            break;
        }
        samples.push(best_v);
        let new_dist = dijkstra_geodesic(mesh, &[best_v]);
        for (v, d) in min_dist.iter_mut().enumerate() {
            *d = d.min(new_dist[v]);
        }
    }
    samples
}
/// Perform FPS using the fast marching method for better distance accuracy.
pub fn farthest_point_sampling_fmm(
    mesh: &GeodesicMesh,
    initial_vertex: usize,
    num_samples: usize,
) -> Vec<usize> {
    let n = mesh.num_vertices();
    if num_samples == 0 || n == 0 {
        return Vec::new();
    }
    let mut samples = Vec::with_capacity(num_samples);
    samples.push(initial_vertex);
    let mut min_dist = fast_marching_geodesic(mesh, &[initial_vertex]);
    for _ in 1..num_samples {
        let mut best_v = 0;
        let mut best_d = f64::NEG_INFINITY;
        for (v, &d) in min_dist.iter().enumerate() {
            if d > best_d && d.is_finite() {
                best_d = d;
                best_v = v;
            }
        }
        if best_d <= 0.0 || !best_d.is_finite() {
            break;
        }
        samples.push(best_v);
        let new_dist = fast_marching_geodesic(mesh, &[best_v]);
        for (v, d) in min_dist.iter_mut().enumerate() {
            *d = d.min(new_dist[v]);
        }
    }
    samples
}
/// Compute the geodesic centroid of a set of vertices on a mesh.
///
/// The geodesic centroid is the vertex that minimizes the sum of squared
/// geodesic distances to all points in the set. Found via iterative evaluation.
///
/// Returns the vertex index of the centroid.
pub fn geodesic_centroid(mesh: &GeodesicMesh, points: &[usize]) -> usize {
    if points.is_empty() {
        return 0;
    }
    let mut best_vertex = points[0];
    let mut best_cost = f64::INFINITY;
    for &candidate in points {
        let dist = dijkstra_geodesic(mesh, &[candidate]);
        let cost: f64 = points.iter().map(|&p| dist[p] * dist[p]).sum();
        if cost < best_cost {
            best_cost = cost;
            best_vertex = candidate;
        }
    }
    best_vertex
}
/// Compute the weighted geodesic centroid.
///
/// Each point has an associated weight; the centroid minimizes the weighted
/// sum of squared geodesic distances.
pub fn weighted_geodesic_centroid(mesh: &GeodesicMesh, points: &[usize], weights: &[f64]) -> usize {
    if points.is_empty() {
        return 0;
    }
    let mut best_vertex = points[0];
    let mut best_cost = f64::INFINITY;
    for &candidate in points {
        let dist = dijkstra_geodesic(mesh, &[candidate]);
        let cost: f64 = points
            .iter()
            .zip(weights.iter())
            .map(|(&p, &w)| w * dist[p] * dist[p])
            .sum();
        if cost < best_cost {
            best_cost = cost;
            best_vertex = candidate;
        }
    }
    best_vertex
}
/// Extract isoline vertices at a given geodesic distance from the source.
///
/// Returns a list of 3D positions on mesh edges where the distance field
/// crosses the specified `iso_distance`.
pub fn geodesic_isolines(
    mesh: &GeodesicMesh,
    dist_field: &[f64],
    iso_distance: f64,
) -> Vec<[f64; 3]> {
    let mut iso_points = Vec::new();
    for f in &mesh.faces {
        for k in 0..3 {
            let i = f[k];
            let j = f[(k + 1) % 3];
            let di = dist_field[i];
            let dj = dist_field[j];
            if (di - iso_distance) * (dj - iso_distance) < 0.0 {
                let t = (iso_distance - di) / (dj - di);
                let p = lerp3(mesh.vertices[i], mesh.vertices[j], t);
                iso_points.push(p);
            }
        }
    }
    iso_points
}
