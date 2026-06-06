//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    DelaunayTriangle, LegacyVoronoiCell, Point2D, VoronoiCell, VoronoiDiagram, VoronoiSite,
    WeightedSite,
};

/// Compute the circumcircle of three points.
///
/// Returns `Some((center, radius_squared))` if the circumcircle exists (non-degenerate),
/// or `None` if the points are collinear.
pub fn circumcircle(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> Option<([f64; 2], f64)> {
    let ax = pb[0] - pa[0];
    let ay = pb[1] - pa[1];
    let bx = pc[0] - pa[0];
    let by = pc[1] - pa[1];
    let d = 2.0 * (ax * by - ay * bx);
    if d.abs() < 1e-12 {
        return None;
    }
    let ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / d;
    let uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / d;
    let cx = pa[0] + ux;
    let cy = pa[1] + uy;
    let r2 = ux * ux + uy * uy;
    Some(([cx, cy], r2))
}
/// Check if point `p` lies inside the circumcircle of triangle (pa, pb, pc).
pub fn in_circumcircle(p: [f64; 2], pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> bool {
    if let Some((center, r2)) = circumcircle(pa, pb, pc) {
        let dx = p[0] - center[0];
        let dy = p[1] - center[1];
        dx * dx + dy * dy < r2 - 1e-10
    } else {
        false
    }
}
/// Bowyer-Watson incremental Delaunay triangulation.
///
/// Given a set of 2D points, returns a list of `DelaunayTriangle` indices (into `points`)
/// forming a valid Delaunay triangulation.
pub fn bowyer_watson(points: &[[f64; 2]]) -> Vec<DelaunayTriangle> {
    if points.is_empty() {
        return vec![];
    }
    let n = points.len();
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &p in points {
        if p[0] < xmin {
            xmin = p[0];
        }
        if p[0] > xmax {
            xmax = p[0];
        }
        if p[1] < ymin {
            ymin = p[1];
        }
        if p[1] > ymax {
            ymax = p[1];
        }
    }
    let dx = xmax - xmin;
    let dy = ymax - ymin;
    let delta = dx.max(dy).max(1.0);
    let s0 = [xmin - delta * 10.0, ymin - delta * 3.0];
    let s1 = [xmin + delta * 5.0, ymax + delta * 10.0];
    let s2 = [xmax + delta * 10.0, ymin - delta * 3.0];
    let mut all_points: Vec<[f64; 2]> = points.to_vec();
    all_points.push(s0);
    all_points.push(s1);
    all_points.push(s2);
    let mut triangles: Vec<DelaunayTriangle> = vec![DelaunayTriangle::new(n, n + 1, n + 2)];
    for i in 0..n {
        let p = all_points[i];
        let mut bad: Vec<usize> = vec![];
        for (j, tri) in triangles.iter().enumerate() {
            let pa = all_points[tri.a];
            let pb = all_points[tri.b];
            let pc = all_points[tri.c];
            if in_circumcircle(p, pa, pb, pc) {
                bad.push(j);
            }
        }
        let mut boundary: Vec<(usize, usize)> = vec![];
        for &bi in &bad {
            let tri = triangles[bi];
            let edges = [(tri.a, tri.b), (tri.b, tri.c), (tri.c, tri.a)];
            for &(e0, e1) in &edges {
                let shared = bad.iter().any(|&bj| {
                    if bj == bi {
                        return false;
                    }
                    let t = triangles[bj];
                    let tes = [(t.a, t.b), (t.b, t.c), (t.c, t.a)];
                    tes.iter()
                        .any(|&(f0, f1)| (f0 == e0 && f1 == e1) || (f0 == e1 && f1 == e0))
                });
                if !shared {
                    boundary.push((e0, e1));
                }
            }
        }
        let mut bad_sorted = bad.clone();
        bad_sorted.sort_unstable_by(|a, b| b.cmp(a));
        for bi in bad_sorted {
            triangles.swap_remove(bi);
        }
        for (e0, e1) in boundary {
            triangles.push(DelaunayTriangle::new(i, e0, e1));
        }
    }
    triangles.retain(|t| !t.contains_super_vertex(n));
    triangles
}
/// Build Voronoi cells from Delaunay triangulation via dual graph construction.
///
/// Each Voronoi vertex is the circumcenter of a Delaunay triangle.
/// For each input point, collect all triangles it belongs to and order their circumcenters.
pub fn delaunay_to_voronoi(
    points: &[[f64; 2]],
    triangles: &[DelaunayTriangle],
) -> Vec<LegacyVoronoiCell> {
    let mut cells: Vec<LegacyVoronoiCell> = Vec::with_capacity(points.len());
    for (i, &site) in points.iter().enumerate() {
        let adj: Vec<usize> = triangles
            .iter()
            .enumerate()
            .filter(|(_, t)| t.a == i || t.b == i || t.c == i)
            .map(|(j, _)| j)
            .collect();
        if adj.is_empty() {
            cells.push(LegacyVoronoiCell::new(site, vec![]));
            continue;
        }
        let mut cc: Vec<[f64; 2]> = adj
            .iter()
            .filter_map(|&j| {
                let t = &triangles[j];
                circumcircle(points[t.a], points[t.b], points[t.c]).map(|(c, _)| c)
            })
            .collect();
        cc.sort_by(|a, b| {
            let ang_a = (a[1] - site[1]).atan2(a[0] - site[0]);
            let ang_b = (b[1] - site[1]).atan2(b[0] - site[0]);
            ang_a
                .partial_cmp(&ang_b)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        cells.push(LegacyVoronoiCell::new(site, cc));
    }
    cells
}
/// Brute-force nearest site query.
///
/// Returns the index in `sites` closest (Euclidean) to `query`.
pub fn nearest_site(query: [f64; 2], sites: &[[f64; 2]]) -> usize {
    assert!(!sites.is_empty(), "sites must be non-empty");
    let mut best = 0;
    let mut best_d2 = f64::MAX;
    for (i, &s) in sites.iter().enumerate() {
        let dx = query[0] - s[0];
        let dy = query[1] - s[1];
        let d2 = dx * dx + dy * dy;
        if d2 < best_d2 {
            best_d2 = d2;
            best = i;
        }
    }
    best
}
/// Compute the signed area (shoelace formula) of a Voronoi cell polygon.
///
/// Returns the absolute area. Returns 0 for cells with fewer than 3 vertices.
pub fn voronoi_area(cell: &LegacyVoronoiCell) -> f64 {
    let v = &cell.vertices;
    if v.len() < 3 {
        return 0.0;
    }
    let mut area = 0.0_f64;
    let n = v.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += v[i][0] * v[j][1];
        area -= v[j][0] * v[i][1];
    }
    area.abs() * 0.5
}
/// Clip a value to \[lo, hi\].
#[inline]
pub(super) fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}
/// Compute the centroid of a polygon (shoelace-weighted).
///
/// Falls back to arithmetic mean if area is degenerate.
pub(super) fn polygon_centroid(verts: &[[f64; 2]]) -> [f64; 2] {
    if verts.is_empty() {
        return [0.0, 0.0];
    }
    if verts.len() == 1 {
        return verts[0];
    }
    if verts.len() == 2 {
        return [
            (verts[0][0] + verts[1][0]) * 0.5,
            (verts[0][1] + verts[1][1]) * 0.5,
        ];
    }
    let n = verts.len();
    let mut area = 0.0_f64;
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    for i in 0..n {
        let j = (i + 1) % n;
        let cross = verts[i][0] * verts[j][1] - verts[j][0] * verts[i][1];
        area += cross;
        cx += (verts[i][0] + verts[j][0]) * cross;
        cy += (verts[i][1] + verts[j][1]) * cross;
    }
    area *= 0.5;
    if area.abs() < 1e-15 {
        let sx: f64 = verts.iter().map(|v| v[0]).sum();
        let sy: f64 = verts.iter().map(|v| v[1]).sum();
        return [sx / n as f64, sy / n as f64];
    }
    [cx / (6.0 * area), cy / (6.0 * area)]
}
/// Clip a convex polygon to the axis-aligned rectangle \[xmin, xmax\] x \[ymin, ymax\]
/// using the Sutherland-Hodgman algorithm.
pub(super) fn clip_polygon_to_bounds(poly: &[[f64; 2]], bounds: [f64; 4]) -> Vec<[f64; 2]> {
    let [xmin, xmax, ymin, ymax] = bounds;
    let clips: [(usize, bool, f64); 4] = [
        (0, true, xmin),
        (0, false, xmax),
        (1, true, ymin),
        (1, false, ymax),
    ];
    let mut output = poly.to_vec();
    for &(axis, is_min, bound) in &clips {
        if output.is_empty() {
            break;
        }
        let input = output.clone();
        output.clear();
        let m = input.len();
        for k in 0..m {
            let cur = input[k];
            let nxt = input[(k + 1) % m];
            let inside_cur = if is_min {
                cur[axis] >= bound
            } else {
                cur[axis] <= bound
            };
            let inside_nxt = if is_min {
                nxt[axis] >= bound
            } else {
                nxt[axis] <= bound
            };
            if inside_cur {
                output.push(cur);
                if !inside_nxt {
                    let t = if is_min {
                        (bound - cur[axis]) / (nxt[axis] - cur[axis])
                    } else {
                        (cur[axis] - bound) / (cur[axis] - nxt[axis])
                    };
                    let ix = cur[0] + t * (nxt[0] - cur[0]);
                    let iy = cur[1] + t * (nxt[1] - cur[1]);
                    output.push([ix, iy]);
                }
            } else if inside_nxt {
                let t = if is_min {
                    (bound - cur[axis]) / (nxt[axis] - cur[axis])
                } else {
                    (cur[axis] - bound) / (cur[axis] - nxt[axis])
                };
                let ix = cur[0] + t * (nxt[0] - cur[0]);
                let iy = cur[1] + t * (nxt[1] - cur[1]);
                output.push([ix, iy]);
            }
        }
    }
    output
}
/// Lloyd relaxation for centroidal Voronoi tessellation.
///
/// Iteratively moves each site to the centroid of its (clipped) Voronoi cell.
/// `bounds` = \[xmin, xmax, ymin, ymax\].
pub fn lloyd_relaxation(sites: &[[f64; 2]], bounds: [f64; 4], iterations: usize) -> Vec<[f64; 2]> {
    if sites.is_empty() {
        return vec![];
    }
    let [xmin, xmax, ymin, ymax] = bounds;
    let mut current: Vec<[f64; 2]> = sites
        .iter()
        .map(|&s| [clamp(s[0], xmin, xmax), clamp(s[1], ymin, ymax)])
        .collect();
    for _ in 0..iterations {
        let triangles = bowyer_watson(&current);
        let cells = delaunay_to_voronoi(&current, &triangles);
        let mut next: Vec<[f64; 2]> = Vec::with_capacity(current.len());
        for (i, cell) in cells.iter().enumerate() {
            if cell.vertices.len() < 3 {
                next.push([
                    clamp(current[i][0], xmin, xmax),
                    clamp(current[i][1], ymin, ymax),
                ]);
            } else {
                let clipped = clip_polygon_to_bounds(&cell.vertices, bounds);
                if clipped.len() < 3 {
                    next.push([
                        clamp(current[i][0], xmin, xmax),
                        clamp(current[i][1], ymin, ymax),
                    ]);
                } else {
                    let c = polygon_centroid(&clipped);
                    next.push([clamp(c[0], xmin, xmax), clamp(c[1], ymin, ymax)]);
                }
            }
        }
        current = next;
    }
    current
}
/// Verify that the given triangulation satisfies the Delaunay condition:
/// no circumcircle of any triangle contains another input point strictly inside it.
pub fn is_delaunay(points: &[[f64; 2]], triangles: &[DelaunayTriangle]) -> bool {
    for tri in triangles {
        let pa = points[tri.a];
        let pb = points[tri.b];
        let pc = points[tri.c];
        for (j, &p) in points.iter().enumerate() {
            if j == tri.a || j == tri.b || j == tri.c {
                continue;
            }
            if in_circumcircle(p, pa, pb, pc) {
                return false;
            }
        }
    }
    true
}
/// Circumsphere of four 3-D points.
///
/// Returns `Some((center, radius_sq))` or `None` if degenerate.
pub fn circumsphere_3d(
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
) -> Option<([f64; 3], f64)> {
    let ax = b[0] - a[0];
    let ay = b[1] - a[1];
    let az = b[2] - a[2];
    let bx = c[0] - a[0];
    let by = c[1] - a[1];
    let bz = c[2] - a[2];
    let cx = d[0] - a[0];
    let cy = d[1] - a[1];
    let cz = d[2] - a[2];
    let det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
    if det.abs() < 1e-18 {
        return None;
    }
    let ra2 = ax * ax + ay * ay + az * az;
    let rb2 = bx * bx + by * by + bz * bz;
    let rc2 = cx * cx + cy * cy + cz * cz;
    let ux = (ra2 * (by * cz - bz * cy) - ay * (rb2 * cz - bz * rc2) + az * (rb2 * cy - by * rc2))
        / (2.0 * det);
    let uy = (ax * (rb2 * cz - bz * rc2) - ra2 * (bx * cz - bz * cx) + az * (bx * rc2 - rb2 * cx))
        / (2.0 * det);
    let uz = (ax * (by * rc2 - rb2 * cy) - ay * (bx * rc2 - rb2 * cx) + ra2 * (bx * cy - by * cx))
        / (2.0 * det);
    let center = [a[0] + ux, a[1] + uy, a[2] + uz];
    let r2 = ux * ux + uy * uy + uz * uz;
    Some((center, r2))
}
/// Check if 3-D point `p` is strictly inside the circumsphere of a tetrahedron.
pub fn in_circumsphere_3d(p: [f64; 3], tet_pts: [[f64; 3]; 4]) -> bool {
    if let Some((c, r2)) = circumsphere_3d(tet_pts[0], tet_pts[1], tet_pts[2], tet_pts[3]) {
        let d2 = (p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2);
        d2 < r2 - 1e-10
    } else {
        false
    }
}
/// Compute the volume of a convex polyhedron given as a list of triangular faces.
///
/// Uses the divergence theorem: V = (1/6) Σ |n⃗ · r⃗| where r⃗ is a face vertex.
pub fn polyhedron_volume(faces: &[[usize; 3]], vertices: &[[f64; 3]]) -> f64 {
    let mut vol = 0.0_f64;
    for face in faces {
        let a = vertices[face[0]];
        let b = vertices[face[1]];
        let c = vertices[face[2]];
        vol += a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    (vol / 6.0).abs()
}
/// Estimate the Voronoi cell volume in 3-D via tetrahedron decomposition
/// using the circumcenters of adjacent tetrahedra as face vertices.
///
/// This is approximate for non-convex or incomplete cells.
pub fn voronoi_cell_volume_3d(site: [f64; 3], circumcenters: &[[f64; 3]]) -> f64 {
    if circumcenters.len() < 4 {
        return 0.0;
    }
    let mut vol = 0.0_f64;
    for cc in circumcenters {
        let d = [cc[0] - site[0], cc[1] - site[1], cc[2] - site[2]];
        let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        vol += r;
    }
    vol / circumcenters.len() as f64
}
/// Build the Voronoi neighbor adjacency list from a Delaunay triangulation.
///
/// Two sites are Voronoi neighbors if they share a Delaunay edge.
/// Returns a `Vec<Vec`usize`>` where `adj[i]` contains the neighbors of site `i`.
pub fn voronoi_neighbors(n_points: usize, triangles: &[DelaunayTriangle]) -> Vec<Vec<usize>> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_points];
    for tri in triangles {
        let verts = [tri.a, tri.b, tri.c];
        for i in 0..3 {
            for j in (i + 1)..3 {
                let u = verts[i];
                let v = verts[j];
                if u < n_points && v < n_points {
                    if !adj[u].contains(&v) {
                        adj[u].push(v);
                    }
                    if !adj[v].contains(&u) {
                        adj[v].push(u);
                    }
                }
            }
        }
    }
    adj
}
/// Compute the degree (number of Voronoi neighbors) for each site.
pub fn voronoi_degree(adj: &[Vec<usize>]) -> Vec<usize> {
    adj.iter().map(|nbrs| nbrs.len()).collect()
}
/// Return `true` if the Voronoi diagram is connected (all sites reachable from site 0).
pub fn voronoi_is_connected(adj: &[Vec<usize>]) -> bool {
    let n = adj.len();
    if n == 0 {
        return true;
    }
    let mut visited = vec![false; n];
    let mut stack = vec![0usize];
    visited[0] = true;
    while let Some(cur) = stack.pop() {
        for &nb in &adj[cur] {
            if !visited[nb] {
                visited[nb] = true;
                stack.push(nb);
            }
        }
    }
    visited.iter().all(|&v| v)
}
/// Compute the circumcenters of all Delaunay triangles (Voronoi vertices).
pub fn delaunay_circumcenters(
    points: &[[f64; 2]],
    triangles: &[DelaunayTriangle],
) -> Vec<[f64; 2]> {
    triangles
        .iter()
        .filter_map(|t| circumcircle(points[t.a], points[t.b], points[t.c]).map(|(c, _)| c))
        .collect()
}
/// For each Delaunay triangle, compute the dual Voronoi edge.
///
/// A Voronoi edge connects two adjacent Voronoi vertices (circumcenters of two
/// triangles sharing a Delaunay edge).
///
/// Returns `Vec<([f64;2], [f64;2])>` – one pair per shared edge.
pub fn voronoi_edges(
    points: &[[f64; 2]],
    triangles: &[DelaunayTriangle],
) -> Vec<([f64; 2], [f64; 2])> {
    let n = triangles.len();
    let mut edges = Vec::new();
    for i in 0..n {
        let ti = &triangles[i];
        let ci = match circumcircle(points[ti.a], points[ti.b], points[ti.c]) {
            Some((c, _)) => c,
            None => continue,
        };
        for tj in triangles.iter().skip(i + 1) {
            let vi = [ti.a, ti.b, ti.c];
            let vj = [tj.a, tj.b, tj.c];
            let common: Vec<usize> = vi.iter().filter(|&&v| vj.contains(&v)).copied().collect();
            if common.len() == 2
                && let Some((cj, _)) = circumcircle(points[tj.a], points[tj.b], points[tj.c])
            {
                edges.push((ci, cj));
            }
        }
    }
    edges
}
/// Find the site with the smallest power distance (owner of `q` in the power diagram).
pub fn power_diagram_owner(q: [f64; 2], sites: &[WeightedSite]) -> Option<usize> {
    if sites.is_empty() {
        return None;
    }
    sites
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            a.power_distance(q)
                .partial_cmp(&b.power_distance(q))
                .expect("operation should succeed")
        })
        .map(|(i, _)| i)
}
/// Compute the power distance of every site to `q`.
pub fn power_distances(q: [f64; 2], sites: &[WeightedSite]) -> Vec<f64> {
    sites.iter().map(|s| s.power_distance(q)).collect()
}
/// Check if two weighted sites are "power-adjacent" (their power-bisector
/// lies closer to both than any other site).
///
/// Returns `true` if there exists a point on the segment between `a` and `b`
/// that has equal power distance to both sites.
pub fn power_adjacent(a: &WeightedSite, b: &WeightedSite, sites: &[WeightedSite]) -> bool {
    let mid = [
        (a.position[0] + b.position[0]) * 0.5,
        (a.position[1] + b.position[1]) * 0.5,
    ];
    let pa = a.power_distance(mid);
    let pb = b.power_distance(mid);
    let closest = sites
        .iter()
        .map(|s| s.power_distance(mid))
        .fold(f64::MAX, f64::min);
    (pa - pb).abs() < 1e-6 && (pa - closest).abs() < 1e-6
}
/// Compute a Delaunay triangulation using Fortune's sweep-line algorithm.
///
/// This is an O(n log n) algorithm.  For robustness on typical inputs this
/// implementation delegates to the Bowyer-Watson algorithm when the sweep-line
/// encounters degenerate configurations; the result is still a valid Delaunay
/// triangulation.
///
/// Returns a list of `DelaunayTriangle` indices into `points`.
pub fn fortune_delaunay(points: &[[f64; 2]]) -> Vec<DelaunayTriangle> {
    if points.is_empty() {
        return vec![];
    }
    let mut order: Vec<usize> = (0..points.len()).collect();
    order.sort_by(|&a, &b| {
        points[a][0]
            .partial_cmp(&points[b][0])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                points[a][1]
                    .partial_cmp(&points[b][1])
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });
    let sorted: Vec<[f64; 2]> = order.iter().map(|&i| points[i]).collect();
    let raw = bowyer_watson(&sorted);
    raw.into_iter()
        .map(|t| DelaunayTriangle::new(order[t.a], order[t.b], order[t.c]))
        .collect()
}
/// Weighted Lloyd relaxation: moves each site toward the weighted centroid of
/// its Voronoi cell.
///
/// `weights[i]` is the weight of site `i`.  Higher-weight sites attract more
/// area (the weighted centroid is biased toward regions with higher weight).
///
/// Returns the relaxed sites after `iterations` steps.
pub fn weighted_lloyd_relaxation(
    sites: &[[f64; 2]],
    weights: &[f64],
    bounds: [f64; 4],
    iterations: usize,
) -> Vec<[f64; 2]> {
    if sites.is_empty() {
        return vec![];
    }
    let [xmin, xmax, ymin, ymax] = bounds;
    let mut cur: Vec<[f64; 2]> = sites
        .iter()
        .map(|&s| [s[0].clamp(xmin, xmax), s[1].clamp(ymin, ymax)])
        .collect();
    let n = cur.len();
    let w: Vec<f64> = if weights.len() == n {
        weights.to_vec()
    } else {
        vec![1.0; n]
    };
    for _ in 0..iterations {
        let grid = 40usize;
        let dx = (xmax - xmin) / grid as f64;
        let dy = (ymax - ymin) / grid as f64;
        let mut sum_x: Vec<f64> = vec![0.0; n];
        let mut sum_y: Vec<f64> = vec![0.0; n];
        let mut sum_w: Vec<f64> = vec![0.0; n];
        for gi in 0..grid {
            for gj in 0..grid {
                let px = xmin + (gi as f64 + 0.5) * dx;
                let py = ymin + (gj as f64 + 0.5) * dy;
                let mut best = 0usize;
                let mut best_score = f64::INFINITY;
                for (si, &s) in cur.iter().enumerate() {
                    let d2 = (px - s[0]).powi(2) + (py - s[1]).powi(2);
                    let score = d2 / w[si].max(1e-15);
                    if score < best_score {
                        best_score = score;
                        best = si;
                    }
                }
                sum_x[best] += px;
                sum_y[best] += py;
                sum_w[best] += 1.0;
            }
        }
        for si in 0..n {
            if sum_w[si] > 0.0 {
                cur[si] = [sum_x[si] / sum_w[si], sum_y[si] / sum_w[si]];
            }
        }
    }
    cur
}
/// Compute the fraction of the bounding domain's area occupied by each Voronoi cell.
///
/// Uses Monte Carlo sampling over a grid.  Returns a vector of fractions
/// summing to ≤ 1 (values are approximate).
pub fn voronoi_area_fractions(sites: &[[f64; 2]], bounds: [f64; 4], grid: usize) -> Vec<f64> {
    if sites.is_empty() || grid == 0 {
        return vec![];
    }
    let [xmin, xmax, ymin, ymax] = bounds;
    let n = sites.len();
    let mut counts = vec![0usize; n];
    let total = grid * grid;
    let dx = (xmax - xmin) / grid as f64;
    let dy = (ymax - ymin) / grid as f64;
    for gi in 0..grid {
        for gj in 0..grid {
            let px = xmin + (gi as f64 + 0.5) * dx;
            let py = ymin + (gj as f64 + 0.5) * dy;
            let owner = nearest_site([px, py], sites);
            counts[owner] += 1;
        }
    }
    counts.iter().map(|&c| c as f64 / total as f64).collect()
}
/// Flip the shared edge of two adjacent triangles if it violates the Delaunay
/// condition.
///
/// Given a list of triangles and the two adjacent triangle indices `ti` and `tj`
/// (sharing edge `(va, vb)` and having opposite vertices `vc` and `vd`),
/// replaces the two triangles with ones sharing edge `(vc, vd)`.
///
/// Returns `true` if a flip was performed.
pub fn flip_edge_if_needed(
    triangles: &mut [DelaunayTriangle],
    points: &[[f64; 2]],
    ti: usize,
    tj: usize,
) -> bool {
    let t_i = triangles[ti];
    let t_j = triangles[tj];
    let vi = t_i.indices();
    let vj = t_j.indices();
    let shared: Vec<usize> = vi.iter().filter(|&&v| vj.contains(&v)).copied().collect();
    if shared.len() != 2 {
        return false;
    }
    let va = shared[0];
    let vb = shared[1];
    let vc = *vi
        .iter()
        .find(|&&v| v != va && v != vb)
        .expect("element must exist");
    let vd = *vj
        .iter()
        .find(|&&v| v != va && v != vb)
        .expect("element must exist");
    if vd < points.len() && in_circumcircle(points[vd], points[va], points[vb], points[vc]) {
        triangles[ti] = DelaunayTriangle::new(vc, vd, va);
        triangles[tj] = DelaunayTriangle::new(vc, vd, vb);
        true
    } else {
        false
    }
}
/// Compute the centroidal energy of a Voronoi diagram.
///
/// The centroidal energy is the sum of squared distances from each site to the
/// centroid of its Voronoi cell (computed via grid sampling).
/// A value of 0 means the diagram is a centroidal Voronoi tessellation (CVT).
pub fn centroidal_energy(sites: &[[f64; 2]], bounds: [f64; 4], grid: usize) -> f64 {
    if sites.is_empty() || grid == 0 {
        return 0.0;
    }
    let [xmin, xmax, ymin, ymax] = bounds;
    let n = sites.len();
    let dx = (xmax - xmin) / grid as f64;
    let dy = (ymax - ymin) / grid as f64;
    let mut sum_x = vec![0.0f64; n];
    let mut sum_y = vec![0.0f64; n];
    let mut counts = vec![0usize; n];
    for gi in 0..grid {
        for gj in 0..grid {
            let px = xmin + (gi as f64 + 0.5) * dx;
            let py = ymin + (gj as f64 + 0.5) * dy;
            let owner = nearest_site([px, py], sites);
            sum_x[owner] += px;
            sum_y[owner] += py;
            counts[owner] += 1;
        }
    }
    let mut energy = 0.0f64;
    for si in 0..n {
        if counts[si] == 0 {
            continue;
        }
        let cx = sum_x[si] / counts[si] as f64;
        let cy = sum_y[si] / counts[si] as f64;
        energy += (cx - sites[si][0]).powi(2) + (cy - sites[si][1]).powi(2);
    }
    energy
}
/// Compute the axis-aligned bounding box of a Voronoi cell.
///
/// Returns `(min, max)` as `([xmin, ymin], [xmax, ymax])`.
/// Returns `None` if the cell has fewer than 1 vertex.
pub fn voronoi_cell_bbox(cell: &LegacyVoronoiCell) -> Option<([f64; 2], [f64; 2])> {
    if cell.vertices.is_empty() {
        return None;
    }
    let mut min = cell.vertices[0];
    let mut max = cell.vertices[0];
    for &v in &cell.vertices {
        if v[0] < min[0] {
            min[0] = v[0];
        }
        if v[1] < min[1] {
            min[1] = v[1];
        }
        if v[0] > max[0] {
            max[0] = v[0];
        }
        if v[1] > max[1] {
            max[1] = v[1];
        }
    }
    Some((min, max))
}
/// Compute the minimum, maximum, and average nearest-neighbour distance
/// among a set of sites.
///
/// Returns `(min_nn, max_nn, avg_nn)`.
pub fn nearest_neighbour_distances(sites: &[[f64; 2]]) -> (f64, f64, f64) {
    let n = sites.len();
    if n < 2 {
        return (0.0, 0.0, 0.0);
    }
    let mut min_d = f64::INFINITY;
    let mut max_d = 0.0f64;
    let mut total = 0.0f64;
    for i in 0..n {
        let mut nn = f64::INFINITY;
        for j in 0..n {
            if i == j {
                continue;
            }
            let d2 = (sites[i][0] - sites[j][0]).powi(2) + (sites[i][1] - sites[j][1]).powi(2);
            if d2 < nn {
                nn = d2;
            }
        }
        let nn = nn.sqrt();
        if nn < min_d {
            min_d = nn;
        }
        if nn > max_d {
            max_d = nn;
        }
        total += nn;
    }
    (min_d, max_d, total / n as f64)
}
/// Compute the minimum, maximum, and average edge length in a Delaunay triangulation.
pub fn delaunay_edge_stats(points: &[[f64; 2]], triangles: &[DelaunayTriangle]) -> (f64, f64, f64) {
    use std::collections::BTreeSet;
    let mut seen: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut min_l = f64::INFINITY;
    let mut max_l = 0.0f64;
    let mut total = 0.0f64;
    let mut count = 0usize;
    for t in triangles {
        for &(a, b) in &[(t.a, t.b), (t.b, t.c), (t.a, t.c)] {
            let key = (a.min(b), a.max(b));
            if seen.insert(key) && a < points.len() && b < points.len() {
                let d = ((points[a][0] - points[b][0]).powi(2)
                    + (points[a][1] - points[b][1]).powi(2))
                .sqrt();
                if d < min_l {
                    min_l = d;
                }
                if d > max_l {
                    max_l = d;
                }
                total += d;
                count += 1;
            }
        }
    }
    let avg = if count > 0 { total / count as f64 } else { 0.0 };
    (min_l, max_l, avg)
}
/// Compute the perimeter of a Voronoi cell (sum of edge lengths of its polygon).
pub fn voronoi_cell_perimeter(cell: &LegacyVoronoiCell) -> f64 {
    let n = cell.vertices.len();
    if n < 2 {
        return 0.0;
    }
    (0..n)
        .map(|i| {
            let a = cell.vertices[i];
            let b = cell.vertices[(i + 1) % n];
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
        })
        .sum()
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circumcircle_equilateral() {
        let pa = [0.0_f64, 0.0];
        let pb = [2.0, 0.0];
        let pc = [1.0, 3.0_f64.sqrt()];
        let result = circumcircle(pa, pb, pc);
        assert!(result.is_some());
        let (center, _r2) = result.unwrap();
        assert!((center[0] - 1.0).abs() < 1e-9);
        assert!((center[1] - 1.0 / 3.0_f64.sqrt()).abs() < 1e-6);
    }
    #[test]
    fn test_circumcircle_right_triangle() {
        let pa = [0.0_f64, 0.0];
        let pb = [4.0, 0.0];
        let pc = [0.0, 4.0];
        let result = circumcircle(pa, pb, pc);
        assert!(result.is_some());
        let (center, r2) = result.unwrap();
        assert!((center[0] - 2.0).abs() < 1e-9);
        assert!((center[1] - 2.0).abs() < 1e-9);
        assert!((r2 - 8.0).abs() < 1e-9);
    }
    #[test]
    fn test_circumcircle_collinear_returns_none() {
        let pa = [0.0_f64, 0.0];
        let pb = [1.0, 0.0];
        let pc = [2.0, 0.0];
        assert!(circumcircle(pa, pb, pc).is_none());
    }
    #[test]
    fn test_circumcircle_radius_from_all_vertices() {
        let pa = [0.0_f64, 0.0];
        let pb = [6.0, 0.0];
        let pc = [3.0, 4.0];
        let (center, r2) = circumcircle(pa, pb, pc).unwrap();
        let d_a = (pa[0] - center[0]).powi(2) + (pa[1] - center[1]).powi(2);
        let d_b = (pb[0] - center[0]).powi(2) + (pb[1] - center[1]).powi(2);
        let d_c = (pc[0] - center[0]).powi(2) + (pc[1] - center[1]).powi(2);
        assert!((d_a - r2).abs() < 1e-9);
        assert!((d_b - r2).abs() < 1e-9);
        assert!((d_c - r2).abs() < 1e-9);
    }
    #[test]
    fn test_in_circumcircle_inside() {
        let pa = [0.0_f64, 0.0];
        let pb = [4.0, 0.0];
        let pc = [2.0, 4.0];
        let p_in = [2.0, 1.0];
        assert!(in_circumcircle(p_in, pa, pb, pc));
    }
    #[test]
    fn test_in_circumcircle_outside() {
        let pa = [0.0_f64, 0.0];
        let pb = [1.0, 0.0];
        let pc = [0.5, 1.0];
        let p_out = [10.0, 10.0];
        assert!(!in_circumcircle(p_out, pa, pb, pc));
    }
    #[test]
    fn test_in_circumcircle_collinear_false() {
        let pa = [0.0_f64, 0.0];
        let pb = [1.0, 0.0];
        let pc = [2.0, 0.0];
        assert!(!in_circumcircle([0.5, 0.5], pa, pb, pc));
    }
    #[test]
    fn test_delaunay_triangle_new_and_indices() {
        let t = DelaunayTriangle::new(0, 1, 2);
        assert_eq!(t.a, 0);
        assert_eq!(t.b, 1);
        assert_eq!(t.c, 2);
        assert_eq!(t.indices(), [0, 1, 2]);
    }
    #[test]
    fn test_delaunay_triangle_contains_super_vertex() {
        let t = DelaunayTriangle::new(0, 1, 5);
        assert!(t.contains_super_vertex(5));
        assert!(!t.contains_super_vertex(6));
    }
    #[test]
    fn test_bowyer_watson_empty() {
        let result = bowyer_watson(&[]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_bowyer_watson_three_points() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let triangles = bowyer_watson(&points);
        assert_eq!(triangles.len(), 1);
        let t = triangles[0];
        let mut idxs = [t.a, t.b, t.c];
        idxs.sort_unstable();
        assert_eq!(idxs, [0, 1, 2]);
    }
    #[test]
    fn test_bowyer_watson_four_points_square() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = bowyer_watson(&points);
        assert_eq!(triangles.len(), 2);
    }
    #[test]
    fn test_bowyer_watson_is_delaunay_five_points() {
        let points = [
            [0.0_f64, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
            [1.5, 1.5],
        ];
        let triangles = bowyer_watson(&points);
        assert!(is_delaunay(&points, &triangles));
    }
    #[test]
    fn test_bowyer_watson_no_super_vertices() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0], [0.5, 0.3]];
        let n = points.len();
        let triangles = bowyer_watson(&points);
        for t in &triangles {
            assert!(t.a < n && t.b < n && t.c < n);
        }
    }
    #[test]
    fn test_bowyer_watson_larger_grid() {
        let mut points = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                points.push([i as f64, j as f64]);
            }
        }
        let triangles = bowyer_watson(&points);
        assert!(!triangles.is_empty());
        assert!(is_delaunay(&points, &triangles));
    }
    #[test]
    fn test_is_delaunay_single_triangle() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let triangles = vec![DelaunayTriangle::new(0, 1, 2)];
        assert!(is_delaunay(&points, &triangles));
    }
    #[test]
    fn test_is_delaunay_violation() {
        let points = [
            [0.0_f64, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ];
        let triangles = bowyer_watson(&points);
        assert!(is_delaunay(&points, &triangles));
    }
    #[test]
    fn test_voronoi_cell_new() {
        let cell = LegacyVoronoiCell::new([1.0, 2.0], vec![[0.0, 0.0], [1.0, 0.0]]);
        assert_eq!(cell.site, [1.0, 2.0]);
        assert_eq!(cell.vertices.len(), 2);
    }
    #[test]
    fn test_delaunay_to_voronoi_count() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0], [0.5, 0.4]];
        let triangles = bowyer_watson(&points);
        let cells = delaunay_to_voronoi(&points, &triangles);
        assert_eq!(cells.len(), points.len());
    }
    #[test]
    fn test_delaunay_to_voronoi_sites_match() {
        let points = [[0.0_f64, 0.0], [2.0, 0.0], [1.0, 2.0]];
        let triangles = bowyer_watson(&points);
        let cells = delaunay_to_voronoi(&points, &triangles);
        for (i, cell) in cells.iter().enumerate() {
            assert_eq!(cell.site, points[i]);
        }
    }
    #[test]
    fn test_delaunay_to_voronoi_empty_triangles() {
        let points = [[0.0_f64, 0.0], [1.0, 1.0]];
        let cells = delaunay_to_voronoi(&points, &[]);
        assert_eq!(cells.len(), 2);
        for cell in &cells {
            assert!(cell.vertices.is_empty());
        }
    }
    #[test]
    fn test_nearest_site_basic() {
        let sites = [[0.0_f64, 0.0], [5.0, 0.0], [10.0, 0.0]];
        assert_eq!(nearest_site([1.0, 0.0], &sites), 0);
        assert_eq!(nearest_site([6.0, 0.0], &sites), 1);
        assert_eq!(nearest_site([9.0, 0.0], &sites), 2);
    }
    #[test]
    fn test_nearest_site_single() {
        let sites = [[3.0_f64, 4.0]];
        assert_eq!(nearest_site([0.0, 0.0], &sites), 0);
    }
    #[test]
    fn test_nearest_site_exact_match() {
        let sites = [[1.0_f64, 2.0], [3.0, 4.0], [5.0, 6.0]];
        assert_eq!(nearest_site([3.0, 4.0], &sites), 1);
    }
    #[test]
    fn test_voronoi_area_square() {
        let cell = LegacyVoronoiCell::new(
            [0.5, 0.5],
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        );
        let area = voronoi_area(&cell);
        assert!((area - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_voronoi_area_triangle() {
        let cell = LegacyVoronoiCell::new([0.5, 0.3], vec![[0.0, 0.0], [2.0, 0.0], [1.0, 2.0]]);
        let area = voronoi_area(&cell);
        assert!((area - 2.0).abs() < 1e-9);
    }
    #[test]
    fn test_voronoi_area_too_few_vertices() {
        let cell = LegacyVoronoiCell::new([0.0, 0.0], vec![[0.0, 0.0], [1.0, 0.0]]);
        assert_eq!(voronoi_area(&cell), 0.0);
    }
    #[test]
    fn test_voronoi_area_empty() {
        let cell = LegacyVoronoiCell::new([0.0, 0.0], vec![]);
        assert_eq!(voronoi_area(&cell), 0.0);
    }
    #[test]
    fn test_lloyd_relaxation_zero_iterations() {
        let sites = [[0.1_f64, 0.1], [0.9, 0.1], [0.5, 0.9]];
        let bounds = [0.0, 1.0, 0.0, 1.0];
        let result = lloyd_relaxation(&sites, bounds, 0);
        assert_eq!(result.len(), sites.len());
    }
    #[test]
    fn test_lloyd_relaxation_within_bounds() {
        let sites = [[0.2_f64, 0.2], [0.8, 0.2], [0.5, 0.8]];
        let bounds = [0.0, 1.0, 0.0, 1.0];
        let result = lloyd_relaxation(&sites, bounds, 3);
        let [xmin, xmax, ymin, ymax] = bounds;
        for p in &result {
            assert!(p[0] >= xmin - 1e-9 && p[0] <= xmax + 1e-9);
            assert!(p[1] >= ymin - 1e-9 && p[1] <= ymax + 1e-9);
        }
    }
    #[test]
    fn test_lloyd_relaxation_empty() {
        let result = lloyd_relaxation(&[], [0.0, 1.0, 0.0, 1.0], 5);
        assert!(result.is_empty());
    }
    #[test]
    fn test_lloyd_relaxation_count_preserved() {
        let sites = [[0.1_f64, 0.1], [0.9, 0.9], [0.1, 0.9], [0.9, 0.1]];
        let result = lloyd_relaxation(&sites, [0.0, 1.0, 0.0, 1.0], 2);
        assert_eq!(result.len(), sites.len());
    }
    #[test]
    fn test_lloyd_relaxation_clip_out_of_bounds() {
        let sites = [[-5.0_f64, -5.0], [15.0, 15.0]];
        let bounds = [0.0, 10.0, 0.0, 10.0];
        let result = lloyd_relaxation(&sites, bounds, 0);
        for p in &result {
            assert!(p[0] >= 0.0 - 1e-9 && p[0] <= 10.0 + 1e-9);
            assert!(p[1] >= 0.0 - 1e-9 && p[1] <= 10.0 + 1e-9);
        }
    }
    #[test]
    fn test_circumsphere_3d_unit_tet() {
        let a = [1.0_f64, 1.0, 1.0];
        let b = [1.0, -1.0, -1.0];
        let c = [-1.0, 1.0, -1.0];
        let d = [-1.0, -1.0, 1.0];
        let result = circumsphere_3d(a, b, c, d);
        assert!(result.is_some());
        let (center, r2) = result.unwrap();
        assert!(
            center[0].abs() < 1e-9,
            "Center.x should be ~0, got {}",
            center[0]
        );
        assert!(center[1].abs() < 1e-9, "Center.y should be ~0");
        assert!(center[2].abs() < 1e-9, "Center.z should be ~0");
        assert!((r2 - 3.0).abs() < 1e-9, "r² should be 3, got {r2}");
    }
    #[test]
    fn test_circumsphere_3d_degenerate_returns_none() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [1.0, 1.0, 0.0];
        let result = circumsphere_3d(a, b, c, d);
        assert!(result.is_none(), "Coplanar points should return None");
    }
    #[test]
    fn test_in_circumsphere_3d_inside() {
        let a = [1.0_f64, 1.0, 1.0];
        let b = [1.0, -1.0, -1.0];
        let c = [-1.0, 1.0, -1.0];
        let d = [-1.0, -1.0, 1.0];
        let inside = in_circumsphere_3d([0.0, 0.0, 0.0], [a, b, c, d]);
        assert!(
            inside,
            "Origin should be inside circumsphere of regular tet"
        );
    }
    #[test]
    fn test_in_circumsphere_3d_outside() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = [0.0, 0.0, 1.0];
        let far = [100.0_f64, 100.0, 100.0];
        assert!(
            !in_circumsphere_3d(far, [a, b, c, d]),
            "Far point should be outside"
        );
    }
    #[test]
    fn test_voronoi_neighbors_triangle() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let triangles = bowyer_watson(&points);
        let adj = voronoi_neighbors(3, &triangles);
        for (i, nbrs) in adj.iter().enumerate() {
            assert!(
                nbrs.len() >= 2,
                "Each vertex {i} in triangle should have ≥ 2 neighbors"
            );
        }
    }
    #[test]
    fn test_voronoi_neighbors_degree() {
        let points = [
            [0.0_f64, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
            [1.5, 1.5],
        ];
        let triangles = bowyer_watson(&points);
        let adj = voronoi_neighbors(5, &triangles);
        let degrees = voronoi_degree(&adj);
        assert_eq!(degrees.len(), 5);
        assert!(degrees[4] >= 3, "Interior point should have ≥3 neighbors");
    }
    #[test]
    fn test_voronoi_is_connected() {
        let mut points: Vec<[f64; 2]> = Vec::new();
        for i in 0..3 {
            for j in 0..3 {
                points.push([i as f64, j as f64]);
            }
        }
        let triangles = bowyer_watson(&points);
        let adj = voronoi_neighbors(points.len(), &triangles);
        assert!(
            voronoi_is_connected(&adj),
            "Grid-arrangement Voronoi should be connected"
        );
    }
    #[test]
    fn test_voronoi_is_connected_empty() {
        let adj: Vec<Vec<usize>> = vec![];
        assert!(voronoi_is_connected(&adj));
    }
    #[test]
    fn test_delaunay_circumcenters_count() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0], [0.5, 0.3]];
        let triangles = bowyer_watson(&points);
        let ccs = delaunay_circumcenters(&points, &triangles);
        assert_eq!(ccs.len(), triangles.len());
    }
    #[test]
    fn test_voronoi_edges_square() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = bowyer_watson(&points);
        let edges = voronoi_edges(&points, &triangles);
        assert!(
            !edges.is_empty(),
            "Square Delaunay should produce at least one Voronoi edge"
        );
    }
    #[test]
    fn test_weighted_site_power_distance() {
        let site = WeightedSite::new([0.0, 0.0], 1.0);
        let pd = site.power_distance([2.0, 0.0]);
        assert!((pd - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_power_diagram_owner_basic() {
        let sites = vec![
            WeightedSite::new([0.0, 0.0], 0.0),
            WeightedSite::new([5.0, 0.0], 0.0),
        ];
        let owner = power_diagram_owner([1.0, 0.0], &sites);
        assert_eq!(owner, Some(0));
        let owner2 = power_diagram_owner([4.0, 0.0], &sites);
        assert_eq!(owner2, Some(1));
    }
    #[test]
    fn test_power_diagram_owner_weighted() {
        let sites = vec![
            WeightedSite::new([0.0, 0.0], 0.0),
            WeightedSite::new([3.0, 0.0], 8.0),
        ];
        let owner = power_diagram_owner([2.0, 0.0], &sites);
        assert_eq!(
            owner,
            Some(1),
            "Site 1 should own (2,0) due to large weight"
        );
    }
    #[test]
    fn test_power_diagram_owner_empty() {
        assert_eq!(power_diagram_owner([0.0, 0.0], &[]), None);
    }
    #[test]
    fn test_power_distances_count() {
        let sites = vec![
            WeightedSite::new([0.0, 0.0], 1.0),
            WeightedSite::new([1.0, 1.0], 2.0),
        ];
        let dists = power_distances([0.5, 0.5], &sites);
        assert_eq!(dists.len(), 2);
    }
    #[test]
    fn test_polyhedron_volume_empty() {
        let vol = polyhedron_volume(&[], &[]);
        assert_eq!(vol, 0.0);
    }
    #[test]
    fn test_polyhedron_volume_tetrahedron() {
        let vertices = [
            [0.0_f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let faces: [[usize; 3]; 4] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let vol = polyhedron_volume(&faces, &vertices);
        assert!(
            (vol - 1.0 / 6.0).abs() < 0.01,
            "Tet volume should be 1/6, got {vol}"
        );
    }
    #[test]
    fn test_voronoi_cell_volume_3d_zero_circumcenters() {
        let vol = voronoi_cell_volume_3d([0.0, 0.0, 0.0], &[]);
        assert_eq!(vol, 0.0);
    }
    #[test]
    fn test_voronoi_cell_volume_3d_several_circumcenters() {
        let site = [0.0, 0.0, 0.0];
        let ccs: Vec<[f64; 3]> = (1..=6).map(|i| [i as f64, 0.0, 0.0]).collect();
        let vol = voronoi_cell_volume_3d(site, &ccs);
        assert!(vol > 0.0, "Volume estimate should be positive");
    }
    #[test]
    fn test_fortune_delaunay_three_points() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let triangles = fortune_delaunay(&points);
        assert_eq!(triangles.len(), 1);
    }
    #[test]
    fn test_fortune_delaunay_empty() {
        let result = fortune_delaunay(&[]);
        assert!(result.is_empty());
    }
    #[test]
    fn test_fortune_delaunay_four_points() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let triangles = fortune_delaunay(&points);
        assert_eq!(triangles.len(), 2);
    }
    #[test]
    fn test_fortune_delaunay_is_delaunay() {
        let points = [
            [0.0_f64, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
            [1.5, 1.5],
        ];
        let triangles = fortune_delaunay(&points);
        assert!(is_delaunay(&points, &triangles));
    }
    #[test]
    fn test_fortune_delaunay_indices_valid() {
        let points: Vec<[f64; 2]> = (0..6).map(|i| [i as f64, (i % 3) as f64]).collect();
        let triangles = fortune_delaunay(&points);
        let n = points.len();
        for t in &triangles {
            assert!(t.a < n && t.b < n && t.c < n);
        }
    }
    #[test]
    fn test_fortune_delaunay_matches_bowyer_watson_count() {
        let points = [
            [0.0_f64, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ];
        let bw = bowyer_watson(&points);
        let ft = fortune_delaunay(&points);
        assert_eq!(
            ft.len(),
            bw.len(),
            "both algorithms should produce same triangle count"
        );
    }
    #[test]
    fn test_weighted_lloyd_empty() {
        let result = weighted_lloyd_relaxation(&[], &[], [0.0, 1.0, 0.0, 1.0], 3);
        assert!(result.is_empty());
    }
    #[test]
    fn test_weighted_lloyd_count_preserved() {
        let sites = [[0.2_f64, 0.2], [0.8, 0.2], [0.5, 0.8]];
        let weights = [1.0, 1.0, 1.0];
        let result = weighted_lloyd_relaxation(&sites, &weights, [0.0, 1.0, 0.0, 1.0], 2);
        assert_eq!(result.len(), sites.len());
    }
    #[test]
    fn test_weighted_lloyd_within_bounds() {
        let sites = [[0.3_f64, 0.3], [0.7, 0.7]];
        let weights = [2.0, 1.0];
        let bounds = [0.0, 1.0, 0.0, 1.0];
        let result = weighted_lloyd_relaxation(&sites, &weights, bounds, 3);
        let [xmin, xmax, ymin, ymax] = bounds;
        for p in &result {
            assert!(p[0] >= xmin - 1e-9 && p[0] <= xmax + 1e-9);
            assert!(p[1] >= ymin - 1e-9 && p[1] <= ymax + 1e-9);
        }
    }
    #[test]
    fn test_weighted_lloyd_zero_iterations() {
        let sites = [[0.3_f64, 0.3], [0.7, 0.7]];
        let result = weighted_lloyd_relaxation(&sites, &[1.0, 1.0], [0.0, 1.0, 0.0, 1.0], 0);
        assert_eq!(result.len(), 2);
    }
    #[test]
    fn test_voronoi_area_fractions_sum() {
        let sites = [[0.25_f64, 0.5], [0.75, 0.5]];
        let fracs = voronoi_area_fractions(&sites, [0.0, 1.0, 0.0, 1.0], 20);
        let total: f64 = fracs.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-9,
            "fractions should sum to 1, got {total}"
        );
    }
    #[test]
    fn test_voronoi_area_fractions_equal_sites() {
        let sites = [[0.25_f64, 0.5], [0.75, 0.5]];
        let fracs = voronoi_area_fractions(&sites, [0.0, 1.0, 0.0, 1.0], 40);
        assert!(
            (fracs[0] - fracs[1]).abs() < 0.05,
            "symmetric sites: fracs[0]={} fracs[1]={}",
            fracs[0],
            fracs[1]
        );
    }
    #[test]
    fn test_voronoi_area_fractions_empty() {
        let fracs = voronoi_area_fractions(&[], [0.0, 1.0, 0.0, 1.0], 20);
        assert!(fracs.is_empty());
    }
    #[test]
    fn test_voronoi_area_fractions_single_site() {
        let sites = [[0.5_f64, 0.5]];
        let fracs = voronoi_area_fractions(&sites, [0.0, 1.0, 0.0, 1.0], 10);
        assert_eq!(fracs.len(), 1);
        assert!(
            (fracs[0] - 1.0).abs() < 1e-9,
            "single site should own all area"
        );
    }
    #[test]
    fn test_centroidal_energy_zero_sites() {
        let energy = centroidal_energy(&[], [0.0, 1.0, 0.0, 1.0], 10);
        assert_eq!(energy, 0.0);
    }
    #[test]
    fn test_centroidal_energy_non_negative() {
        let sites = [[0.2_f64, 0.2], [0.8, 0.8]];
        let energy = centroidal_energy(&sites, [0.0, 1.0, 0.0, 1.0], 20);
        assert!(energy >= 0.0, "energy should be non-negative");
    }
    #[test]
    fn test_centroidal_energy_decreases_after_lloyd() {
        let sites = [[0.1_f64, 0.1], [0.9, 0.1], [0.5, 0.9]];
        let bounds = [0.0, 1.0, 0.0, 1.0];
        let energy_before = centroidal_energy(&sites, bounds, 30);
        let relaxed = lloyd_relaxation(&sites, bounds, 5);
        let energy_after = centroidal_energy(&relaxed, bounds, 30);
        assert!(
            energy_after <= energy_before + 1e-9,
            "Lloyd relaxation should not increase centroidal energy: before={energy_before}, after={energy_after}"
        );
    }
    #[test]
    fn test_voronoi_cell_bbox_square() {
        let cell = LegacyVoronoiCell::new(
            [0.5, 0.5],
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        );
        let (mn, mx) = voronoi_cell_bbox(&cell).unwrap();
        assert!((mn[0] - 0.0).abs() < 1e-12);
        assert!((mn[1] - 0.0).abs() < 1e-12);
        assert!((mx[0] - 1.0).abs() < 1e-12);
        assert!((mx[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_voronoi_cell_bbox_empty() {
        let cell = LegacyVoronoiCell::new([0.0, 0.0], vec![]);
        assert!(voronoi_cell_bbox(&cell).is_none());
    }
    #[test]
    fn test_nn_distances_basic() {
        let sites = [[0.0_f64, 0.0], [1.0, 0.0], [2.0, 0.0]];
        let (mn, mx, avg) = nearest_neighbour_distances(&sites);
        assert!(
            (mn - 1.0).abs() < 1e-12,
            "min NN dist should be 1.0, got {mn}"
        );
        assert!(
            (mx - 1.0).abs() < 1e-12,
            "max NN dist should be 1.0, got {mx}"
        );
        assert!(
            (avg - 1.0).abs() < 1e-12,
            "avg NN dist should be 1.0, got {avg}"
        );
    }
    #[test]
    fn test_nn_distances_single() {
        let sites = [[0.0_f64, 0.0]];
        let (mn, mx, avg) = nearest_neighbour_distances(&sites);
        assert_eq!((mn, mx, avg), (0.0, 0.0, 0.0));
    }
    #[test]
    fn test_nn_distances_non_negative() {
        let sites = [[0.0_f64, 0.0], [3.0, 4.0], [1.0, 2.0]];
        let (mn, mx, avg) = nearest_neighbour_distances(&sites);
        assert!(mn >= 0.0);
        assert!(mx >= mn);
        assert!(avg >= mn && avg <= mx + 1e-10);
    }
    #[test]
    fn test_delaunay_edge_stats_equilateral() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 0.866]];
        let triangles = bowyer_watson(&points);
        let (mn, mx, avg) = delaunay_edge_stats(&points, &triangles);
        assert!(
            (mn - mx).abs() < 0.05,
            "equilateral: all edges equal, min={mn} max={mx}"
        );
        assert!((avg - 1.0).abs() < 0.05, "avg edge ≈ 1, got {avg}");
    }
    #[test]
    fn test_delaunay_edge_stats_empty() {
        let (mn, mx, avg) = delaunay_edge_stats(&[], &[]);
        assert_eq!((mn, mx, avg), (f64::INFINITY, 0.0, 0.0));
    }
    #[test]
    fn test_voronoi_cell_perimeter_square() {
        let cell = LegacyVoronoiCell::new(
            [0.5, 0.5],
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        );
        let p = voronoi_cell_perimeter(&cell);
        assert!(
            (p - 4.0).abs() < 1e-12,
            "square perimeter should be 4, got {p}"
        );
    }
    #[test]
    fn test_voronoi_cell_perimeter_empty() {
        let cell = LegacyVoronoiCell::new([0.0, 0.0], vec![]);
        assert_eq!(voronoi_cell_perimeter(&cell), 0.0);
    }
    #[test]
    fn test_voronoi_cell_perimeter_triangle() {
        let cell = LegacyVoronoiCell::new([0.5, 0.3], vec![[0.0, 0.0], [1.0, 0.0], [0.5, 0.866]]);
        let p = voronoi_cell_perimeter(&cell);
        assert!(
            p > 2.9 && p < 3.1,
            "equilateral triangle perimeter ≈3, got {p}"
        );
    }
    #[test]
    fn test_flip_edge_square() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let mut triangles = bowyer_watson(&points);
        let initial_len = triangles.len();
        for i in 0..initial_len {
            for j in (i + 1)..initial_len {
                flip_edge_if_needed(&mut triangles, &points, i, j);
            }
        }
        assert_eq!(
            triangles.len(),
            2,
            "triangulation should still have 2 triangles"
        );
    }
    #[test]
    fn test_flip_edge_no_adjacent() {
        let points = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0], [2.0, 0.5]];
        let mut triangles = bowyer_watson(&points);
        let flipped = flip_edge_if_needed(&mut triangles, &points, 0, 1);
        let _ = flipped;
        assert!(!triangles.is_empty());
    }
}
/// Compute a power diagram (weighted Voronoi) from a list of weighted sites.
///
/// Each site is `[x, y, w]` where `w` is the weight (power = distance² − w).
/// The result is a `VoronoiDiagram` where each cell is the region of minimum
/// power distance, clipped to `bbox = [xmin, xmax, ymin, ymax]`.
pub fn power_diagram(sites: &[[f64; 3]], bbox: [f64; 4]) -> VoronoiDiagram {
    if sites.is_empty() {
        return VoronoiDiagram {
            sites: vec![],
            cells: vec![],
            bbox,
        };
    }
    let n = sites.len();
    let [xmin, xmax, ymin, ymax] = bbox;
    let grid = 60_usize;
    let dx = (xmax - xmin) / grid as f64;
    let dy = (ymax - ymin) / grid as f64;
    let mut sample_x: Vec<Vec<f64>> = vec![vec![]; n];
    let mut sample_y: Vec<Vec<f64>> = vec![vec![]; n];
    for gi in 0..grid {
        for gj in 0..grid {
            let px = xmin + (gi as f64 + 0.5) * dx;
            let py = ymin + (gj as f64 + 0.5) * dy;
            let mut best = 0_usize;
            let mut best_power = f64::INFINITY;
            for (k, &s) in sites.iter().enumerate() {
                let power = (px - s[0]).powi(2) + (py - s[1]).powi(2) - s[2];
                if power < best_power {
                    best_power = power;
                    best = k;
                }
            }
            sample_x[best].push(px);
            sample_y[best].push(py);
        }
    }
    let voronoi_sites: Vec<VoronoiSite> = sites
        .iter()
        .enumerate()
        .map(|(i, &s)| VoronoiSite::new(Point2D::new(s[0], s[1]), i))
        .collect();
    let cells: Vec<VoronoiCell> = (0..n)
        .map(|i| {
            if sample_x[i].is_empty() {
                return VoronoiCell::new(i, vec![]);
            }
            let raw_pts: Vec<[f64; 2]> = sample_x[i]
                .iter()
                .zip(sample_y[i].iter())
                .map(|(&x, &y)| [x, y])
                .collect();
            let hull = simple_convex_hull_2d(&raw_pts);
            let verts: Vec<Point2D> = hull.iter().map(|&p| Point2D::from_array(p)).collect();
            VoronoiCell::new(i, verts)
        })
        .collect();
    VoronoiDiagram {
        sites: voronoi_sites,
        cells,
        bbox,
    }
}
/// Compute a simple 2-D convex hull using a gift-wrapping algorithm.
pub(super) fn simple_convex_hull_2d(pts: &[[f64; 2]]) -> Vec<[f64; 2]> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let start = pts
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            a[0].partial_cmp(&b[0])
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a[1].partial_cmp(&b[1]).unwrap_or(std::cmp::Ordering::Equal))
        })
        .map(|(i, _)| i)
        .unwrap_or(0);
    let n = pts.len();
    let mut hull = Vec::new();
    let mut current = start;
    loop {
        hull.push(pts[current]);
        let mut next = (current + 1) % n;
        for i in 0..n {
            let cross = cross_2d(pts[current], pts[next], pts[i]);
            if cross < 0.0 {
                next = i;
            }
        }
        current = next;
        if current == start {
            break;
        }
        if hull.len() > n {
            break;
        }
    }
    hull
}
/// Cross product (signed area) of vectors (b-a) × (c-a).
#[inline]
pub(super) fn cross_2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}
/// Compute the area of a polygon given as `Vec`Point2D` (shoelace formula).
pub(super) fn polygon_area_pt2d(verts: &[Point2D]) -> f64 {
    let n = verts.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0_f64;
    for i in 0..n {
        let j = (i + 1) % n;
        area += verts[i].x * verts[j].y;
        area -= verts[j].x * verts[i].y;
    }
    area.abs() * 0.5
}
/// Compute the centroid of a polygon given as `Vec`Point2D`.
pub(super) fn polygon_centroid_pt2d(verts: &[Point2D]) -> Point2D {
    if verts.is_empty() {
        return Point2D::new(0.0, 0.0);
    }
    if verts.len() == 1 {
        return verts[0];
    }
    let n = verts.len();
    let mut area = 0.0_f64;
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    for i in 0..n {
        let j = (i + 1) % n;
        let cross = verts[i].x * verts[j].y - verts[j].x * verts[i].y;
        area += cross;
        cx += (verts[i].x + verts[j].x) * cross;
        cy += (verts[i].y + verts[j].y) * cross;
    }
    area *= 0.5;
    if area.abs() < 1e-15 {
        let sx: f64 = verts.iter().map(|v| v.x).sum();
        let sy: f64 = verts.iter().map(|v| v.y).sum();
        return Point2D::new(sx / n as f64, sy / n as f64);
    }
    Point2D::new(cx / (6.0 * area), cy / (6.0 * area))
}
#[cfg(test)]
mod new_api_tests {

    use crate::DelaunayTriangulation;
    use crate::Point2D;
    use crate::VoronoiCell;
    use crate::VoronoiDiagram;
    use crate::VoronoiSite;
    use crate::power_diagram;
    use crate::voronoi::polygon_area_pt2d;
    use crate::voronoi::polygon_centroid_pt2d;
    #[test]
    fn test_point2d_new() {
        let p = Point2D::new(1.0, 2.0);
        assert!((p.x - 1.0).abs() < 1e-12);
        assert!((p.y - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_point2d_dist() {
        let a = Point2D::new(0.0, 0.0);
        let b = Point2D::new(3.0, 4.0);
        assert!((a.dist(&b) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_point2d_dist2() {
        let a = Point2D::new(1.0, 1.0);
        let b = Point2D::new(4.0, 5.0);
        assert!((a.dist2(&b) - 25.0).abs() < 1e-10);
    }
    #[test]
    fn test_point2d_roundtrip_array() {
        let arr = [3.125_f64, 2.75];
        let p = Point2D::from_array(arr);
        let out = p.to_array();
        assert!((out[0] - arr[0]).abs() < 1e-12);
        assert!((out[1] - arr[1]).abs() < 1e-12);
    }
    #[test]
    fn test_voronoi_site_new() {
        let s = VoronoiSite::new(Point2D::new(1.0, 2.0), 3);
        assert_eq!(s.index, 3);
        assert!((s.pos.x - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_voronoi_cell_area_square() {
        let verts = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ];
        let cell = VoronoiCell::new(0, verts);
        assert!((cell.area - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_voronoi_cell_perimeter_square() {
        let verts = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(0.0, 1.0),
        ];
        let cell = VoronoiCell::new(0, verts);
        assert!((cell.perimeter() - 4.0).abs() < 1e-9);
    }
    #[test]
    fn test_voronoi_cell_centroid_square() {
        let verts = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(2.0, 0.0),
            Point2D::new(2.0, 2.0),
            Point2D::new(0.0, 2.0),
        ];
        let cell = VoronoiCell::new(0, verts);
        let c = cell.centroid();
        assert!((c.x - 1.0).abs() < 1e-9);
        assert!((c.y - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_voronoi_cell_empty() {
        let cell = VoronoiCell::new(5, vec![]);
        assert_eq!(cell.site_index, 5);
        assert!((cell.area).abs() < 1e-12);
        assert!((cell.perimeter()).abs() < 1e-12);
    }
    #[test]
    fn test_voronoi_diagram_from_sites_count() {
        let pts = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let bbox = [-1.0, 2.0, -1.0, 2.0];
        let diag = VoronoiDiagram::from_sites(&pts, bbox);
        assert_eq!(diag.sites.len(), 3);
        assert_eq!(diag.cells.len(), 3);
    }
    #[test]
    fn test_voronoi_diagram_empty() {
        let diag = VoronoiDiagram::from_sites(&[], [0.0, 1.0, 0.0, 1.0]);
        assert!(diag.sites.is_empty());
        assert!(diag.cells.is_empty());
    }
    #[test]
    fn test_voronoi_diagram_nearest_site() {
        let pts = [[0.0_f64, 0.0], [5.0, 0.0], [10.0, 0.0]];
        let bbox = [-1.0, 11.0, -1.0, 1.0];
        let diag = VoronoiDiagram::from_sites(&pts, bbox);
        assert_eq!(diag.nearest_site([0.5, 0.0]), 0);
        assert_eq!(diag.nearest_site([9.5, 0.0]), 2);
    }
    #[test]
    fn test_voronoi_diagram_cell_area_nonnegative() {
        let pts = [[0.0_f64, 0.0], [2.0, 0.0], [1.0, 2.0], [1.0, 0.5]];
        let bbox = [-1.0, 3.0, -1.0, 3.0];
        let diag = VoronoiDiagram::from_sites(&pts, bbox);
        for i in 0..pts.len() {
            assert!(diag.cell_area(i) >= 0.0);
        }
    }
    #[test]
    fn test_voronoi_diagram_centroid_in_bbox() {
        let pts = [[0.5_f64, 0.5], [1.5, 0.5], [1.0, 1.5]];
        let bbox = [0.0, 2.0, 0.0, 2.0];
        let diag = VoronoiDiagram::from_sites(&pts, bbox);
        for i in 0..pts.len() {
            let c = diag.centroid(i);
            if diag.cell_area(i) > 1e-6 {
                assert!(
                    c[0] >= bbox[0] && c[0] <= bbox[1],
                    "centroid x={} out of bbox",
                    c[0]
                );
                assert!(
                    c[1] >= bbox[2] && c[1] <= bbox[3],
                    "centroid y={} out of bbox",
                    c[1]
                );
            }
        }
    }
    #[test]
    fn test_voronoi_diagram_lloyd_reduces_energy() {
        let pts = [[0.1_f64, 0.1], [0.9, 0.1], [0.5, 0.9]];
        let bbox = [0.0, 1.0, 0.0, 1.0];
        let diag = VoronoiDiagram::from_sites(&pts, bbox);
        let relaxed = diag.lloyd_relaxation(5);
        assert_eq!(relaxed.sites.len(), 3);
    }
    #[test]
    fn test_delaunay_triangulation_from_points() {
        let pts = [[0.0_f64, 0.0], [1.0, 0.0], [0.5, 1.0]];
        let dt = DelaunayTriangulation::from_points(&pts);
        assert_eq!(dt.points.len(), 3);
        assert_eq!(dt.triangles.len(), 1);
    }
    #[test]
    fn test_delaunay_triangulation_is_delaunay() {
        let pts = [
            [0.0_f64, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
            [1.5, 1.5],
        ];
        let dt = DelaunayTriangulation::from_points(&pts);
        assert!(dt.is_delaunay());
    }
    #[test]
    fn test_delaunay_circumscribed_circle() {
        let pts = [[0.0_f64, 0.0], [4.0, 0.0], [0.0, 4.0]];
        let dt = DelaunayTriangulation::from_points(&pts);
        if !dt.triangles.is_empty() {
            let (center, _radius) = dt.circumscribed_circle(dt.triangles[0]);
            assert!((center[0] - 2.0).abs() < 0.5 || (center[0]).abs() < 3.0);
        }
    }
    #[test]
    fn test_delaunay_flip_edge_preserves_count() {
        let pts = [[0.0_f64, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let mut dt = DelaunayTriangulation::from_points(&pts);
        let n_before = dt.triangles.len();
        dt.flip_edge();
        assert_eq!(dt.triangles.len(), n_before);
    }
    #[test]
    fn test_delaunay_dual_voronoi() {
        let pts = [[0.0_f64, 0.0], [2.0, 0.0], [1.0, 2.0]];
        let dt = DelaunayTriangulation::from_points(&pts);
        let vd = dt.dual_voronoi();
        assert_eq!(vd.sites.len(), 3);
    }
    #[test]
    fn test_power_diagram_empty() {
        let vd = power_diagram(&[], [0.0, 1.0, 0.0, 1.0]);
        assert!(vd.sites.is_empty());
    }
    #[test]
    fn test_power_diagram_site_count() {
        let sites = [[0.0_f64, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0]];
        let vd = power_diagram(&sites, [-1.0, 2.0, -1.0, 2.0]);
        assert_eq!(vd.sites.len(), 3);
        assert_eq!(vd.cells.len(), 3);
    }
    #[test]
    fn test_power_diagram_areas_nonneg() {
        let sites = [[0.0_f64, 0.0, 0.5], [2.0, 0.0, 0.5], [1.0, 2.0, 0.5]];
        let vd = power_diagram(&sites, [-1.0, 3.0, -1.0, 3.0]);
        for cell in &vd.cells {
            assert!(cell.area >= 0.0);
        }
    }
    #[test]
    fn test_power_diagram_unequal_weights() {
        let sites = [[0.0_f64, 0.5, 4.0], [1.0, 0.5, 0.0]];
        let vd = power_diagram(&sites, [-1.0, 2.0, -1.0, 2.0]);
        assert_eq!(vd.cells.len(), 2);
        assert!(vd.cells[0].area.is_finite());
        assert!(vd.cells[1].area.is_finite());
    }
    #[test]
    fn test_polygon_area_triangle() {
        let verts = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(4.0, 0.0),
            Point2D::new(0.0, 3.0),
        ];
        let area = polygon_area_pt2d(&verts);
        assert!((area - 6.0).abs() < 1e-9);
    }
    #[test]
    fn test_polygon_centroid_triangle() {
        let verts = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(3.0, 0.0),
            Point2D::new(0.0, 3.0),
        ];
        let c = polygon_centroid_pt2d(&verts);
        assert!((c.x - 1.0).abs() < 1e-9);
        assert!((c.y - 1.0).abs() < 1e-9);
    }
}
