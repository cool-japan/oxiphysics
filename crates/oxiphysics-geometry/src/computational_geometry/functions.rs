//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::{
    ArtGalleryResult, ConvexFace3D, ConvexHull3D, DelaunayTri, Line2D, Point2, VoronoiCell2D,
};

/// A point in 3D space represented as a plain array.
pub type Point3 = [f64; 3];
/// Subtract two 3D points.
pub fn sub3(a: Point3, b: Point3) -> Point3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Add two 3D vectors.
pub fn add3(a: Point3, b: Point3) -> Point3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3D vector.
pub fn scale3(a: Point3, t: f64) -> Point3 {
    [a[0] * t, a[1] * t, a[2] * t]
}
/// Dot product of two 3D vectors.
pub fn dot3(a: Point3, b: Point3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3D vectors.
pub fn cross3(a: Point3, b: Point3) -> Point3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Squared magnitude of a 3D vector.
pub fn mag3_sq(a: Point3) -> f64 {
    dot3(a, a)
}
/// Magnitude of a 3D vector.
pub fn mag3(a: Point3) -> f64 {
    mag3_sq(a).sqrt()
}
/// Normalize a 3D vector; returns the zero vector if the magnitude is too small.
pub fn normalize3(a: Point3) -> Point3 {
    let m = mag3(a);
    if m < 1e-15 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / m)
    }
}
/// Index of a half-edge in a `HalfEdgeMesh`.
pub type HalfEdgeId = usize;
/// Index of a vertex in a `HalfEdgeMesh`.
pub type VertexId = usize;
/// Index of a face in a `HalfEdgeMesh`.
pub type FaceId = usize;
/// Compute the signed volume of a tetrahedron formed by three points and the origin.
pub(super) fn signed_tet_volume(a: Point3, b: Point3, c: Point3) -> f64 {
    dot3(a, cross3(b, c)) / 6.0
}
/// Check if a point `p` is above (in the direction of the outward normal) a face defined
/// by vertices `a`, `b`, `c` (counter-clockwise order from outside).
pub(super) fn above_face(a: Point3, b: Point3, c: Point3, p: Point3) -> bool {
    let normal = cross3(sub3(b, a), sub3(c, a));
    dot3(normal, sub3(p, a)) > 1e-12
}
/// Build the initial tetrahedron from the first four non-coplanar points.
pub(super) fn initial_tetrahedron(pts: &[Point3]) -> Option<ConvexHull3D> {
    if pts.len() < 4 {
        return None;
    }
    let p0 = pts[0];
    let mut p1_idx = None;
    for i in 1..pts.len() {
        if mag3(sub3(pts[i], p0)) > 1e-10 {
            p1_idx = Some(i);
            break;
        }
    }
    let p1_idx = p1_idx?;
    let p1 = pts[p1_idx];
    let mut p2_idx = None;
    for i in (p1_idx + 1)..pts.len() {
        let v = cross3(sub3(p1, p0), sub3(pts[i], p0));
        if mag3(v) > 1e-10 {
            p2_idx = Some(i);
            break;
        }
    }
    let p2_idx = p2_idx?;
    let p2 = pts[p2_idx];
    let mut p3_idx = None;
    let normal = cross3(sub3(p1, p0), sub3(p2, p0));
    for i in (p2_idx + 1)..pts.len() {
        if dot3(normal, sub3(pts[i], p0)).abs() > 1e-10 {
            p3_idx = Some(i);
            break;
        }
    }
    let p3_idx = p3_idx?;
    let p3 = pts[p3_idx];
    let verts = vec![p0, p1, p2, p3];
    let cx = (p0[0] + p1[0] + p2[0] + p3[0]) / 4.0;
    let cy = (p0[1] + p1[1] + p2[1] + p3[1]) / 4.0;
    let cz = (p0[2] + p1[2] + p2[2] + p3[2]) / 4.0;
    let center: Point3 = [cx, cy, cz];
    let face_triples: [[usize; 3]; 4] = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
    let mut faces = Vec::new();
    for tri in &face_triples {
        let a = verts[tri[0]];
        let b = verts[tri[1]];
        let c = verts[tri[2]];
        let mut normal = normalize3(cross3(sub3(b, a), sub3(c, a)));
        if dot3(normal, sub3(a, center)) < 0.0 {
            normal = scale3(normal, -1.0);
            faces.push(ConvexFace3D {
                verts: [tri[0], tri[2], tri[1]],
                normal,
            });
        } else {
            faces.push(ConvexFace3D {
                verts: *tri,
                normal,
            });
        }
    }
    Some(ConvexHull3D {
        vertices: verts,
        faces,
    })
}
/// Compute the 3D convex hull of a set of points using the incremental algorithm.
///
/// Returns `None` if fewer than 4 non-coplanar points are provided.
///
/// # Algorithm
/// 1. Build an initial tetrahedron from 4 non-coplanar seed points.
/// 2. For each remaining point, find all faces visible from it.
/// 3. Remove visible faces, identify the horizon edge loop, and connect horizon
///    edges to the new point to form new faces.
///
/// Time complexity: O(n²) worst case; O(n log n) expected for random points.
pub fn convex_hull_3d(points: &[Point3]) -> Option<ConvexHull3D> {
    if points.len() < 4 {
        return None;
    }
    let mut hull = initial_tetrahedron(points)?;
    let skip_eps = 1e-10;
    let mut remaining: Vec<Point3> = Vec::new();
    'outer: for &p in points.iter() {
        for v in &hull.vertices {
            if mag3(sub3(p, *v)) < skip_eps {
                continue 'outer;
            }
        }
        remaining.push(p);
    }
    for &pt in &remaining {
        let visible: Vec<bool> = hull
            .faces
            .iter()
            .map(|f| {
                let a = hull.vertices[f.verts[0]];
                above_face(a, hull.vertices[f.verts[1]], hull.vertices[f.verts[2]], pt)
            })
            .collect();
        let any_visible = visible.iter().any(|&v| v);
        if !any_visible {
            continue;
        }
        let mut horizon: Vec<(usize, usize)> = Vec::new();
        for (fi, face) in hull.faces.iter().enumerate() {
            if !visible[fi] {
                continue;
            }
            let tri = face.verts;
            let edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
            for (ea, eb) in edges {
                let on_horizon = hull.faces.iter().enumerate().any(|(fj, g)| {
                    if visible[fj] {
                        return false;
                    }
                    let t = g.verts;
                    [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
                        .iter()
                        .any(|&(ga, gb)| ga == eb && gb == ea)
                });
                if on_horizon {
                    horizon.push((ea, eb));
                }
            }
        }
        let new_vid = hull.vertices.len();
        hull.vertices.push(pt);
        let n_v = hull.vertices.len() as f64;
        let cx: f64 = hull.vertices.iter().map(|v| v[0]).sum::<f64>() / n_v;
        let cy: f64 = hull.vertices.iter().map(|v| v[1]).sum::<f64>() / n_v;
        let cz: f64 = hull.vertices.iter().map(|v| v[2]).sum::<f64>() / n_v;
        let center = [cx, cy, cz];
        let new_faces: Vec<ConvexFace3D> = hull
            .faces
            .iter()
            .enumerate()
            .filter(|(fi, _)| !visible[*fi])
            .map(|(_, f)| f.clone())
            .collect();
        hull.faces = new_faces;
        for (ha, hb) in horizon {
            let a = hull.vertices[ha];
            let b = hull.vertices[hb];
            let mut normal = normalize3(cross3(sub3(b, a), sub3(pt, a)));
            let verts = if dot3(normal, sub3(a, center)) < 0.0 {
                normal = scale3(normal, -1.0);
                [ha, new_vid, hb]
            } else {
                [ha, hb, new_vid]
            };
            hull.faces.push(ConvexFace3D { verts, normal });
        }
    }
    Some(hull)
}
/// Clip a polygon against an infinite half-plane defined by the directed edge from
/// `edge_a` to `edge_b` (the inside is to the left).
///
/// Returns the clipped polygon vertices.
pub fn clip_polygon_by_edge(polygon: &[Point2], edge_a: Point2, edge_b: Point2) -> Vec<Point2> {
    if polygon.is_empty() {
        return vec![];
    }
    let n = polygon.len();
    let mut output = Vec::with_capacity(n + 1);
    for i in 0..n {
        let cur = polygon[i];
        let next = polygon[(i + 1) % n];
        let d_cur = Point2::cross2(edge_a, edge_b, cur);
        let d_next = Point2::cross2(edge_a, edge_b, next);
        if d_cur >= 0.0 {
            output.push(cur);
            if d_next < 0.0 {
                output.push(edge_intersect(cur, next, edge_a, edge_b));
            }
        } else {
            if d_next >= 0.0 {
                output.push(edge_intersect(cur, next, edge_a, edge_b));
            }
        }
    }
    output
}
/// Compute the intersection of segment (p, q) with the line through (a, b).
pub(super) fn edge_intersect(p: Point2, q: Point2, a: Point2, b: Point2) -> Point2 {
    let pq = q.sub(p);
    let ab = b.sub(a);
    let ap = a.sub(p);
    let denom = pq.cross(ab);
    if denom.abs() < 1e-15 {
        return p;
    }
    let t = ap.cross(ab) / denom;
    p.add(pq.scale(t))
}
/// Clip a subject polygon against a convex clipping polygon using the
/// Sutherland-Hodgman algorithm.
///
/// Both polygons must be given as counter-clockwise vertex lists.
/// Returns the vertices of the clipped (intersection) polygon, which is
/// convex if the clipping polygon is convex.
pub fn sutherland_hodgman(subject: &[Point2], clip: &[Point2]) -> Vec<Point2> {
    if subject.is_empty() || clip.is_empty() {
        return vec![];
    }
    let mut output = subject.to_vec();
    let n = clip.len();
    for i in 0..n {
        if output.is_empty() {
            break;
        }
        let a = clip[i];
        let b = clip[(i + 1) % n];
        output = clip_polygon_by_edge(&output, a, b);
    }
    output
}
/// Compute the area of a simple polygon (positive for CCW).
pub fn polygon_area(poly: &[Point2]) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut area = 0.0f64;
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly[i].x * poly[j].y;
        area -= poly[j].x * poly[i].y;
    }
    area / 2.0
}
/// Compute the Minkowski sum of two convex polygons P and Q.
///
/// Both polygons must be given as vertices in counter-clockwise order.
/// The result is the convex polygon P ⊕ Q.
///
/// # Algorithm
/// Rotating calipers / edge-angle merge:
/// 1. Rotate both polygons so the bottom-most (then left-most) vertex is first.
/// 2. Merge the edge direction sequences of P and Q in angular order.
/// 3. Walk the merged sequence to build the sum polygon.
pub fn minkowski_sum(p: &[Point2], q: &[Point2]) -> Vec<Point2> {
    if p.is_empty() || q.is_empty() {
        return vec![];
    }
    let start_p = bottom_vertex(p);
    let start_q = bottom_vertex(q);
    let n = p.len();
    let m = q.len();
    let p_rot: Vec<Point2> = (0..n).map(|i| p[(start_p + i) % n]).collect();
    let q_rot: Vec<Point2> = (0..m).map(|i| q[(start_q + i) % m]).collect();
    let edges_p: Vec<Point2> = (0..n).map(|i| p_rot[(i + 1) % n].sub(p_rot[i])).collect();
    let edges_q: Vec<Point2> = (0..m).map(|i| q_rot[(i + 1) % m].sub(q_rot[i])).collect();
    let mut result = Vec::with_capacity(n + m);
    let mut cur = p_rot[0].add(q_rot[0]);
    result.push(cur);
    let mut i = 0;
    let mut j = 0;
    while i < n || j < m {
        let ep = if i < n {
            edges_p[i]
        } else {
            Point2::new(1e18, 1e18)
        };
        let eq = if j < m {
            edges_q[j]
        } else {
            Point2::new(1e18, 1e18)
        };
        let cross = ep.cross(eq);
        let next_edge = if cross > 0.0 || j >= m {
            i += 1;
            ep
        } else if cross < 0.0 || i >= n {
            j += 1;
            eq
        } else {
            i += 1;
            j += 1;
            ep.add(eq)
        };
        if i <= n || j <= m {
            cur = cur.add(next_edge);
            result.push(cur);
        }
    }
    if result.len() > 1 {
        let last = *result.last().expect("collection should not be empty");
        let first = result[0];
        if last.dist_sq(first) < 1e-20 {
            result.pop();
        }
    }
    result
}
/// Find the index of the bottom-most (then left-most) vertex.
pub(super) fn bottom_vertex(poly: &[Point2]) -> usize {
    let mut best = 0;
    for i in 1..poly.len() {
        let b = poly[best];
        let c = poly[i];
        if c.y < b.y || (c.y == b.y && c.x < b.x) {
            best = i;
        }
    }
    best
}
/// Intersection point of two lines (if they are not parallel).
pub fn line_intersect_2d(l1: &Line2D, l2: &Line2D) -> Option<Point2> {
    let det = l1.a * l2.b - l2.a * l1.b;
    if det.abs() < 1e-15 {
        return None;
    }
    let x = (l1.c * l2.b - l2.c * l1.b) / det;
    let y = (l1.a * l2.c - l2.a * l1.c) / det;
    Some(Point2::new(x, y))
}
/// Test whether the open segment (p, q) is visible (does not cross any blocking edge).
///
/// Endpoint coincidences are permitted (two vertices of the same polygon can see
/// each other if they are adjacent, handled by the caller).
pub(super) fn segment_visible(p: Point2, q: Point2, blocking: &[(Point2, Point2)]) -> bool {
    for &(a, b) in blocking {
        if segments_properly_intersect(p, q, a, b) {
            return false;
        }
    }
    true
}
/// Test whether two open segments properly intersect (not at shared endpoints).
pub(super) fn segments_properly_intersect(p: Point2, q: Point2, a: Point2, b: Point2) -> bool {
    let eps = 1e-12;
    if p.dist_sq(a) < eps || p.dist_sq(b) < eps || q.dist_sq(a) < eps || q.dist_sq(b) < eps {
        return false;
    }
    let d1 = Point2::cross2(a, b, p);
    let d2 = Point2::cross2(a, b, q);
    let d3 = Point2::cross2(p, q, a);
    let d4 = Point2::cross2(p, q, b);
    if d1 * d2 < 0.0 && d3 * d4 < 0.0 {
        return true;
    }
    false
}
/// Greedy vertex-guard approximation for the art gallery problem.
///
/// Given a simple polygon (vertices in CCW order), iteratively places a guard
/// at the vertex that covers the maximum number of currently uncovered vertices,
/// until all vertices are covered or `max_guards` is reached.
///
/// This is a set-cover greedy approximation with ratio O(log n).
///
/// Returns a [`ArtGalleryResult`] describing the guards and coverage.
pub fn art_gallery_greedy(poly: &[Point2], max_guards: usize) -> ArtGalleryResult {
    let n = poly.len();
    if n == 0 {
        return ArtGalleryResult {
            guards: vec![],
            covered: vec![],
        };
    }
    let edges: Vec<(Point2, Point2)> = (0..n).map(|i| (poly[i], poly[(i + 1) % n])).collect();
    let vis: Vec<Vec<bool>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        true
                    } else {
                        segment_visible(poly[i], poly[j], &edges)
                    }
                })
                .collect()
        })
        .collect();
    let mut covered = vec![false; n];
    let mut guards = Vec::new();
    for _ in 0..max_guards {
        if covered.iter().all(|&c| c) {
            break;
        }
        let best = (0..n)
            .max_by_key(|&i| (0..n).filter(|&j| !covered[j] && vis[i][j]).count())
            .expect("operation should succeed");
        guards.push(best);
        for j in 0..n {
            if vis[best][j] {
                covered[j] = true;
            }
        }
    }
    ArtGalleryResult { guards, covered }
}
/// Check if every vertex of a polygon is covered by the given set of guards.
///
/// Returns `true` if full coverage is achieved.
pub fn check_full_coverage(poly: &[Point2], guards: &[usize]) -> bool {
    let n = poly.len();
    let edges: Vec<(Point2, Point2)> = (0..n).map(|i| (poly[i], poly[(i + 1) % n])).collect();
    for j in 0..n {
        let seen = guards
            .iter()
            .any(|&g| g == j || segment_visible(poly[g], poly[j], &edges));
        if !seen {
            return false;
        }
    }
    true
}
/// Circumcircle of three 2D points. Returns (center, radius²) or None if degenerate.
pub fn circumcircle_2d(pa: Point2, pb: Point2, pc: Point2) -> Option<(Point2, f64)> {
    let ax = pb.x - pa.x;
    let ay = pb.y - pa.y;
    let bx = pc.x - pa.x;
    let by = pc.y - pa.y;
    let d = 2.0 * (ax * by - ay * bx);
    if d.abs() < 1e-12 {
        return None;
    }
    let ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / d;
    let uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / d;
    let cx = pa.x + ux;
    let cy = pa.y + uy;
    let center = Point2::new(cx, cy);
    let r2 = center.dist_sq(pa);
    Some((center, r2))
}
/// Test if point `p` lies inside the circumcircle of triangle (a, b, c).
pub fn in_circumcircle_2d(pa: Point2, pb: Point2, pc: Point2, p: Point2) -> bool {
    match circumcircle_2d(pa, pb, pc) {
        Some((center, r2)) => center.dist_sq(p) < r2 - 1e-12,
        None => false,
    }
}
/// Bowyer-Watson incremental Delaunay triangulation for 2D points.
///
/// Returns the list of [`DelaunayTri`] triangles (indices into `points`).
/// Points are assumed to be in general position (no four co-circular).
pub fn delaunay_2d(points: &[Point2]) -> Vec<DelaunayTri> {
    let n = points.len();
    if n < 3 {
        return vec![];
    }
    let min_x = points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let max_x = points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let min_y = points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let max_y = points.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let dx = max_x - min_x;
    let dy = max_y - min_y;
    let delta = dx.max(dy) * 3.0 + 1.0;
    let mx = (min_x + max_x) / 2.0;
    let my = (min_y + max_y) / 2.0;
    let mut all_pts = points.to_vec();
    all_pts.push(Point2::new(mx - 20.0 * delta, my - delta));
    all_pts.push(Point2::new(mx, my + 20.0 * delta));
    all_pts.push(Point2::new(mx + 20.0 * delta, my - delta));
    let s0 = n;
    let s1 = n + 1;
    let s2 = n + 2;
    let mut triangles: Vec<DelaunayTri> = vec![DelaunayTri {
        a: s0,
        b: s1,
        c: s2,
    }];
    for i in 0..n {
        let pt = all_pts[i];
        let (bad, good): (Vec<_>, Vec<_>) = triangles
            .into_iter()
            .partition(|t| in_circumcircle_2d(all_pts[t.a], all_pts[t.b], all_pts[t.c], pt));
        let mut boundary: Vec<(usize, usize)> = Vec::new();
        for t in &bad {
            let edges = [(t.a, t.b), (t.b, t.c), (t.c, t.a)];
            for (ea, eb) in edges {
                let shared = bad.iter().any(|u| {
                    u != t
                        && ([(u.a, u.b), (u.b, u.c), (u.c, u.a)]
                            .iter()
                            .any(|&(ua, ub)| ua == eb && ub == ea))
                });
                if !shared {
                    boundary.push((ea, eb));
                }
            }
        }
        triangles = good;
        for (ea, eb) in boundary {
            triangles.push(DelaunayTri { a: ea, b: eb, c: i });
        }
    }
    triangles.retain(|t| t.a < n && t.b < n && t.c < n);
    triangles
}
/// Compute the dual Voronoi diagram from a Delaunay triangulation.
///
/// For each site point, collects all Delaunay triangles incident to that site
/// and returns their circumcenters as the Voronoi polygon vertices.
pub fn voronoi_from_delaunay(points: &[Point2], tris: &[DelaunayTri]) -> Vec<VoronoiCell2D> {
    let n = points.len();
    let mut cells: Vec<VoronoiCell2D> = (0..n)
        .map(|i| VoronoiCell2D {
            site: i,
            circumcenters: Vec::new(),
        })
        .collect();
    for t in tris {
        if let Some((center, _)) = circumcircle_2d(points[t.a], points[t.b], points[t.c]) {
            cells[t.a].circumcenters.push(center);
            cells[t.b].circumcenters.push(center);
            cells[t.c].circumcenters.push(center);
        }
    }
    for cell in &mut cells {
        let site = points[cell.site];
        cell.circumcenters.sort_by(|a, b| {
            let ta = (a.y - site.y).atan2(a.x - site.x);
            let tb = (b.y - site.y).atan2(b.x - site.x);
            ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    cells
}
/// Test if three 2D points are in counter-clockwise orientation.
pub fn ccw(a: Point2, b: Point2, c: Point2) -> bool {
    Point2::cross2(a, b, c) > 0.0
}
/// Test if three 2D points are collinear (within tolerance).
pub fn collinear(a: Point2, b: Point2, c: Point2) -> bool {
    Point2::cross2(a, b, c).abs() < 1e-10
}
/// Test if point `p` lies inside the triangle (a, b, c) using barycentric coordinates.
///
/// Returns `true` if strictly inside, `false` on boundary or outside.
pub fn point_in_triangle(p: Point2, a: Point2, b: Point2, c: Point2) -> bool {
    let d1 = Point2::cross2(p, a, b);
    let d2 = Point2::cross2(p, b, c);
    let d3 = Point2::cross2(p, c, a);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    !(has_neg && has_pos)
}
/// Test if point `p` lies inside a convex polygon (given in CCW order).
pub fn point_in_convex_polygon(p: Point2, poly: &[Point2]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        if Point2::cross2(a, b, p) < 0.0 {
            return false;
        }
    }
    true
}
/// Compute the 2D convex hull of a point set (Graham scan).
///
/// Returns vertices in CCW order.
pub fn convex_hull_2d(points: &[Point2]) -> Vec<Point2> {
    let mut pts = points.to_vec();
    if pts.len() < 3 {
        return pts;
    }
    pts.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .expect("operation should succeed")
            .then(a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
    });
    pts.dedup_by(|a, b| a.dist_sq(*b) < 1e-20);
    let n = pts.len();
    let mut hull: Vec<Point2> = Vec::with_capacity(2 * n);
    for &p in &pts {
        while hull.len() >= 2 {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            if Point2::cross2(a, b, p) <= 0.0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(p);
    }
    let lower_len = hull.len();
    for &p in pts.iter().rev() {
        while hull.len() > lower_len {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            if Point2::cross2(a, b, p) <= 0.0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(p);
    }
    hull.pop();
    hull
}
