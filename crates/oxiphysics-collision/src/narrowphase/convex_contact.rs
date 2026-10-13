// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Penetration queries between support-mapped convex shapes (GJK + EPA).
//!
//! [`penetration`] takes the support functions of two convex shapes and decides
//! whether they overlap. If they do it returns the minimum translation that
//! separates them: a unit normal, the penetration depth and one witness point on
//! each shape.
//!
//! The query runs on the Minkowski difference `A - B`:
//!
//! 1. **GJK** (distance form) refines the point of the current simplex that is
//!    closest to the origin until either the origin is enclosed (overlap), a
//!    separating plane is certified (no contact) or the closest feature cannot
//!    be improved any further (touching).
//! 2. The terminating simplex is completed to a non-degenerate tetrahedron with
//!    extra support queries, so that EPA can start from any overlap GJK reports.
//! 3. **EPA** grows a polytope inside `A - B` towards its boundary. The face
//!    closest to the origin gives the contact normal and depth. Every polytope
//!    vertex remembers the support points of `A` and `B` it was built from, so
//!    the witness points are the barycentric mix of those points on the winning
//!    face (not the shapes' origins or an axis-aligned guess).
//!
//! The normal points from shape A towards shape B and `point_a - point_b ==
//! normal * depth`: moving A by `-normal * depth` (or B by `normal * depth`)
//! separates the shapes. This is the orientation documented for
//! [`NarrowPhaseContact`].

use super::types::NarrowPhaseContact;

type Vec3 = [f64; 3];

/// Maximum number of GJK refinement steps.
const GJK_MAX_ITERATIONS: usize = 96;
/// GJK distance tolerance relative to the size of the shapes.
const GJK_TOLERANCE: f64 = 1e-10;
/// Maximum number of EPA expansion steps.
const EPA_MAX_ITERATIONS: usize = 256;
/// EPA convergence tolerance relative to the size of the shapes.
const EPA_TOLERANCE: f64 = 1e-9;
/// Upper bound on the number of EPA faces (guards against runaway growth).
const EPA_MAX_FACES: usize = 4096;
/// Minimum offset from the current affine hull, relative to the size of the
/// shapes, for a support point to count as a new affine dimension.
const INDEPENDENCE_TOLERANCE: f64 = 1e-8;
/// Relative tolerance below which a polytope face counts as visible.
const VISIBILITY_TOLERANCE: f64 = 1e-12;

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn neg(a: Vec3) -> Vec3 {
    [-a[0], -a[1], -a[2]]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// A vertex of the Minkowski difference `A - B` together with the points on the
/// two shapes it was built from.
#[derive(Debug, Clone, Copy)]
struct Vertex {
    /// `a - b`.
    w: Vec3,
    /// Support point on shape A.
    a: Vec3,
    /// Support point on shape B.
    b: Vec3,
}

/// The two support functions of a shape pair.
struct SupportPair<'a> {
    support_a: &'a dyn Fn(Vec3) -> Vec3,
    support_b: &'a dyn Fn(Vec3) -> Vec3,
}

impl SupportPair<'_> {
    /// Vertex of `A - B` that is extreme in direction `dir`.
    fn vertex(&self, dir: Vec3) -> Vertex {
        let a = (self.support_a)(dir);
        let b = (self.support_b)(neg(dir));
        Vertex { w: sub(a, b), a, b }
    }
}

/// Result of the GJK stage.
enum GjkOutcome {
    /// A separating plane was certified: the shapes do not touch.
    Separated,
    /// The shapes overlap or touch. `simplex` spans the origin (within the
    /// tolerance); `scale` is the largest Minkowski vertex norm seen (at least 1).
    Overlap { simplex: Vec<Vertex>, scale: f64 },
}

/// The point of a simplex closest to the origin and the simplex vertices that
/// carry it.
pub(super) struct Closest {
    /// The closest point.
    pub(super) point: Vec3,
    /// Indices of the simplex vertices whose convex hull contains `point`; only
    /// the first `len` entries are meaningful.
    pub(super) keep: [usize; 4],
    /// Number of vertices carrying `point` (1 to 4).
    pub(super) len: usize,
}

impl Closest {
    fn vertex(point: Vec3, i: usize) -> Self {
        Self {
            point,
            keep: [i, 0, 0, 0],
            len: 1,
        }
    }

    fn edge(point: Vec3, i: usize, j: usize) -> Self {
        Self {
            point,
            keep: [i, j, 0, 0],
            len: 2,
        }
    }

    fn face(point: Vec3, i: usize, j: usize, k: usize) -> Self {
        Self {
            point,
            keep: [i, j, k, 0],
            len: 3,
        }
    }

    fn dist_sq(&self) -> f64 {
        dot(self.point, self.point)
    }
}

/// Closest point of segment `ab` to the origin.
fn closest_on_segment(a: Vec3, b: Vec3, ia: usize, ib: usize) -> Closest {
    let ab = sub(b, a);
    let len_sq = dot(ab, ab);
    if len_sq <= 0.0 {
        return Closest::vertex(a, ia);
    }
    let t = -dot(a, ab) / len_sq;
    if t <= 0.0 {
        Closest::vertex(a, ia)
    } else if t >= 1.0 {
        Closest::vertex(b, ib)
    } else {
        Closest::edge(add(a, scale(ab, t)), ia, ib)
    }
}

/// Closest point of triangle `abc` to the origin (Ericson, "Real-Time Collision
/// Detection", 5.1.5), reporting the vertices of the closest feature.
fn closest_on_triangle(a: Vec3, b: Vec3, c: Vec3, ia: usize, ib: usize, ic: usize) -> Closest {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let normal = cross(ab, ac);
    // Collinear vertices have no interior; fall back to the three edges.
    if dot(normal, normal) <= 1e-28 * dot(ab, ab) * dot(ac, ac) {
        let candidates = [
            closest_on_segment(a, b, ia, ib),
            closest_on_segment(a, c, ia, ic),
            closest_on_segment(b, c, ib, ic),
        ];
        let mut best = Closest::vertex(a, ia);
        let mut best_sq = f64::INFINITY;
        for candidate in candidates {
            let d = candidate.dist_sq();
            if d < best_sq {
                best_sq = d;
                best = candidate;
            }
        }
        return best;
    }
    let ap = neg(a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return Closest::vertex(a, ia);
    }
    let bp = neg(b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return Closest::vertex(b, ib);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return Closest::edge(add(a, scale(ab, v)), ia, ib);
    }
    let cp = neg(c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return Closest::vertex(c, ic);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return Closest::edge(add(a, scale(ac, w)), ia, ic);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return Closest::edge(add(b, scale(sub(c, b), w)), ib, ic);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    Closest::face(add(add(a, scale(ab, v)), scale(ac, w)), ia, ib, ic)
}

/// Closest point of the solid tetrahedron `p` to the origin. The origin itself is
/// returned (with all four vertices kept) when it lies inside.
fn closest_on_tetrahedron(p: [Vec3; 4]) -> Closest {
    const FACES: [(usize, usize, usize, usize); 4] =
        [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 3, 1), (1, 2, 3, 0)];
    let mut longest_sq = 0.0_f64;
    for (i, pi) in p.iter().enumerate() {
        for pj in &p[i + 1..] {
            let e = sub(*pj, *pi);
            longest_sq = longest_sq.max(dot(e, e));
        }
    }
    // A flat tetrahedron has no inside: every face is then tested directly.
    let flat_tolerance = 1e-14 * longest_sq * longest_sq.sqrt();
    let mut best: Option<Closest> = None;
    let mut best_sq = f64::INFINITY;
    for (i, j, k, l) in FACES {
        let n = cross(sub(p[j], p[i]), sub(p[k], p[i]));
        let side_of_origin = -dot(n, p[i]);
        let side_of_vertex = dot(n, sub(p[l], p[i]));
        let origin_outside = if side_of_vertex.abs() <= flat_tolerance {
            true
        } else {
            side_of_origin * side_of_vertex < 0.0
        };
        if origin_outside {
            let candidate = closest_on_triangle(p[i], p[j], p[k], i, j, k);
            let d = candidate.dist_sq();
            if d < best_sq {
                best_sq = d;
                best = Some(candidate);
            }
        }
    }
    match best {
        Some(closest) => closest,
        None => Closest {
            point: [0.0; 3],
            keep: [0, 1, 2, 3],
            len: 4,
        },
    }
}

/// Closest point to the origin of the convex hull of `points` (the first four are
/// used), together with the points that carry it.
///
/// The origin itself is returned when it lies inside a tetrahedron. This is the
/// sub-algorithm shared by the GJK distance iterations of this module and the
/// GJK ray cast.
pub(super) fn closest_to_origin(points: &[Vec3]) -> Closest {
    match points {
        [] => Closest::vertex([0.0; 3], 0),
        [a] => Closest::vertex(*a, 0),
        [a, b] => closest_on_segment(*a, *b, 0, 1),
        [a, b, c] => closest_on_triangle(*a, *b, *c, 0, 1, 2),
        [a, b, c, d, ..] => closest_on_tetrahedron([*a, *b, *c, *d]),
    }
}

/// Closest point to the origin of the convex hull of `simplex` (1 to 4 vertices).
fn closest_on_simplex(simplex: &[Vertex]) -> Closest {
    let mut points = [[0.0; 3]; 4];
    let count = simplex.len().min(4);
    for (slot, vertex) in points.iter_mut().zip(simplex) {
        *slot = vertex.w;
    }
    closest_to_origin(&points[..count])
}

/// GJK in distance form on the Minkowski difference of the shape pair.
fn gjk(pair: &SupportPair<'_>, hint: Vec3) -> GjkOutcome {
    let mut simplex: Vec<Vertex> = vec![pair.vertex(hint)];
    let mut size = 1.0_f64.max(norm(simplex[0].w));
    for _ in 0..GJK_MAX_ITERATIONS {
        let closest = closest_on_simplex(&simplex);
        simplex = closest.keep[..closest.len]
            .iter()
            .map(|&i| simplex[i])
            .collect();
        let v = closest.point;
        let v_sq = dot(v, v);
        let tolerance = GJK_TOLERANCE * size;
        if v_sq <= tolerance * tolerance {
            return GjkOutcome::Overlap {
                simplex,
                scale: size,
            };
        }
        let v_len = v_sq.sqrt();
        let w = pair.vertex(neg(v));
        size = size.max(norm(w.w));
        let tolerance = GJK_TOLERANCE * size;
        // `w` is extreme in direction -v, so the plane through it with normal v
        // bounds A - B away from the origin: the distance is at least `lower`.
        let lower = dot(w.w, v) / v_len;
        if lower > tolerance {
            return GjkOutcome::Separated;
        }
        // The closest point can no longer be improved and lies within the
        // tolerance of the origin: the shapes touch.
        if v_len - lower <= tolerance {
            return GjkOutcome::Overlap {
                simplex,
                scale: size,
            };
        }
        simplex.push(w);
    }
    let closest = closest_on_simplex(&simplex);
    if norm(closest.point) <= 1e-6 * size {
        GjkOutcome::Overlap {
            simplex,
            scale: size,
        }
    } else {
        GjkOutcome::Separated
    }
}

/// Distance of `p` from the affine hull of `pts` (1 to 3 affinely independent
/// points).
fn distance_from_hull(pts: &[Vertex], p: Vec3) -> f64 {
    match pts {
        [] => 0.0,
        [a] => norm(sub(p, a.w)),
        [a, b] => {
            let e = sub(b.w, a.w);
            let r = sub(p, a.w);
            let along = dot(r, e) / dot(e, e);
            norm(sub(r, scale(e, along)))
        }
        [a, b, c, ..] => {
            let n = cross(sub(b.w, a.w), sub(c.w, a.w));
            dot(sub(p, a.w), n).abs() / norm(n)
        }
    }
}

/// Directions along which a support query can leave the affine hull of `pts`.
fn expansion_directions(pts: &[Vertex]) -> Vec<Vec3> {
    match pts {
        [] | [_] => vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
        [a, b] => {
            let e = sub(b.w, a.w);
            // The coordinate axis least aligned with the edge gives a stable
            // perpendicular.
            let mut axis = [0.0; 3];
            let mut least = f64::INFINITY;
            for k in 0..3 {
                if e[k].abs() < least {
                    least = e[k].abs();
                    axis = [0.0; 3];
                    axis[k] = 1.0;
                }
            }
            let first = cross(e, axis);
            let second = cross(e, first);
            vec![first, neg(first), second, neg(second)]
        }
        [a, b, c, ..] => {
            let n = cross(sub(b.w, a.w), sub(c.w, a.w));
            vec![n, neg(n)]
        }
    }
}

/// Complete `simplex` to a tetrahedron of affinely independent Minkowski
/// vertices, or `None` if `A - B` is flat (has no volume).
fn complete_tetrahedron(
    pair: &SupportPair<'_>,
    simplex: &[Vertex],
    size: f64,
) -> Option<[Vertex; 4]> {
    let independence = INDEPENDENCE_TOLERANCE * size;
    let mut pts: Vec<Vertex> = Vec::with_capacity(4);
    for v in simplex {
        if pts.len() == 4 {
            break;
        }
        if pts.is_empty() || distance_from_hull(&pts, v.w) > independence {
            pts.push(*v);
        }
    }
    if pts.is_empty() {
        pts.push(pair.vertex([1.0, 0.0, 0.0]));
    }
    while pts.len() < 4 {
        let mut best: Option<(f64, Vertex)> = None;
        for dir in expansion_directions(&pts) {
            let candidate = pair.vertex(dir);
            let offset = distance_from_hull(&pts, candidate.w);
            if offset > independence && best.is_none_or(|(best_offset, _)| offset > best_offset) {
                best = Some((offset, candidate));
            }
        }
        let (_, vertex) = best?;
        pts.push(vertex);
    }
    Some([pts[0], pts[1], pts[2], pts[3]])
}

/// A triangular face of the EPA polytope.
#[derive(Debug, Clone, Copy)]
struct Face {
    /// Vertex indices, counter-clockwise when seen from outside.
    v: [usize; 3],
    /// Outward unit normal.
    normal: Vec3,
    /// Signed distance of the face plane from the origin (negative while the
    /// origin is still outside the polytope).
    dist: f64,
}

/// The EPA polytope: a convex hull of Minkowski vertices that always contains
/// the interior point it was seeded with.
struct Polytope {
    verts: Vec<Vertex>,
    faces: Vec<Face>,
    interior: Vec3,
    /// Smallest face area² accepted as non-degenerate.
    area_sq_floor: f64,
    visibility: f64,
}

impl Polytope {
    fn new(tetra: [Vertex; 4], size: f64) -> Option<Self> {
        let sum = add(add(tetra[0].w, tetra[1].w), add(tetra[2].w, tetra[3].w));
        let mut polytope = Self {
            verts: tetra.to_vec(),
            faces: Vec::with_capacity(64),
            interior: scale(sum, 0.25),
            area_sq_floor: (1e-14 * size * size).powi(2),
            visibility: VISIBILITY_TOLERANCE * size,
        };
        for [i, j, k] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let face = polytope.make_face(i, j, k)?;
            polytope.faces.push(face);
        }
        Some(polytope)
    }

    /// Face through vertices `i`, `j`, `k`, wound and oriented outward.
    fn make_face(&self, i: usize, j: usize, k: usize) -> Option<Face> {
        let (pi, pj, pk) = (self.verts[i].w, self.verts[j].w, self.verts[k].w);
        let n = cross(sub(pj, pi), sub(pk, pi));
        let area_sq = dot(n, n);
        if area_sq.is_nan() || area_sq <= self.area_sq_floor {
            return None;
        }
        let mut normal = scale(n, 1.0 / area_sq.sqrt());
        let mut v = [i, j, k];
        if dot(normal, sub(pi, self.interior)) < 0.0 {
            normal = neg(normal);
            v = [i, k, j];
        }
        Some(Face {
            v,
            normal,
            dist: dot(normal, pi),
        })
    }

    /// The face with the smallest signed distance to the origin.
    fn closest_face(&self) -> Face {
        let mut best = self.faces[0];
        for face in &self.faces[1..] {
            if face.dist < best.dist {
                best = *face;
            }
        }
        best
    }

    /// Grow the polytope to include `w`. Returns `false`, leaving the polytope
    /// untouched, when `w` cannot be added without degenerate faces.
    fn add_vertex(&mut self, w: Vertex) -> bool {
        let visible: Vec<bool> = self
            .faces
            .iter()
            .map(|f| dot(f.normal, w.w) - f.dist > self.visibility)
            .collect();
        if !visible.iter().any(|&v| v) {
            return false;
        }
        // Horizon: directed edges of visible faces whose reverse edge belongs to
        // a face that stays.
        let mut horizon: Vec<(usize, usize)> = Vec::new();
        for (face, &is_visible) in self.faces.iter().zip(&visible) {
            if !is_visible {
                continue;
            }
            let [a, b, c] = face.v;
            for (start, end) in [(a, b), (b, c), (c, a)] {
                if let Some(pos) = horizon.iter().position(|&(s, e)| s == end && e == start) {
                    horizon.swap_remove(pos);
                } else {
                    horizon.push((start, end));
                }
            }
        }
        let closed = horizon.iter().all(|&(start, end)| {
            horizon.iter().filter(|&&(_, e)| e == start).count() == 1
                && horizon.iter().filter(|&&(s, _)| s == end).count() == 1
        });
        if horizon.len() < 3 || !closed {
            return false;
        }
        let new_index = self.verts.len();
        self.verts.push(w);
        let mut new_faces = Vec::with_capacity(horizon.len());
        for &(start, end) in &horizon {
            match self.make_face(start, end, new_index) {
                Some(face) => new_faces.push(face),
                None => {
                    self.verts.pop();
                    return false;
                }
            }
        }
        let mut index = 0;
        self.faces.retain(|_| {
            let keep = !visible[index];
            index += 1;
            keep
        });
        self.faces.extend(new_faces);
        true
    }

    /// Contact of the polytope face nearest to the origin.
    ///
    /// A flat facet of the polytope may be split into several coplanar
    /// triangles that all share the minimum distance; the contact is taken from
    /// the one that contains the origin's projection, so the witness points
    /// reproduce `normal * depth` exactly.
    fn nearest_contact(&self, tolerance: f64) -> NarrowPhaseContact {
        let nearest = self.closest_face();
        let mut best = nearest;
        let mut best_score = f64::NEG_INFINITY;
        for face in self
            .faces
            .iter()
            .filter(|f| f.dist <= nearest.dist + tolerance)
        {
            let score = self.containment(face);
            if score > best_score {
                best_score = score;
                best = *face;
            }
        }
        self.contact(&best)
    }

    /// Smallest barycentric coordinate of the origin's projection onto `face`;
    /// non-negative when the projection lies inside the triangle.
    fn containment(&self, face: &Face) -> f64 {
        let [a, b, c] = face.v.map(|i| self.verts[i].w);
        let (u, v, w) = barycentric_raw(scale(face.normal, face.dist), a, b, c);
        u.min(v).min(w)
    }

    /// Contact described by `face`.
    fn contact(&self, face: &Face) -> NarrowPhaseContact {
        let [a, b, c] = face.v.map(|i| self.verts[i]);
        // Barycentric coordinates of the origin's projection onto the face plane.
        let (u, v, w) = barycentric(scale(face.normal, face.dist), a.w, b.w, c.w);
        let mix = |pa: Vec3, pb: Vec3, pc: Vec3| add(add(scale(pa, u), scale(pb, v)), scale(pc, w));
        NarrowPhaseContact {
            normal: face.normal,
            depth: face.dist.max(0.0),
            point_a: mix(a.a, b.a, c.a),
            point_b: mix(a.b, b.b, c.b),
        }
    }
}

/// Barycentric coordinates `(u, v, w)` of `p` with respect to triangle `abc`
/// (`p = u a + v b + w c` for `p` in the triangle's plane). A degenerate
/// triangle yields `(NEG_INFINITY, ..)` so that it never counts as containing `p`.
fn barycentric_raw(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> (f64, f64, f64) {
    let v0 = sub(b, a);
    let v1 = sub(c, a);
    let v2 = sub(p, a);
    let d00 = dot(v0, v0);
    let d01 = dot(v0, v1);
    let d11 = dot(v1, v1);
    let d20 = dot(v2, v0);
    let d21 = dot(v2, v1);
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() <= f64::MIN_POSITIVE {
        return (f64::NEG_INFINITY, 0.0, 0.0);
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    (1.0 - v - w, v, w)
}

/// Barycentric coordinates of `p` in triangle `abc`, clamped to the triangle and
/// renormalised so that they sum to one.
fn barycentric(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> (f64, f64, f64) {
    let (u, v, w) = barycentric_raw(p, a, b, c);
    if u == f64::NEG_INFINITY {
        return (1.0, 0.0, 0.0);
    }
    let (u, v, w) = (u.max(0.0), v.max(0.0), w.max(0.0));
    let total = u + v + w;
    if total <= 0.0 {
        return (1.0, 0.0, 0.0);
    }
    (u / total, v / total, w / total)
}

/// Expanding Polytope Algorithm from a tetrahedron of Minkowski vertices.
///
/// The face with the smallest *signed* distance is expanded first, so a seed
/// tetrahedron that does not yet enclose the origin (touching shapes, or a
/// completed lower-dimensional GJK simplex) is grown until it does.
fn epa(pair: &SupportPair<'_>, tetra: [Vertex; 4], size: f64) -> Option<NarrowPhaseContact> {
    let mut polytope = Polytope::new(tetra, size)?;
    let tolerance = EPA_TOLERANCE * size;
    for _ in 0..EPA_MAX_ITERATIONS {
        if polytope.faces.len() > EPA_MAX_FACES {
            break;
        }
        let face = polytope.closest_face();
        let w = pair.vertex(face.normal);
        let improvement = dot(w.w, face.normal) - face.dist;
        if improvement <= tolerance || !polytope.add_vertex(w) {
            break;
        }
    }
    Some(polytope.nearest_contact(tolerance))
}

/// Contact for a pair whose Minkowski difference is flat (zero volume), e.g. two
/// coplanar polygons. Such shapes touch without penetrating: the depth is zero
/// and the normal is the plane normal, oriented from A towards B.
fn flat_contact(simplex: &[Vertex]) -> NarrowPhaseContact {
    let count = simplex.len().max(1) as f64;
    let mut centre_a = [0.0; 3];
    let mut centre_b = [0.0; 3];
    for v in simplex {
        centre_a = add(centre_a, v.a);
        centre_b = add(centre_b, v.b);
    }
    let centre_a = scale(centre_a, 1.0 / count);
    let centre_b = scale(centre_b, 1.0 / count);
    let mut normal = match simplex {
        [a, b, c, ..] => cross(sub(b.w, a.w), sub(c.w, a.w)),
        [a, b] => {
            let e = sub(b.w, a.w);
            let axis = if e[0].abs() <= e[1].abs() && e[0].abs() <= e[2].abs() {
                [1.0, 0.0, 0.0]
            } else if e[1].abs() <= e[2].abs() {
                [0.0, 1.0, 0.0]
            } else {
                [0.0, 0.0, 1.0]
            };
            cross(e, axis)
        }
        _ => [0.0, 1.0, 0.0],
    };
    let length = norm(normal);
    normal = if length > 0.0 {
        scale(normal, 1.0 / length)
    } else {
        [0.0, 1.0, 0.0]
    };
    if dot(normal, sub(centre_b, centre_a)) < 0.0 {
        normal = neg(normal);
    }
    NarrowPhaseContact {
        normal,
        depth: 0.0,
        point_a: centre_a,
        point_b: centre_b,
    }
}

/// Penetration query between two convex shapes given by their support functions
/// (world-space point of the shape furthest along a direction).
///
/// `direction_hint` seeds the first GJK support query; pass the vector from
/// shape A's centre to shape B's centre. Returns `None` when the shapes are
/// separated; shapes that merely touch yield a contact of depth zero.
pub(super) fn penetration(
    support_a: &dyn Fn(Vec3) -> Vec3,
    support_b: &dyn Fn(Vec3) -> Vec3,
    direction_hint: Vec3,
) -> Option<NarrowPhaseContact> {
    let pair = SupportPair {
        support_a,
        support_b,
    };
    let hint = if direction_hint.iter().all(|c| c.is_finite()) && norm(direction_hint) > 1e-12 {
        direction_hint
    } else {
        [1.0, 0.0, 0.0]
    };
    match gjk(&pair, hint) {
        GjkOutcome::Separated => None,
        GjkOutcome::Overlap { simplex, scale } => {
            let contact = complete_tetrahedron(&pair, &simplex, scale)
                .and_then(|tetra| epa(&pair, tetra, scale))
                .unwrap_or_else(|| flat_contact(&simplex));
            Some(contact)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sphere(center: Vec3, radius: f64) -> impl Fn(Vec3) -> Vec3 {
        move |d| {
            let len = norm(d);
            if len <= 0.0 {
                return add(center, [0.0, radius, 0.0]);
            }
            add(center, scale(d, radius / len))
        }
    }

    fn hull(points: Vec<Vec3>) -> impl Fn(Vec3) -> Vec3 {
        move |d| {
            let mut best = points[0];
            let mut best_dot = dot(best, d);
            for p in &points[1..] {
                let value = dot(*p, d);
                if value > best_dot {
                    best_dot = value;
                    best = *p;
                }
            }
            best
        }
    }

    fn cube_points(center: Vec3, half: f64) -> Vec<Vec3> {
        let mut pts = Vec::new();
        for sx in [-1.0, 1.0] {
            for sy in [-1.0, 1.0] {
                for sz in [-1.0, 1.0] {
                    pts.push(add(center, [sx * half, sy * half, sz * half]));
                }
            }
        }
        pts
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn overlapping_spheres_match_the_closed_form() {
        let a = sphere([0.0, 0.0, 0.0], 1.0);
        let b = sphere([1.5, 0.0, 0.0], 1.0);
        let c = penetration(&a, &b, [1.5, 0.0, 0.0]).expect("spheres overlap");
        assert!(close(c.depth, 0.5, 1e-6), "depth {}", c.depth);
        assert!(close(c.normal[0], 1.0, 1e-4), "normal {:?}", c.normal);
        assert!(close(c.point_a[0], 1.0, 1e-4), "point_a {:?}", c.point_a);
        assert!(close(c.point_b[0], 0.5, 1e-4), "point_b {:?}", c.point_b);
    }

    #[test]
    fn separated_spheres_report_no_contact() {
        let a = sphere([0.0, 0.0, 0.0], 1.0);
        let b = sphere([2.001, 0.0, 0.0], 1.0);
        assert!(penetration(&a, &b, [1.0, 0.0, 0.0]).is_none());
    }

    #[test]
    fn touching_cubes_give_a_zero_depth_contact() {
        let a = hull(cube_points([0.0; 3], 1.0));
        let b = hull(cube_points([2.0, 0.0, 0.0], 1.0));
        let c = penetration(&a, &b, [2.0, 0.0, 0.0]).expect("touching counts as contact");
        assert!(c.depth.abs() < 1e-9, "depth {}", c.depth);
        assert!(close(c.normal[0], 1.0, 1e-9), "normal {:?}", c.normal);
    }

    #[test]
    fn axis_aligned_cubes_report_axis_depth_and_a_to_b_normal() {
        let a = hull(cube_points([0.0; 3], 1.0));
        let b = hull(cube_points([1.5, 0.2, -0.1], 1.0));
        let c = penetration(&a, &b, [1.5, 0.2, -0.1]).expect("cubes overlap");
        assert!(close(c.depth, 0.5, 1e-9), "depth {}", c.depth);
        assert!(close(c.normal[0], 1.0, 1e-9), "normal {:?}", c.normal);
        assert!(close(c.point_a[0], 1.0, 1e-9), "point_a {:?}", c.point_a);
        assert!(close(c.point_b[0], 0.5, 1e-9), "point_b {:?}", c.point_b);
        let gap = sub(c.point_a, c.point_b);
        assert!(norm(sub(gap, scale(c.normal, c.depth))) < 1e-9);
    }

    #[test]
    fn identical_concentric_cubes_resolve_through_a_face() {
        let a = hull(cube_points([0.0; 3], 1.0));
        let b = hull(cube_points([0.0; 3], 1.0));
        let c = penetration(&a, &b, [0.0; 3]).expect("cubes overlap");
        assert!(close(c.depth, 2.0, 1e-9), "depth {}", c.depth);
        let axis_aligned = c.normal.iter().filter(|n| n.abs() > 1.0 - 1e-9).count();
        assert_eq!(axis_aligned, 1, "normal {:?}", c.normal);
    }

    #[test]
    fn flat_pair_degrades_to_a_touching_contact() {
        let tri_a = hull(vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]);
        let tri_b = hull(vec![[0.2, 0.2, 0.0], [1.0, 0.2, 0.0], [0.2, 1.0, 0.0]]);
        let c = penetration(&tri_a, &tri_b, [0.0, 0.0, 1.0]).expect("coplanar triangles overlap");
        assert!(c.depth.abs() < 1e-12);
        assert!(
            close(c.normal[2].abs(), 1.0, 1e-12),
            "normal {:?}",
            c.normal
        );
    }

    /// Tiny deterministic generator for the randomised checks.
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
        }

        fn in_range(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.next_f64()
        }

        fn direction(&mut self) -> Vec3 {
            loop {
                let v = [
                    self.in_range(-1.0, 1.0),
                    self.in_range(-1.0, 1.0),
                    self.in_range(-1.0, 1.0),
                ];
                let len = norm(v);
                if len > 0.1 && len <= 1.0 {
                    return scale(v, 1.0 / len);
                }
            }
        }

        fn cloud(&mut self, center: Vec3, spread: f64, count: usize) -> Vec<Vec3> {
            (0..count)
                .map(|_| {
                    add(
                        center,
                        [
                            self.in_range(-spread, spread),
                            self.in_range(-spread, spread),
                            self.in_range(-spread, spread),
                        ],
                    )
                })
                .collect()
        }
    }

    fn overlaps(points_a: &[Vec3], points_b: &[Vec3], shift_a: Vec3) -> bool {
        let moved: Vec<Vec3> = points_a.iter().map(|p| add(*p, shift_a)).collect();
        let a = hull(moved);
        let b = hull(points_b.to_vec());
        penetration(&a, &b, [1.0, 0.3, 0.2]).is_some()
    }

    /// Smallest translation of A along `dir` (unit) that separates it from B.
    fn separation_along(points_a: &[Vec3], points_b: &[Vec3], dir: Vec3) -> f64 {
        let mut lo = 0.0;
        let mut hi = 8.0;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if overlaps(points_a, points_b, scale(dir, mid)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        hi
    }

    #[test]
    fn random_hull_pairs_return_the_minimum_translation() {
        let mut rng = Lcg(0x9E37_79B9_7F4A_7C15);
        let mut checked = 0;
        for _ in 0..60 {
            let points_a = rng.cloud([0.0; 3], 1.0, 9);
            let offset = [
                rng.in_range(-0.9, 0.9),
                rng.in_range(-0.9, 0.9),
                rng.in_range(-0.9, 0.9),
            ];
            let points_b = rng.cloud(offset, 1.0, 9);
            let a = hull(points_a.clone());
            let b = hull(points_b.clone());
            let Some(c) = penetration(&a, &b, offset) else {
                continue;
            };
            checked += 1;
            assert!(close(norm(c.normal), 1.0, 1e-9), "normal {:?}", c.normal);
            // Witness points lie on the supporting planes of the contact normal.
            let a_max = points_a
                .iter()
                .map(|p| dot(*p, c.normal))
                .fold(f64::NEG_INFINITY, f64::max);
            let b_min = points_b
                .iter()
                .map(|p| dot(*p, c.normal))
                .fold(f64::INFINITY, f64::min);
            assert!(close(dot(c.point_a, c.normal), a_max, 1e-7));
            assert!(close(dot(c.point_b, c.normal), b_min, 1e-7));
            assert!(close(a_max - b_min, c.depth, 1e-7));
            assert!(norm(sub(sub(c.point_a, c.point_b), scale(c.normal, c.depth))) < 1e-7);
            // Moving A along -normal by the depth separates the pair, and no
            // other direction separates them with a shorter move.
            let along_normal = separation_along(&points_a, &points_b, neg(c.normal));
            assert!(
                close(along_normal, c.depth, 1e-6),
                "separation {along_normal} vs depth {}",
                c.depth
            );
            for _ in 0..24 {
                let dir = rng.direction();
                let s = separation_along(&points_a, &points_b, dir);
                assert!(
                    s >= c.depth - 1e-6,
                    "direction {dir:?} separates after {s} < depth {}",
                    c.depth
                );
            }
        }
        assert!(checked >= 30, "only {checked} random pairs overlapped");
    }
}
