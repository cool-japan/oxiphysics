// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact swept sphere vs triangle time of impact.
//!
//! The sphere centre `c(t) = c0 + t * (c1 - c0)` touches a triangle when it
//! reaches the boundary of the triangle's Minkowski sum with the sphere, which
//! consists of
//!
//! 1. two offset triangles at distance `radius` above and below the face,
//! 2. three cylinders of `radius` around the edges, and
//! 3. three spheres of `radius` around the vertices
//!
//! (Fauerby, "Improved Collision Detection and Response", 2003; Ericson,
//! "Real-Time Collision Detection", 5.5.6). The time of impact is the earliest
//! time at which the centre enters any of these pieces; a sphere that already
//! touches the triangle at `t = 0` reports `t = 0`.

use super::point_cast::{
    closest_point_on_segment, cylinder_entry, earliest, sphere_entry, within_contact,
};
use super::{SweepResult, add, cross, dot, len, len_sq, scale, sub};

/// Tolerance on the barycentric coordinates of the face-hit test. A centre that
/// descends exactly onto an edge is also reported by the edge cylinder at the
/// same time; the tolerance only keeps rounding from opening a gap between the
/// face and the edge.
const BARYCENTRIC_SLACK: f64 = 1e-9;

/// Relative threshold below which the area of a triangle is treated as zero.
const DEGENERATE_AREA: f64 = 1e-12;

/// Result of a swept sphere vs triangle mesh query.
#[derive(Debug, Clone)]
pub struct TriangleMeshSweepResult {
    /// Time of impact.
    pub toi: f64,
    /// Triangle index that was hit.
    pub tri_index: usize,
    /// Contact normal.
    pub normal: [f64; 3],
    /// Contact point on the sphere.
    pub point_a: [f64; 3],
    /// Contact point on the triangle.
    pub point_b: [f64; 3],
}

/// Closest point on the segment edges of a triangle to `p`.
fn closest_point_on_edges(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    let mut best = v0;
    let mut best_dist_sq = f64::INFINITY;
    for (a, b) in [(v0, v1), (v1, v2), (v2, v0)] {
        let (candidate, _) = closest_point_on_segment(a, b, p);
        let dist_sq = len_sq(sub(p, candidate));
        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            best = candidate;
        }
    }
    best
}

/// Closest point on the triangle `v0 v1 v2` to `p` (Ericson, 5.1.5).
///
/// A triangle without area (collinear or coincident vertices) is handled as the
/// union of its edges.
pub(super) fn closest_point_on_triangle(
    v0: [f64; 3],
    v1: [f64; 3],
    v2: [f64; 3],
    p: [f64; 3],
) -> [f64; 3] {
    let ab = sub(v1, v0);
    let ac = sub(v2, v0);
    let ap = sub(p, v0);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return v0;
    }
    let bp = sub(p, v1);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return v1;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return add(v0, scale(ab, v));
    }
    let cp = sub(p, v2);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return v2;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return add(v0, scale(ac, w));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(v1, scale(sub(v2, v1), w));
    }
    let sum = va + vb + vc;
    if sum <= 0.0 || !sum.is_finite() {
        return closest_point_on_edges(v0, v1, v2, p);
    }
    let v = vb / sum;
    let w = vc / sum;
    add(add(v0, scale(ab, v)), scale(ac, w))
}

/// Unit normal of the triangle, or `None` when its area is (numerically) zero.
fn triangle_normal(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> Option<[f64; 3]> {
    let e1 = sub(v1, v0);
    let e2 = sub(v2, v0);
    let n = cross(e1, e2);
    let n_len = len(n);
    if n_len <= DEGENERATE_AREA * len(e1) * len(e2) {
        return None;
    }
    Some(scale(n, 1.0 / n_len))
}

/// True when the point `q` of the triangle's plane lies inside the triangle.
fn plane_point_in_triangle(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3], q: [f64; 3]) -> bool {
    let e1 = sub(v1, v0);
    let e2 = sub(v2, v0);
    let w = sub(q, v0);
    let d11 = dot(e1, e1);
    let d12 = dot(e1, e2);
    let d22 = dot(e2, e2);
    let dw1 = dot(w, e1);
    let dw2 = dot(w, e2);
    let denom = d11 * d22 - d12 * d12;
    if denom <= 0.0 {
        return false;
    }
    let v = (d22 * dw1 - d12 * dw2) / denom;
    let w_coord = (d11 * dw2 - d12 * dw1) / denom;
    let u = 1.0 - v - w_coord;
    u >= -BARYCENTRIC_SLACK && v >= -BARYCENTRIC_SLACK && w_coord >= -BARYCENTRIC_SLACK
}

/// First time in `[0, 1]` at which the centre reaches one of the two offset
/// triangles at signed distance `+-radius` from the face plane, with its
/// projection inside the triangle.
///
/// A centre that starts inside the slab `|distance| <= radius` never crosses an
/// offset plane; the contact is then made with an edge or a vertex.
fn face_entry(
    c0: [f64; 3],
    motion: [f64; 3],
    radius: f64,
    tri: [[f64; 3]; 3],
    unit_normal: [f64; 3],
) -> Option<f64> {
    let [v0, v1, v2] = tri;
    let side0 = dot(unit_normal, sub(c0, v0));
    let (sign, offset) = if side0 > radius {
        (1.0, radius)
    } else if side0 < -radius {
        (-1.0, -radius)
    } else {
        return None;
    };
    let rate = dot(unit_normal, motion);
    if sign * rate >= 0.0 {
        // Not approaching the offset plane the centre starts on the far side of.
        return None;
    }
    let t = (offset - side0) / rate;
    if t > 1.0 {
        return None;
    }
    let centre = add(c0, scale(motion, t));
    let on_plane = sub(centre, scale(unit_normal, offset));
    plane_point_in_triangle(v0, v1, v2, on_plane).then_some(t)
}

/// Contact normal (from the triangle towards the sphere) when the sphere centre
/// coincides with the triangle, so that the offset between them has no
/// direction: the face normal on the side the sphere came from, otherwise the
/// direction opposing the motion.
fn fallback_normal(
    c0: [f64; 3],
    motion: [f64; 3],
    v0: [f64; 3],
    unit_normal: Option<[f64; 3]>,
) -> [f64; 3] {
    if let Some(n) = unit_normal {
        let side = dot(n, sub(c0, v0));
        let along_motion = dot(n, motion) > 0.0;
        let faces_sphere = side > 0.0 || (side == 0.0 && !along_motion);
        return if faces_sphere { n } else { scale(n, -1.0) };
    }
    let motion_len = len(motion);
    if motion_len > 0.0 {
        scale(motion, -1.0 / motion_len)
    } else {
        [0.0, 1.0, 0.0]
    }
}

/// Assemble the contact data for the sphere centre at time `toi`.
fn contact_at(
    c0: [f64; 3],
    motion: [f64; 3],
    toi: f64,
    radius: f64,
    tri: [[f64; 3]; 3],
    unit_normal: Option<[f64; 3]>,
) -> SweepResult {
    let [v0, v1, v2] = tri;
    let centre = add(c0, scale(motion, toi));
    let on_triangle = closest_point_on_triangle(v0, v1, v2, centre);
    let offset = sub(centre, on_triangle);
    let dist = len(offset);
    let normal = if dist > 1e-12 * (1.0 + radius) {
        scale(offset, 1.0 / dist)
    } else {
        fallback_normal(c0, motion, v0, unit_normal)
    };
    SweepResult {
        toi,
        normal,
        point_a: sub(centre, scale(normal, radius)),
        point_b: on_triangle,
    }
}

/// Sweep a sphere against one triangle and return the first impact.
///
/// The sphere of `radius` moves with its centre from `c0` to `c1`. The triangle
/// `[v0, v1, v2]` is static and double sided. The result carries
/// `toi` in `[0, 1]`, the unit contact `normal` pointing from the triangle
/// towards the sphere, the contact point `point_a` on the sphere surface and
/// the contact point `point_b` on the triangle (a face, edge or vertex point).
///
/// The impact is exact for all contact features: a face hit, an edge hit (the
/// sphere's centre reaching the cylinder of `radius` around an edge) and a
/// vertex hit. A sphere that already touches or overlaps the triangle at
/// `c0` reports `toi = 0`. Returns `None` when the sphere never touches the
/// triangle during the step, and for non-finite input or a negative `radius`.
pub fn swept_sphere_triangle(
    c0: [f64; 3],
    c1: [f64; 3],
    radius: f64,
    triangle: &[[f64; 3]; 3],
) -> Option<SweepResult> {
    let [v0, v1, v2] = *triangle;
    let finite = c0
        .iter()
        .chain(c1.iter())
        .chain(triangle.iter().flatten())
        .all(|x| x.is_finite());
    if !finite || !radius.is_finite() || radius < 0.0 {
        return None;
    }
    let motion = sub(c1, c0);
    let unit_normal = triangle_normal(v0, v1, v2);

    let at_start = closest_point_on_triangle(v0, v1, v2, c0);
    if within_contact(len_sq(sub(c0, at_start)), radius) {
        return Some(contact_at(c0, motion, 0.0, radius, *triangle, unit_normal));
    }

    let mut best = None;
    if let Some(n) = unit_normal {
        best = earliest(best, face_entry(c0, motion, radius, *triangle, n));
    }
    for vertex in [v0, v1, v2] {
        best = earliest(best, sphere_entry(c0, motion, vertex, radius));
    }
    for (a, b) in [(v0, v1), (v1, v2), (v2, v0)] {
        best = earliest(best, cylinder_entry(c0, motion, a, b, radius));
    }
    let toi = best?;
    Some(contact_at(c0, motion, toi, radius, *triangle, unit_normal))
}

/// Axis-aligned bounds of a triangle.
fn triangle_bounds(tri: &[[f64; 3]; 3]) -> ([f64; 3], [f64; 3]) {
    let mut lo = tri[0];
    let mut hi = tri[0];
    for vertex in &tri[1..] {
        for axis in 0..3 {
            lo[axis] = lo[axis].min(vertex[axis]);
            hi[axis] = hi[axis].max(vertex[axis]);
        }
    }
    (lo, hi)
}

/// Sweep a sphere against a list of triangles and return the earliest impact.
///
/// Each triangle is `[v0, v1, v2]` as `[[f64; 3]; 3]` and is treated as double
/// sided. The sphere of `radius` moves with its centre from `c0` to `c1`. The
/// time of impact is exact for face, edge and vertex contacts (see
/// [`swept_sphere_triangle`]); the normal points from the triangle towards the
/// sphere. When several triangles are hit at the same time the one with the
/// lowest index is reported.
pub fn swept_sphere_triangle_mesh(
    c0: [f64; 3],
    c1: [f64; 3],
    radius: f64,
    triangles: &[[[f64; 3]; 3]],
) -> Option<TriangleMeshSweepResult> {
    // Bounds of the swept sphere, used to skip triangles that cannot be reached.
    let reach = radius.max(0.0) * (1.0 + 1e-9) + 1e-9;
    let mut sweep_lo = [0.0; 3];
    let mut sweep_hi = [0.0; 3];
    for axis in 0..3 {
        sweep_lo[axis] = c0[axis].min(c1[axis]) - reach;
        sweep_hi[axis] = c0[axis].max(c1[axis]) + reach;
    }

    let mut best: Option<TriangleMeshSweepResult> = None;
    for (tri_index, tri) in triangles.iter().enumerate() {
        let (lo, hi) = triangle_bounds(tri);
        let reachable = (0..3).all(|axis| lo[axis] <= sweep_hi[axis] && hi[axis] >= sweep_lo[axis]);
        if !reachable {
            continue;
        }
        let Some(hit) = swept_sphere_triangle(c0, c1, radius, tri) else {
            continue;
        };
        if best.as_ref().is_none_or(|current| hit.toi < current.toi) {
            best = Some(TriangleMeshSweepResult {
                toi: hit.toi,
                tri_index,
                normal: hit.normal,
                point_a: hit.point_a,
                point_b: hit.point_b,
            });
        }
    }
    best
}
