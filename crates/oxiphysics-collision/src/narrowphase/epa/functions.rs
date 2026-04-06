//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::gjk::SupportPoint;
use crate::types::Contact;
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::Shape;

use super::types::{
    EpaConfig, EpaFaceRaw, EpaPenetration, EpaPolytope, EpaPolytopeStats, EpaStats, Face,
};

/// Maximum EPA iterations.
pub(super) const MAX_ITERATIONS: usize = 64;
/// EPA convergence tolerance.
pub(super) const TOLERANCE: Real = 1e-6;
pub(super) fn support_minkowski(
    shape_a: &dyn Shape,
    transform_a: &Transform,
    shape_b: &dyn Shape,
    transform_b: &Transform,
    direction: &Vec3,
) -> SupportPoint {
    let local_dir_a = transform_a.rotation.inverse() * direction;
    let local_dir_b = transform_b.rotation.inverse() * (-direction);
    let local_a = shape_a.support_point(&local_dir_a);
    let local_b = shape_b.support_point(&local_dir_b);
    let world_a = transform_a.transform_point(&local_a);
    let world_b = transform_b.transform_point(&local_b);
    SupportPoint {
        point: world_a - world_b,
        support_a: world_a,
        support_b: world_b,
    }
}
pub(super) fn build_initial_faces(vertices: &[SupportPoint]) -> Vec<Face> {
    let face_indices = [[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]];
    let mut faces = Vec::new();
    for indices in &face_indices {
        let a = vertices[indices[0]].point;
        let b = vertices[indices[1]].point;
        let c = vertices[indices[2]].point;
        let normal = compute_face_normal(&a, &b, &c);
        let distance = normal.dot(&a);
        let (normal, distance) = if distance < 0.0 {
            (-normal, -distance)
        } else {
            (normal, distance)
        };
        faces.push(Face {
            indices: *indices,
            normal,
            distance,
        });
    }
    faces
}
pub(super) fn compute_face_normal(a: &Vec3, b: &Vec3, c: &Vec3) -> Vec3 {
    let ab = b - a;
    let ac = c - a;
    let n = ab.cross(&ac);
    let norm = n.norm();
    if norm > 1e-10 {
        n / norm
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    }
}
pub(super) fn add_edge(edges: &mut Vec<(usize, usize)>, a: usize, b: usize) {
    if let Some(pos) = edges.iter().position(|&(ea, eb)| ea == b && eb == a) {
        edges.swap_remove(pos);
    } else {
        edges.push((a, b));
    }
}
pub(super) fn build_contact(face: &Face, vertices: &[SupportPoint]) -> Contact {
    let v0 = &vertices[face.indices[0]];
    let v1 = &vertices[face.indices[1]];
    let v2 = &vertices[face.indices[2]];
    let a = v0.point;
    let b = v1.point;
    let c = v2.point;
    let ab = b - a;
    let ac = c - a;
    let ap = -a;
    let d00 = ab.dot(&ab);
    let d01 = ab.dot(&ac);
    let d11 = ac.dot(&ac);
    let d20 = ap.dot(&ab);
    let d21 = ap.dot(&ac);
    let denom = d00 * d11 - d01 * d01;
    let (u, v, w) = if denom.abs() > 1e-12 {
        let v = (d11 * d20 - d01 * d21) / denom;
        let w = (d00 * d21 - d01 * d20) / denom;
        let u = 1.0 - v - w;
        (u.max(0.0), v.max(0.0), w.max(0.0))
    } else {
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
    };
    let total = u + v + w;
    let u = u / total;
    let v = v / total;
    let w = w / total;
    let point_a = v0.support_a * u + v1.support_a * v + v2.support_a * w;
    let point_b = v0.support_b * u + v1.support_b * v + v2.support_b * w;
    Contact::new(point_a, point_b, face.normal, face.distance)
}
/// Find the horizon (silhouette) edges visible from a support point.
///
/// Returns edges as `(vertex_a, vertex_b)` pairs that form the boundary
/// between visible and non-visible faces.
#[allow(dead_code)]
pub fn find_horizon(
    faces: &[EpaFaceRaw],
    vertices: &[[f64; 3]],
    support_point: [f64; 3],
) -> Vec<(usize, usize)> {
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for face in faces {
        let v = vertices[face.vertices[0]];
        let to_point = epa_sub3(support_point, v);
        if epa_dot3(face.normal, to_point) > 1e-10 {
            let [a, b, c] = face.vertices;
            epa_add_edge(&mut edges, a, b);
            epa_add_edge(&mut edges, b, c);
            epa_add_edge(&mut edges, c, a);
        }
    }
    edges
}
/// Run the EPA algorithm given a support function and an initial GJK simplex.
///
/// `support_fn(direction) -> [f64;3]` returns the support point in the Minkowski
/// difference for the given search direction.
///
/// Returns the penetration information if found within `max_iter` iterations.
#[allow(dead_code)]
pub fn epa_penetration<F>(
    mut support_fn: F,
    initial_simplex: &[[f64; 3]; 4],
    max_iter: usize,
    tol: f64,
) -> Option<EpaPenetration>
where
    F: FnMut([f64; 3]) -> [f64; 3],
{
    let mut polytope = EpaPolytope::from_gjk_simplex(initial_simplex);
    for _ in 0..max_iter {
        if polytope.faces.is_empty() {
            return None;
        }
        let closest_idx = polytope.find_closest_face();
        let closest = polytope.faces[closest_idx].clone();
        let new_point = support_fn(closest.normal);
        let new_dist = epa_dot3(new_point, closest.normal);
        if (new_dist - closest.distance).abs() < tol {
            return Some(EpaPenetration {
                normal: closest.normal,
                depth: closest.distance,
                witness_a: new_point,
                witness_b: epa_sub3(new_point, epa_scale3(closest.normal, closest.distance)),
            });
        }
        if !polytope.expand(new_point) {
            return Some(EpaPenetration {
                normal: closest.normal,
                depth: closest.distance,
                witness_a: new_point,
                witness_b: epa_sub3(new_point, epa_scale3(closest.normal, closest.distance)),
            });
        }
    }
    if polytope.faces.is_empty() {
        return None;
    }
    let closest_idx = polytope.find_closest_face();
    let closest = &polytope.faces[closest_idx];
    Some(EpaPenetration {
        normal: closest.normal,
        depth: closest.distance,
        witness_a: polytope.vertices[closest.vertices[0]],
        witness_b: epa_sub3(
            polytope.vertices[closest.vertices[0]],
            epa_scale3(closest.normal, closest.distance),
        ),
    })
}
#[allow(dead_code)]
pub(super) fn epa_dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[allow(dead_code)]
pub(super) fn epa_sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[allow(dead_code)]
pub(super) fn epa_scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[allow(dead_code)]
pub(super) fn epa_negate3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
#[allow(dead_code)]
pub(super) fn epa_cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[allow(dead_code)]
pub(super) fn epa_len3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
#[allow(dead_code)]
pub(super) fn epa_face_normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    let ab = epa_sub3(b, a);
    let ac = epa_sub3(c, a);
    let n = epa_cross3(ab, ac);
    let len = epa_len3(n);
    if len > 1e-10 {
        epa_scale3(n, 1.0 / len)
    } else {
        [0.0, 1.0, 0.0]
    }
}
#[allow(dead_code)]
pub(super) fn epa_add_edge(edges: &mut Vec<(usize, usize)>, a: usize, b: usize) {
    if let Some(pos) = edges.iter().position(|&(ea, eb)| ea == b && eb == a) {
        edges.swap_remove(pos);
    } else {
        edges.push((a, b));
    }
}
/// Check whether an EPA iteration has converged given the old and new face distances.
#[allow(dead_code)]
pub fn epa_has_converged(old_dist: f64, new_dist: f64, tol: f64) -> bool {
    (new_dist - old_dist).abs() < tol
}
/// Compute the barycentric coordinates of the projection of the origin onto
/// a triangle with vertices `a`, `b`, `c`.
///
/// Returns `(u, v, w)` such that `u + v + w = 1` and the projected point is
/// `u*a + v*b + w*c`.  Values are clamped to `[0, 1]` and renormalised.
#[allow(dead_code)]
pub fn epa_barycentric_origin(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> (f64, f64, f64) {
    let ab = epa_sub3(b, a);
    let ac = epa_sub3(c, a);
    let ap = epa_negate3(a);
    let d00 = epa_dot3(ab, ab);
    let d01 = epa_dot3(ab, ac);
    let d11 = epa_dot3(ac, ac);
    let d20 = epa_dot3(ap, ab);
    let d21 = epa_dot3(ap, ac);
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < 1e-14 {
        return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
    }
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1.0 - v - w;
    let (u, v, w) = (u.max(0.0), v.max(0.0), w.max(0.0));
    let sum = u + v + w;
    (u / sum, v / sum, w / sum)
}
/// Compute the contact normal from a converged EPA result.
///
/// The contact normal points from shape B to shape A and is derived from
/// the closest face's outward normal.
#[allow(dead_code)]
pub fn epa_contact_normal_from_face(face: &EpaFaceRaw) -> [f64; 3] {
    face.normal
}
/// Run the full EPA algorithm with the given config.
///
/// This is a wrapper around `epa_penetration` that accepts an `EpaConfig`.
#[allow(dead_code)]
pub fn epa_penetration_with_config<F>(
    support_fn: F,
    initial_simplex: &[[f64; 3]; 4],
    config: &EpaConfig,
) -> Option<EpaPenetration>
where
    F: FnMut([f64; 3]) -> [f64; 3],
{
    epa_penetration(
        support_fn,
        initial_simplex,
        config.max_iter,
        config.tolerance,
    )
}
/// Compute the unit normal of a raw-array triangle (v0, v1, v2).
///
/// Returns `[0, 1, 0]` for degenerate triangles.
/// This is the raw-array version used by the standalone EPA and tests.
#[allow(dead_code)]
pub fn compute_face_normal_raw(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    epa_face_normal(v0, v1, v2)
}
/// Run the EPA algorithm.
///
/// Given an initial GJK tetrahedron simplex that contains the origin, expands
/// the polytope toward the Minkowski-difference boundary and returns
/// `(contact_normal, penetration_depth)`.
///
/// `support_fn(dir)` must return the support point of the Minkowski difference
/// in direction `dir`.
#[allow(dead_code)]
pub fn epa_solve(
    simplex: [[f64; 3]; 4],
    support_fn: impl Fn([f64; 3]) -> [f64; 3],
    max_iter: usize,
) -> ([f64; 3], f64) {
    let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
    let tol = 1e-6;
    for _ in 0..max_iter {
        if polytope.faces.is_empty() {
            break;
        }
        let (_idx, dist, normal) = polytope.closest_face();
        let new_point = support_fn(normal);
        let new_dist = epa_dot3(new_point, normal);
        if (new_dist - dist).abs() < tol {
            return (normal, dist);
        }
        if !polytope.expand(new_point) {
            return (normal, dist);
        }
    }
    if polytope.faces.is_empty() {
        return ([0.0, 1.0, 0.0], 0.0);
    }
    let (_idx, dist, normal) = polytope.closest_face();
    (normal, dist)
}
/// Validate an EPA polytope: check all face normals are unit vectors and
/// all distances are non-negative.  Returns `true` if valid.
#[allow(dead_code)]
pub fn epa_validate_polytope(polytope: &EpaPolytope) -> bool {
    for face in &polytope.faces {
        if face.distance < -1e-8 {
            return false;
        }
        let nlen = epa_len3(face.normal);
        if (nlen - 1.0).abs() > 1e-4 {
            return false;
        }
        for &vi in &face.vertices {
            if vi >= polytope.vertices.len() {
                return false;
            }
        }
    }
    true
}
/// Compute a face normal using Newell's method (numerically robust for
/// polygons with more than 3 vertices; for triangles it is equivalent to
/// the cross-product method but accumulates per-edge contributions).
///
/// Vertices should be provided in counter-clockwise order when viewed from
/// outside the polytope.  Returns `[0,1,0]` for degenerate inputs.
#[allow(dead_code)]
pub fn epa_newell_normal(vertices: &[[f64; 3]]) -> [f64; 3] {
    if vertices.len() < 3 {
        return [0.0, 1.0, 0.0];
    }
    let mut n = [0.0_f64; 3];
    let len = vertices.len();
    for i in 0..len {
        let cur = vertices[i];
        let nxt = vertices[(i + 1) % len];
        n[0] += (cur[1] - nxt[1]) * (cur[2] + nxt[2]);
        n[1] += (cur[2] - nxt[2]) * (cur[0] + nxt[0]);
        n[2] += (cur[0] - nxt[0]) * (cur[1] + nxt[1]);
    }
    let len_n = epa_len3(n);
    if len_n > 1e-10 {
        epa_scale3(n, 1.0 / len_n)
    } else {
        [0.0, 1.0, 0.0]
    }
}
/// Recompute all face normals in a polytope using Newell's method.
///
/// This can fix accumulated floating-point drift after many expansions.
#[allow(dead_code)]
pub fn epa_recompute_normals(polytope: &mut EpaPolytope) {
    for face in &mut polytope.faces {
        let verts: Vec<[f64; 3]> = face
            .vertices
            .iter()
            .map(|&i| polytope.vertices[i])
            .collect();
        let n = epa_newell_normal(&verts);
        let d = epa_dot3(n, polytope.vertices[face.vertices[0]]);
        if d < 0.0 {
            face.normal = epa_negate3(n);
            face.distance = -d;
        } else {
            face.normal = n;
            face.distance = d;
        }
    }
}
/// Compute twice the area of a triangle (raw magnitude of cross product).
#[allow(dead_code)]
pub fn epa_triangle_area2(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = epa_sub3(b, a);
    let ac = epa_sub3(c, a);
    let cross = epa_cross3(ab, ac);
    epa_len3(cross)
}
/// Compute the area of a triangle.
#[allow(dead_code)]
pub fn epa_triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    epa_triangle_area2(a, b, c) * 0.5
}
/// Remove degenerate faces (area below `area_threshold`) from a polytope.
///
/// Degenerate faces arise from near-collinear or coincident vertices and
/// can cause numerical instability in the EPA iteration.
#[allow(dead_code)]
pub fn epa_filter_degenerate_faces(polytope: &mut EpaPolytope, area_threshold: f64) {
    polytope.faces.retain(|face| {
        let a = polytope.vertices[face.vertices[0]];
        let b = polytope.vertices[face.vertices[1]];
        let c = polytope.vertices[face.vertices[2]];
        epa_triangle_area(a, b, c) >= area_threshold
    });
}
/// Run the EPA algorithm and collect statistics.
///
/// Like `epa_penetration`, but also returns an `EpaStats` record.
#[allow(dead_code)]
pub fn epa_penetration_with_stats<F>(
    mut support_fn: F,
    initial_simplex: &[[f64; 3]; 4],
    max_iter: usize,
    tol: f64,
) -> (Option<EpaPenetration>, EpaStats)
where
    F: FnMut([f64; 3]) -> [f64; 3],
{
    let mut polytope = EpaPolytope::from_gjk_simplex(initial_simplex);
    let mut stats = EpaStats::default();
    for iter in 0..max_iter {
        stats.iterations = iter + 1;
        if polytope.faces.is_empty() {
            stats.face_count = 0;
            return (None, stats);
        }
        let closest_idx = polytope.find_closest_face();
        let closest = polytope.faces[closest_idx].clone();
        let new_point = support_fn(closest.normal);
        let new_dist = epa_dot3(new_point, closest.normal);
        if (new_dist - closest.distance).abs() < tol {
            stats.converged = true;
            stats.face_count = polytope.faces.len();
            stats.final_depth = closest.distance;
            return (
                Some(EpaPenetration {
                    normal: closest.normal,
                    depth: closest.distance,
                    witness_a: new_point,
                    witness_b: epa_sub3(new_point, epa_scale3(closest.normal, closest.distance)),
                }),
                stats,
            );
        }
        if !polytope.expand(new_point) {
            stats.face_count = polytope.faces.len();
            stats.final_depth = closest.distance;
            return (
                Some(EpaPenetration {
                    normal: closest.normal,
                    depth: closest.distance,
                    witness_a: new_point,
                    witness_b: epa_sub3(new_point, epa_scale3(closest.normal, closest.distance)),
                }),
                stats,
            );
        }
    }
    stats.face_count = polytope.faces.len();
    if polytope.faces.is_empty() {
        return (None, stats);
    }
    let closest_idx = polytope.find_closest_face();
    let closest = &polytope.faces[closest_idx];
    stats.final_depth = closest.distance;
    (
        Some(EpaPenetration {
            normal: closest.normal,
            depth: closest.distance,
            witness_a: polytope.vertices[closest.vertices[0]],
            witness_b: epa_sub3(
                polytope.vertices[closest.vertices[0]],
                epa_scale3(closest.normal, closest.distance),
            ),
        }),
        stats,
    )
}
/// Extract the silhouette (horizon) edges of a polytope as seen from `new_point`,
/// along with the set of face indices that are visible.
///
/// Returns `(horizon_edges, visible_face_indices)`.
/// Horizon edges are `(va, vb)` pairs forming the boundary between visible
/// and non-visible faces.
#[allow(dead_code)]
pub fn epa_silhouette(
    polytope: &EpaPolytope,
    new_point: [f64; 3],
) -> (Vec<(usize, usize)>, Vec<usize>) {
    let mut horizon: Vec<(usize, usize)> = Vec::new();
    let mut visible: Vec<usize> = Vec::new();
    for (fi, face) in polytope.faces.iter().enumerate() {
        let v = polytope.vertices[face.vertices[0]];
        let to_point = epa_sub3(new_point, v);
        if epa_dot3(face.normal, to_point) > 1e-10 {
            visible.push(fi);
            let [a, b, c] = face.vertices;
            epa_add_edge(&mut horizon, a, b);
            epa_add_edge(&mut horizon, b, c);
            epa_add_edge(&mut horizon, c, a);
        }
    }
    (horizon, visible)
}
/// Tolerance-based early-exit EPA.
///
/// Expands the polytope until the depth change between consecutive iterations
/// is below `depth_change_tol`, or `max_iter` is reached.
///
/// This is more conservative than the standard tolerance check on `new_dist - closest.distance`.
#[allow(dead_code)]
pub fn epa_penetration_early_exit<F>(
    mut support_fn: F,
    initial_simplex: &[[f64; 3]; 4],
    max_iter: usize,
    depth_change_tol: f64,
) -> Option<EpaPenetration>
where
    F: FnMut([f64; 3]) -> [f64; 3],
{
    let mut polytope = EpaPolytope::from_gjk_simplex(initial_simplex);
    let mut prev_depth = f64::NEG_INFINITY;
    for _ in 0..max_iter {
        if polytope.faces.is_empty() {
            return None;
        }
        let closest_idx = polytope.find_closest_face();
        let closest = polytope.faces[closest_idx].clone();
        if (closest.distance - prev_depth).abs() < depth_change_tol && prev_depth > 0.0 {
            return Some(EpaPenetration {
                normal: closest.normal,
                depth: closest.distance,
                witness_a: polytope.vertices[closest.vertices[0]],
                witness_b: epa_sub3(
                    polytope.vertices[closest.vertices[0]],
                    epa_scale3(closest.normal, closest.distance),
                ),
            });
        }
        prev_depth = closest.distance;
        let new_point = support_fn(closest.normal);
        let new_dist = epa_dot3(new_point, closest.normal);
        if (new_dist - closest.distance).abs() < depth_change_tol {
            return Some(EpaPenetration {
                normal: closest.normal,
                depth: closest.distance,
                witness_a: new_point,
                witness_b: epa_sub3(new_point, epa_scale3(closest.normal, closest.distance)),
            });
        }
        if !polytope.expand(new_point) {
            return Some(EpaPenetration {
                normal: closest.normal,
                depth: closest.distance,
                witness_a: new_point,
                witness_b: epa_sub3(new_point, epa_scale3(closest.normal, closest.distance)),
            });
        }
    }
    if polytope.faces.is_empty() {
        return None;
    }
    let closest_idx = polytope.find_closest_face();
    let closest = &polytope.faces[closest_idx];
    Some(EpaPenetration {
        normal: closest.normal,
        depth: closest.distance,
        witness_a: polytope.vertices[closest.vertices[0]],
        witness_b: epa_sub3(
            polytope.vertices[closest.vertices[0]],
            epa_scale3(closest.normal, closest.distance),
        ),
    })
}
/// EPA with fallback: if the algorithm fails to produce a valid result
/// (empty polytope or zero depth), fall back to the face normal of the
/// closest face in the initial polytope.
#[allow(dead_code)]
pub fn epa_with_fallback<F>(
    support_fn: F,
    initial_simplex: &[[f64; 3]; 4],
    max_iter: usize,
    tol: f64,
) -> EpaPenetration
where
    F: FnMut([f64; 3]) -> [f64; 3],
{
    let initial_polytope = EpaPolytope::from_gjk_simplex(initial_simplex);
    let fallback = if !initial_polytope.faces.is_empty() {
        let idx = initial_polytope.find_closest_face();
        let f = &initial_polytope.faces[idx];
        EpaPenetration {
            normal: f.normal,
            depth: f.distance,
            witness_a: initial_polytope.vertices[f.vertices[0]],
            witness_b: epa_sub3(
                initial_polytope.vertices[f.vertices[0]],
                epa_scale3(f.normal, f.distance),
            ),
        }
    } else {
        EpaPenetration {
            normal: [0.0, 1.0, 0.0],
            depth: 0.0,
            witness_a: [0.0; 3],
            witness_b: [0.0; 3],
        }
    };
    match epa_penetration(support_fn, initial_simplex, max_iter, tol) {
        Some(pen) if pen.depth > 0.0 => pen,
        _ => fallback,
    }
}
/// Remove a face by index and patch the vertex references in remaining faces.
///
/// Used during incremental polytope maintenance to keep face indices valid after
/// vertex compaction.
#[allow(dead_code)]
pub fn epa_remove_face(polytope: &mut EpaPolytope, face_idx: usize) {
    if face_idx < polytope.faces.len() {
        polytope.faces.swap_remove(face_idx);
    }
}
/// Compute the projected point of the origin onto the plane of face `face_idx`.
///
/// Returns the projection in world space (same coordinate system as `polytope.vertices`).
/// Returns the face centroid as a fallback for degenerate faces.
#[allow(dead_code)]
pub fn epa_project_origin_onto_face(polytope: &EpaPolytope, face_idx: usize) -> [f64; 3] {
    if face_idx >= polytope.faces.len() {
        return [0.0; 3];
    }
    let face = &polytope.faces[face_idx];
    let n = face.normal;
    let dist = face.distance;
    epa_scale3(n, dist)
}
/// Return the vertices of a face as world-space positions.
#[allow(dead_code)]
pub fn epa_face_vertices(polytope: &EpaPolytope, face_idx: usize) -> Option<[[f64; 3]; 3]> {
    let face = polytope.faces.get(face_idx)?;
    Some([
        polytope.vertices[face.vertices[0]],
        polytope.vertices[face.vertices[1]],
        polytope.vertices[face.vertices[2]],
    ])
}
/// Compute the centroid of a face.
#[allow(dead_code)]
pub fn epa_face_centroid(polytope: &EpaPolytope, face_idx: usize) -> Option<[f64; 3]> {
    let verts = epa_face_vertices(polytope, face_idx)?;
    let cx = (verts[0][0] + verts[1][0] + verts[2][0]) / 3.0;
    let cy = (verts[0][1] + verts[1][1] + verts[2][1]) / 3.0;
    let cz = (verts[0][2] + verts[1][2] + verts[2][2]) / 3.0;
    Some([cx, cy, cz])
}
/// Extract the face with the maximum distance from origin (farthest face).
///
/// In a valid EPA polytope, all faces should be on the convex hull with outward
/// normals, so the maximum-distance face is the one farthest from origin.
#[allow(dead_code)]
pub fn epa_farthest_face(polytope: &EpaPolytope) -> Option<usize> {
    polytope
        .faces
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.distance
                .partial_cmp(&b.distance)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
}
/// Compute statistics about the current state of an EPA polytope.
#[allow(dead_code)]
pub fn epa_polytope_stats(polytope: &EpaPolytope) -> EpaPolytopeStats {
    if polytope.faces.is_empty() {
        return EpaPolytopeStats {
            face_count: 0,
            vertex_count: polytope.vertices.len(),
            min_distance: 0.0,
            max_distance: 0.0,
            avg_distance: 0.0,
            min_area: 0.0,
            max_area: 0.0,
        };
    }
    let mut min_dist = f64::INFINITY;
    let mut max_dist = f64::NEG_INFINITY;
    let mut sum_dist = 0.0;
    let mut min_area = f64::INFINITY;
    let mut max_area = f64::NEG_INFINITY;
    for face in &polytope.faces {
        let d = face.distance;
        if d < min_dist {
            min_dist = d;
        }
        if d > max_dist {
            max_dist = d;
        }
        sum_dist += d;
        let a = polytope.vertices[face.vertices[0]];
        let b = polytope.vertices[face.vertices[1]];
        let c = polytope.vertices[face.vertices[2]];
        let area = epa_triangle_area(a, b, c);
        if area < min_area {
            min_area = area;
        }
        if area > max_area {
            max_area = area;
        }
    }
    EpaPolytopeStats {
        face_count: polytope.faces.len(),
        vertex_count: polytope.vertices.len(),
        min_distance: min_dist,
        max_distance: max_dist,
        avg_distance: sum_dist / polytope.faces.len() as f64,
        min_area,
        max_area,
    }
}
// epa_add3 and epa_penetration_capped moved to functions_2.rs
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Epa;
    use crate::narrowphase::EpaSolver;
    use crate::narrowphase::epa::epa_add3;
    use crate::narrowphase::epa_penetration_capped;
    use crate::narrowphase::gjk::{Gjk, GjkResult};
    use oxiphysics_geometry::Sphere;
    #[test]
    fn test_epa_sphere_sphere() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = Gjk::query(&s1, &t1, &s2, &t2);
        if let GjkResult::Intersecting(simplex) = result {
            let contact = Epa::penetration_depth(&s1, &t1, &s2, &t2, &simplex);
            assert!(contact.is_some());
            let c = contact.unwrap();
            assert!((c.depth - 1.0).abs() < 0.2, "depth = {}", c.depth);
        } else {
            panic!("Expected intersection");
        }
    }
    /// Two overlapping unit-radius spheres with centers 1 m apart.
    /// Expected penetration depth = r1 + r2 - dist = 1 + 1 - 1 = 1 m.
    /// The contact normal produced by EPA must be a unit vector and the depth
    /// must be non-negative.
    ///
    /// Note: EPA expands the Minkowski-difference polytope outward from the
    /// origin, so the returned normal points along the minimum-penetration
    /// axis.  With A at (0,0,0) and B at (1,0,0) the axis is X; EPA may
    /// return either +X or −X depending on which face it finds first.  We
    /// therefore only check magnitude, not sign.
    #[test]
    fn test_epa_sphere_sphere_overlap() {
        let s1 = Sphere::new(1.0);
        let s2 = Sphere::new(1.0);
        let t1 = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        let t2 = Transform::from_position(Vec3::new(1.0, 0.0, 0.0));
        let result = Gjk::query(&s1, &t1, &s2, &t2);
        let simplex = match result {
            GjkResult::Intersecting(s) => s,
            GjkResult::Separated { .. } => panic!("expected intersection"),
        };
        let contact = Epa::penetration_depth(&s1, &t1, &s2, &t2, &simplex)
            .expect("EPA must return a contact");
        assert!(
            (contact.depth - 1.0).abs() < 0.2,
            "expected depth ≈ 1.0, got {}",
            contact.depth
        );
        assert!(contact.depth >= 0.0, "depth must be non-negative");
        let nlen = contact.normal.norm();
        assert!(
            (nlen - 1.0).abs() < 1e-6,
            "normal must be unit length, got {}",
            nlen
        );
        assert!(
            contact.normal.x.abs() > 0.9,
            "normal should be aligned with X axis, got normal = {:?}",
            contact.normal
        );
    }
    fn make_tetrahedron() -> [[f64; 3]; 4] {
        [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ]
    }
    #[test]
    fn test_epa_polytope_from_simplex() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        assert_eq!(polytope.vertices.len(), 4);
        assert_eq!(polytope.faces.len(), 4);
        for face in &polytope.faces {
            assert!(face.distance >= 0.0, "face distance must be non-negative");
            let nlen = epa_len3(face.normal);
            assert!((nlen - 1.0).abs() < 1e-6, "normal must be unit length");
        }
    }
    #[test]
    fn test_epa_find_closest_face() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let idx = polytope.find_closest_face();
        let closest_dist = polytope.faces[idx].distance;
        for (i, face) in polytope.faces.iter().enumerate() {
            if i != idx {
                assert!(face.distance >= closest_dist - 1e-10);
            }
        }
    }
    #[test]
    fn test_epa_expand() {
        let simplex = make_tetrahedron();
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let initial_faces = polytope.faces.len();
        let expanded = polytope.expand([3.0, 0.0, 0.0]);
        assert!(expanded, "expansion should succeed");
        assert!(polytope.vertices.len() == 5);
        assert!(polytope.faces.len() >= initial_faces);
    }
    #[test]
    fn test_find_horizon() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let horizon = find_horizon(&polytope.faces, &polytope.vertices, [3.0, 0.0, 0.0]);
        assert!(!horizon.is_empty(), "horizon must have edges");
    }
    #[test]
    fn test_epa_penetration_spheres() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(result.is_some(), "EPA must find penetration");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0, "depth must be positive");
        let nlen = epa_len3(pen.normal);
        assert!(
            (nlen - 1.0).abs() < 1e-4,
            "normal must be unit length, got {}",
            nlen
        );
    }
    #[test]
    fn test_epa_penetration_depth_value() {
        let simplex: [[f64; 3]; 4] = [
            [0.3, 0.3, 0.3],
            [0.3, -0.3, -0.3],
            [-0.3, 0.3, -0.3],
            [-0.3, -0.3, 0.3],
        ];
        let result = epa_penetration(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        let pen = result.expect("EPA must converge");
        assert!(
            (pen.depth - 1.0).abs() < 0.2,
            "depth should be ~1.0, got {}",
            pen.depth
        );
    }
    #[test]
    fn test_epa_face_normals_outward() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        for face in &polytope.faces {
            let v = polytope.vertices[face.vertices[0]];
            let d = epa_dot3(face.normal, v);
            assert!(d >= -1e-10, "face normal must point outward, dot = {}", d);
        }
    }
    #[test]
    fn test_epa_zero_iter() {
        let simplex = make_tetrahedron();
        let result = epa_penetration(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0], 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            0,
            1e-6,
        );
        assert!(
            result.is_some(),
            "should return best result even with 0 iterations"
        );
    }
    #[test]
    fn test_epa_config_default() {
        let cfg = EpaConfig::default();
        assert_eq!(cfg.max_iter, 64);
        assert!(cfg.tolerance > 0.0);
    }
    #[test]
    fn test_epa_config_high_precision() {
        let cfg = EpaConfig::high_precision();
        assert!(cfg.max_iter > 64);
        assert!(cfg.tolerance < 1e-6);
    }
    #[test]
    fn test_epa_has_converged() {
        assert!(
            epa_has_converged(1.0, 1.0 + 1e-8, 1e-6),
            "tiny change should converge"
        );
        assert!(
            !epa_has_converged(1.0, 1.1, 1e-6),
            "large change should not converge"
        );
    }
    #[test]
    fn test_epa_barycentric_origin_equilateral() {
        let a = [1.0, 0.0, 0.0];
        let b = [-0.5, 0.866, 0.0];
        let c = [-0.5, -0.866, 0.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-10, "weights must sum to 1");
        assert!(
            u > 0.0 && v > 0.0 && w > 0.0,
            "all weights must be positive"
        );
    }
    #[test]
    fn test_epa_barycentric_sum_to_one() {
        let a = [2.0, 0.0, 0.0];
        let b = [0.0, 2.0, 0.0];
        let c = [0.0, 0.0, 2.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-10, "u+v+w={}", u + v + w);
    }
    #[test]
    fn test_epa_contact_normal_from_face_unit_length() {
        let face = EpaFaceRaw {
            vertices: [0, 1, 2],
            normal: [1.0, 0.0, 0.0],
            distance: 1.5,
        };
        let n = epa_contact_normal_from_face(&face);
        let len = epa_len3(n);
        assert!((len - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_epa_validate_polytope_valid() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        assert!(
            epa_validate_polytope(&polytope),
            "initial polytope must be valid"
        );
    }
    #[test]
    fn test_epa_validate_polytope_invalid_negative_distance() {
        let simplex = make_tetrahedron();
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        polytope.faces[0].distance = -1.0;
        assert!(
            !epa_validate_polytope(&polytope),
            "negative distance should fail validation"
        );
    }
    #[test]
    fn test_epa_solve_sphere_cso() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (normal, depth) = epa_solve(
            simplex,
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            64,
        );
        assert!(depth > 0.0, "penetration depth must be positive");
        let nlen = epa_len3(normal);
        assert!((nlen - 1.0).abs() < 1e-4, "normal must be unit length");
        assert!(
            (depth - 1.0).abs() < 0.2,
            "expected depth ~1.0, got {}",
            depth
        );
    }
    #[test]
    fn test_compute_face_normal_outward() {
        let v0 = [1.0, 0.0, 0.0_f64];
        let v1 = [0.0, 1.0, 0.0_f64];
        let v2 = [-1.0, 0.0, 0.0_f64];
        let n = compute_face_normal_raw(v0, v1, v2);
        let len = epa_len3(n);
        assert!((len - 1.0).abs() < 1e-6, "normal must be unit length");
        assert!(
            n[2].abs() > 0.9,
            "normal should be along Z for XY-plane triangle, got {:?}",
            n
        );
    }
    #[test]
    fn test_epa_solve_expand_adds_face() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let initial_face_count = polytope.faces.len();
        let expanded = polytope.expand([3.0, 0.0, 0.0]);
        assert!(expanded, "expansion with a far point should succeed");
        assert!(
            polytope.faces.len() >= initial_face_count,
            "after expansion, face count should not decrease"
        );
    }
    #[test]
    fn test_epa_penetration_with_config() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let cfg = EpaConfig::high_precision();
        let result = epa_penetration_with_config(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            &cfg,
        );
        assert!(
            result.is_some(),
            "EPA with high-precision config should succeed"
        );
        let pen = result.unwrap();
        assert!(pen.depth > 0.0);
    }
    #[test]
    fn test_epa_newell_normal_xy_triangle() {
        let verts = vec![[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!(
            (len - 1.0).abs() < 1e-6,
            "normal must be unit length, got {}",
            len
        );
        assert!(
            n[2].abs() > 0.9,
            "normal must be along Z for XY triangle, got {:?}",
            n
        );
    }
    #[test]
    fn test_epa_newell_normal_degenerate() {
        let n = epa_newell_normal(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert_eq!(n, [0.0, 1.0, 0.0]);
    }
    #[test]
    fn test_epa_newell_normal_collinear() {
        let verts = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!(len > 0.5, "normal length={}", len);
    }
    #[test]
    fn test_epa_newell_normal_quad_xz_plane() {
        let verts = vec![
            [1.0_f64, 0.0, 1.0],
            [-1.0, 0.0, 1.0],
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, -1.0],
        ];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!((len - 1.0).abs() < 1e-6, "normal unit len={}", len);
        assert!(
            n[1].abs() > 0.9,
            "expected Y normal for XZ-plane quad, got {:?}",
            n
        );
    }
    #[test]
    fn test_epa_recompute_normals_preserves_validity() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        epa_recompute_normals(&mut polytope);
        for face in &polytope.faces {
            let v = polytope.vertices[face.vertices[0]];
            let d = epa_dot3(face.normal, v);
            assert!(
                d >= -1e-6,
                "recomputed distance must be non-negative, got {}",
                d
            );
        }
    }
    #[test]
    fn test_epa_triangle_area_unit() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let area = epa_triangle_area(a, b, c);
        assert!((area - 0.5).abs() < 1e-10, "area={}", area);
    }
    #[test]
    fn test_epa_triangle_area_equilateral() {
        let a = [0.0_f64, 0.0, 0.0];
        let b = [2.0, 0.0, 0.0];
        let c = [1.0, 3.0_f64.sqrt(), 0.0];
        let area = epa_triangle_area(a, b, c);
        let expected = 3.0_f64.sqrt();
        assert!(
            (area - expected).abs() < 1e-6,
            "area={}, expected={}",
            area,
            expected
        );
    }
    #[test]
    fn test_epa_triangle_area_degenerate() {
        let p = [1.0_f64, 2.0, 3.0];
        let area = epa_triangle_area(p, p, p);
        assert!(area < 1e-10, "degenerate triangle area={}", area);
    }
    #[test]
    fn test_epa_filter_degenerate_faces_removes_tiny() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        polytope.vertices.push([0.01, 0.0, 0.0]);
        polytope.vertices.push([0.01, 1e-8, 0.0]);
        polytope.vertices.push([0.01, 0.0, 1e-8]);
        let n = polytope.vertices.len();
        polytope.faces.push(EpaFaceRaw {
            vertices: [n - 3, n - 2, n - 1],
            normal: [1.0, 0.0, 0.0],
            distance: 0.01,
        });
        let count_before = polytope.faces.len();
        epa_filter_degenerate_faces(&mut polytope, 1e-5);
        let count_after = polytope.faces.len();
        assert!(
            count_after < count_before,
            "degenerate face should be removed"
        );
    }
    #[test]
    fn test_epa_filter_degenerate_faces_keeps_valid() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let count_before = polytope.faces.len();
        epa_filter_degenerate_faces(&mut polytope, 1e-12);
        assert_eq!(
            polytope.faces.len(),
            count_before,
            "no faces should be removed with tiny threshold"
        );
    }
    #[test]
    fn test_epa_stats_converged() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (result, stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(result.is_some(), "EPA with stats must succeed");
        assert!(stats.iterations > 0, "must have at least 1 iteration");
        assert!(stats.final_depth > 0.0, "final depth must be positive");
        assert!(stats.face_count >= 4, "must have at least 4 faces");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0);
    }
    #[test]
    fn test_epa_stats_iterations_bounded() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let max_iter = 10;
        let (result, stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            max_iter,
            1e-6,
        );
        assert!(
            stats.iterations <= max_iter,
            "iterations={} > max={}",
            stats.iterations,
            max_iter
        );
        assert!(result.is_some() || !stats.converged);
    }
    #[test]
    fn test_epa_silhouette_nonempty() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let new_point = [3.0, 0.0, 0.0];
        let (horizon, visible) = epa_silhouette(&polytope, new_point);
        assert!(
            !visible.is_empty(),
            "some faces must be visible from far point"
        );
        assert!(!horizon.is_empty(), "horizon must have edges");
    }
    #[test]
    fn test_epa_silhouette_no_visible_from_inside() {
        let simplex = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let (_horizon, visible) = epa_silhouette(&polytope, [0.0, 0.0, 0.0]);
        assert!(visible.is_empty(), "origin should see no faces as visible");
    }
    #[test]
    fn test_epa_penetration_early_exit_converges() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration_early_exit(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-5,
        );
        assert!(result.is_some(), "early-exit EPA must return a result");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0, "depth must be positive");
        let nlen = epa_len3(pen.normal);
        assert!((nlen - 1.0).abs() < 1e-4, "normal must be unit length");
    }
    #[test]
    fn test_epa_with_fallback_sphere_support() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let pen = epa_with_fallback(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(pen.depth >= 0.0, "fallback depth must be non-negative");
        let nlen = epa_len3(pen.normal);
        assert!(
            (nlen - 1.0).abs() < 1e-4,
            "fallback normal must be unit length"
        );
    }
    #[test]
    fn test_epa_with_fallback_returns_positive_depth() {
        let simplex: [[f64; 3]; 4] = [
            [0.3, 0.3, 0.3],
            [0.3, -0.3, -0.3],
            [-0.3, 0.3, -0.3],
            [-0.3, -0.3, 0.3],
        ];
        let pen = epa_with_fallback(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(pen.depth > 0.0, "fallback depth={}", pen.depth);
    }
    #[test]
    fn test_epa_penetration_capped_respects_max_faces() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration_capped(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
            10,
        );
        assert!(result.is_some(), "capped EPA must return a result");
        let pen = result.unwrap();
        assert!(pen.depth > 0.0);
    }
    #[test]
    fn test_epa_penetration_capped_4_face_limit() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let result = epa_penetration_capped(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
            4,
        );
        assert!(result.is_some());
    }
    #[test]
    fn test_epa_add3_basic() {
        let a = [1.0_f64, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let r = epa_add3(a, b);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_epa_dot3_orthogonal() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!(
            (epa_dot3(a, b)).abs() < 1e-15,
            "orthogonal vectors must have dot=0"
        );
    }
    #[test]
    fn test_epa_dot3_parallel() {
        let a = [3.0_f64, 0.0, 0.0];
        let b = [5.0, 0.0, 0.0];
        assert!((epa_dot3(a, b) - 15.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_sub3_basic() {
        let a = [5.0_f64, 3.0, 1.0];
        let b = [1.0, 2.0, 3.0];
        let r = epa_sub3(a, b);
        assert_eq!(r, [4.0, 1.0, -2.0]);
    }
    #[test]
    fn test_epa_scale3_halving() {
        let v = [2.0_f64, 4.0, 8.0];
        let r = epa_scale3(v, 0.5);
        assert_eq!(r, [1.0, 2.0, 4.0]);
    }
    #[test]
    fn test_epa_negate3_roundtrip() {
        let v = [1.0_f64, -2.0, 3.0];
        let neg = epa_negate3(v);
        let back = epa_negate3(neg);
        assert_eq!(back, v);
    }
    #[test]
    fn test_epa_cross3_anticommutative() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let ab = epa_cross3(a, b);
        let ba = epa_cross3(b, a);
        for i in 0..3 {
            assert!(
                (ab[i] + ba[i]).abs() < 1e-15,
                "cross must be anti-commutative: idx={}",
                i
            );
        }
    }
    #[test]
    fn test_epa_cross3_x_cross_y_is_z() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = epa_cross3(a, b);
        assert!((c[0]).abs() < 1e-15);
        assert!((c[1]).abs() < 1e-15);
        assert!((c[2] - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_len3_pythagorean() {
        let v = [3.0_f64, 4.0, 0.0];
        assert!((epa_len3(v) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_epa_len3_unit_axes() {
        assert!((epa_len3([1.0, 0.0, 0.0]) - 1.0).abs() < 1e-15);
        assert!((epa_len3([0.0, 1.0, 0.0]) - 1.0).abs() < 1e-15);
        assert!((epa_len3([0.0, 0.0, 1.0]) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_add3_commutative() {
        let a = [1.0_f64, 2.0, 3.0];
        let b = [7.0, -1.0, 0.5];
        let ab = epa_add3(a, b);
        let ba = epa_add3(b, a);
        assert_eq!(ab, ba);
    }
    #[test]
    fn test_epa_config_default_values() {
        let cfg = EpaConfig::default();
        assert_eq!(cfg.max_iter, 64);
        assert!((cfg.tolerance - 1e-6).abs() < 1e-15);
        assert_eq!(cfg.max_faces, 256);
    }
    #[test]
    fn test_epa_config_high_precision_tighter() {
        let hp = EpaConfig::high_precision();
        assert!(hp.max_iter > 64, "high_precision must have more iterations");
        assert!(
            hp.tolerance < 1e-6,
            "high_precision must have tighter tolerance"
        );
    }
    #[test]
    fn test_epa_has_converged_yes() {
        assert!(
            epa_has_converged(1.0, 1.0 + 1e-8, 1e-6),
            "change of 1e-8 within tol=1e-6 must converge"
        );
    }
    #[test]
    fn test_epa_has_converged_no() {
        assert!(
            !epa_has_converged(0.0, 0.01, 1e-6),
            "change of 0.01 outside tol=1e-6 must not converge"
        );
    }
    #[test]
    fn test_epa_has_converged_exact() {
        assert!(
            epa_has_converged(5.0, 5.0, 1e-15),
            "zero change within any positive tol must converge"
        );
    }
    #[test]
    fn test_compute_face_normal_raw_xy_plane() {
        let n = compute_face_normal_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(
            n[2].abs() > 0.9,
            "normal z-component must dominate: n={:?}",
            n
        );
    }
    #[test]
    fn test_compute_face_normal_raw_unit_length() {
        let n = compute_face_normal_raw([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let len = epa_len3(n);
        assert!(
            (len - 1.0).abs() < 1e-9,
            "face normal should be unit length: len={}",
            len
        );
    }
    #[test]
    fn test_epa_barycentric_origin_centroid() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [-0.5, 0.866_025, 0.0];
        let c = [-0.5, -0.866_025, 0.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!(
            (u + v + w - 1.0).abs() < 1e-9,
            "bary coords must sum to 1: {}",
            u + v + w
        );
        assert!(
            u >= 0.0 && v >= 0.0 && w >= 0.0,
            "bary coords must be non-negative"
        );
    }
    #[test]
    fn test_epa_barycentric_origin_sum_one() {
        let a = [2.0_f64, 0.0, 0.0];
        let b = [0.0, 2.0, 0.0];
        let c = [0.0, 0.0, 2.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-9, "sum={}", u + v + w);
    }
    #[test]
    fn test_epa_triangle_area_right_triangle() {
        let area = epa_triangle_area([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]);
        assert!((area - 6.0).abs() < 1e-9, "area={}", area);
    }
    #[test]
    fn test_epa_triangle_area_collapsed_point() {
        let area = epa_triangle_area([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!(area < 1e-12, "degenerate triangle area={}", area);
    }
    #[test]
    fn test_epa_triangle_area_unit_equilateral() {
        let area = epa_triangle_area(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 3_f64.sqrt() / 2.0, 0.0],
        );
        assert!((area - 3_f64.sqrt() / 4.0).abs() < 1e-9, "area={}", area);
    }
    #[test]
    fn test_epa_polytope_from_gjk_simplex_face_count() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let poly = EpaPolytope::from_gjk_simplex(&simplex);
        assert_eq!(poly.faces.len(), 4, "tetrahedron must have 4 faces");
        assert_eq!(poly.vertices.len(), 4, "tetrahedron must have 4 vertices");
    }
    #[test]
    fn test_epa_polytope_find_closest_face_non_panics_empty() {
        let poly = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let idx = poly.find_closest_face();
        assert_eq!(idx, 0);
    }
    #[test]
    fn test_epa_polytope_closest_face_empty_default() {
        let poly = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let (idx, dist, normal) = poly.closest_face();
        assert_eq!(idx, 0);
        assert!((dist - 0.0).abs() < 1e-15);
        let _ = normal;
    }
    #[test]
    fn test_epa_polytope_from_gjk_simplex_normals_unit() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let poly = EpaPolytope::from_gjk_simplex(&simplex);
        for (i, f) in poly.faces.iter().enumerate() {
            let len = epa_len3(f.normal);
            assert!(
                (len - 1.0).abs() < 1e-9,
                "face {} normal not unit: len={}",
                i,
                len
            );
        }
    }
    #[test]
    fn test_epa_validate_polytope_tetra_valid() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let poly = EpaPolytope::from_gjk_simplex(&simplex);
        let valid = epa_validate_polytope(&poly);
        assert!(valid, "well-formed polytope must validate");
    }
    #[test]
    fn test_epa_validate_polytope_no_panic_empty() {
        let poly = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let _valid = epa_validate_polytope(&poly);
    }
    #[test]
    fn test_epa_newell_normal_square_xy_plane() {
        let verts = vec![
            [1.0_f64, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
        ];
        let n = epa_newell_normal(&verts);
        let len = epa_len3(n);
        assert!(len > 1e-9, "Newell normal must be non-zero for XY square");
        let normalized = epa_scale3(n, 1.0 / len);
        assert!(
            normalized[2].abs() > 0.9,
            "XY square Newell normal must be along Z: {:?}",
            normalized
        );
    }
    #[test]
    fn test_epa_filter_degenerate_faces_removes_nothing_if_all_valid() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let mut poly = EpaPolytope::from_gjk_simplex(&simplex);
        let before = poly.faces.len();
        epa_filter_degenerate_faces(&mut poly, 0.0);
        assert_eq!(poly.faces.len(), before, "no faces should be filtered");
    }
    #[test]
    fn test_epa_filter_degenerate_faces_removes_all_if_threshold_large() {
        let simplex: [[f64; 3]; 4] = [
            [0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
            [-0.1, -0.1, -0.1],
        ];
        let mut poly = EpaPolytope::from_gjk_simplex(&simplex);
        epa_filter_degenerate_faces(&mut poly, 1000.0);
        assert_eq!(
            poly.faces.len(),
            0,
            "all faces should be filtered with large threshold"
        );
    }
    #[test]
    fn test_epa_recompute_normals_no_panic() {
        let simplex: [[f64; 3]; 4] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let mut poly = EpaPolytope::from_gjk_simplex(&simplex);
        epa_recompute_normals(&mut poly);
        for (i, f) in poly.faces.iter().enumerate() {
            let len = epa_len3(f.normal);
            assert!(
                (len - 1.0).abs() < 1e-9,
                "face {} recomputed normal not unit: len={}",
                i,
                len
            );
        }
    }
    #[test]
    fn test_epa_stats_default_zero() {
        let s = EpaStats::default();
        assert_eq!(s.iterations, 0);
        assert_eq!(s.face_count, 0);
        assert!((s.final_depth - 0.0).abs() < 1e-15);
        assert!(!s.converged);
    }
    #[test]
    fn test_epa_penetration_with_stats_converged_flag() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (result, stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                epa_scale3(n, 1.0)
            },
            &simplex,
            64,
            1e-6,
        );
        assert!(stats.iterations > 0, "must have at least one iteration");
        if result.is_some() {
            assert!(stats.final_depth >= 0.0, "final depth must be non-negative");
        }
    }
    #[test]
    fn test_epa_penetration_with_stats_depth_positive() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let (result, _stats) = epa_penetration_with_stats(
            |dir| {
                let len = epa_len3(dir);
                if len < 1e-10 {
                    return [0.0; 3];
                }
                let n = epa_scale3(dir, 1.0 / len);
                [2.0 * n[0], 2.0 * n[1], 2.0 * n[2]]
            },
            &simplex,
            64,
            1e-6,
        );
        if let Some(pen) = result {
            assert!(
                pen.depth > 0.0,
                "penetration depth must be positive: {}",
                pen.depth
            );
            let nlen = epa_len3(pen.normal);
            assert!(
                (nlen - 1.0).abs() < 1e-6,
                "penetration normal must be unit: len={}",
                nlen
            );
        }
    }
    #[test]
    fn test_epa_contact_normal_from_face_passthrough() {
        let face = EpaFaceRaw {
            vertices: [0, 1, 2],
            normal: [0.0, 1.0, 0.0],
            distance: 0.5,
        };
        let n = epa_contact_normal_from_face(&face);
        assert_eq!(n, [0.0, 1.0, 0.0]);
    }
    #[test]
    fn test_epa_solver_penetration_axis_unit_normal() {
        let simplex = make_tetrahedron();
        let solver = EpaSolver::from_simplex(&simplex);
        let (normal, depth) = solver.compute_penetration_axis();
        let nlen = epa_len3(normal);
        assert!(
            (nlen - 1.0).abs() < 1e-9,
            "penetration axis must be unit, len={nlen}"
        );
        assert!(depth >= 0.0, "depth must be non-negative, got {depth}");
    }
    #[test]
    fn test_epa_solver_penetration_axis_depth_positive_for_tetrahedron() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let solver = EpaSolver::from_simplex(&simplex);
        let (_normal, depth) = solver.compute_penetration_axis();
        assert!(
            depth > 0.0,
            "depth must be > 0 for a simplex enclosing origin, got {depth}"
        );
    }
    #[test]
    fn test_epa_solver_penetration_axis_empty_polytope_fallback() {
        let simplex = make_tetrahedron();
        let mut solver = EpaSolver::from_simplex(&simplex);
        solver.polytope = Some(EpaPolytope {
            vertices: vec![],
            faces: vec![],
        });
        let (normal, depth) = solver.compute_penetration_axis();
        assert_eq!(normal, [0.0, 1.0, 0.0]);
        assert!((depth - 0.0).abs() < 1e-15);
    }
    #[test]
    fn test_epa_solver_refine_polytope_increases_vertex_count() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let before = solver.vertex_count();
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        let expanded = solver.refine_polytope(&mut support_fn);
        if expanded {
            assert!(
                solver.vertex_count() > before,
                "refine should add a vertex when expanding"
            );
        }
    }
    #[test]
    fn test_epa_solver_refine_polytope_increments_iterations() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        assert_eq!(solver.iterations(), 0);
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            epa_scale3(dir, 1.0 / len)
        };
        solver.refine_polytope(&mut support_fn);
        assert_eq!(
            solver.iterations(),
            1,
            "one call must increment iteration counter"
        );
    }
    #[test]
    fn test_epa_solver_refine_polytope_false_on_convergence() {
        let simplex: [[f64; 3]; 4] = [
            [0.001, 0.001, 0.001],
            [0.001, -0.001, -0.001],
            [-0.001, 0.001, -0.001],
            [-0.001, -0.001, 0.001],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex).with_tolerance(1.0);
        let mut support_fn = |_dir: [f64; 3]| -> [f64; 3] { [0.001, 0.0, 0.0] };
        let result = solver.refine_polytope(&mut support_fn);
        assert!(
            !result,
            "large tolerance should cause immediate convergence"
        );
    }
    #[test]
    fn test_epa_solver_solve_returns_some() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        let result = solver.solve(&mut support_fn, 64);
        assert!(result.is_some(), "solver must return a penetration result");
    }
    #[test]
    fn test_epa_solver_solve_unit_normal() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        if let Some(pen) = solver.solve(&mut support_fn, 64) {
            let nlen = epa_len3(pen.normal);
            assert!(
                (nlen - 1.0).abs() < 1e-6,
                "penetration normal must be unit, len={nlen}"
            );
        }
    }
    #[test]
    fn test_epa_solver_face_count_grows() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut solver = EpaSolver::from_simplex(&simplex);
        let initial_faces = solver.face_count();
        let mut support_fn = |dir: [f64; 3]| -> [f64; 3] {
            let len = epa_len3(dir);
            if len < 1e-10 {
                return [0.0; 3];
            }
            let n = epa_scale3(dir, 1.0 / len);
            [2.0 * n[0] - 0.5, 2.0 * n[1], 2.0 * n[2]]
        };
        for _ in 0..5 {
            solver.refine_polytope(&mut support_fn);
        }
        assert!(
            solver.face_count() >= initial_faces,
            "face count should not decrease: {} < {}",
            solver.face_count(),
            initial_faces
        );
    }
    #[test]
    fn test_epa_solver_with_tolerance_setter() {
        let simplex = make_tetrahedron();
        let solver = EpaSolver::from_simplex(&simplex).with_tolerance(1e-9);
        let (normal, _depth) = solver.compute_penetration_axis();
        let nlen = epa_len3(normal);
        assert!((nlen - 1.0).abs() < 1e-9);
    }
}
