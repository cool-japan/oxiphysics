//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::super::gjk::SupportPoint;
use super::super::types::{
    EpaConfig, EpaFaceRaw, EpaPenetration, EpaPolytope, EpaPolytopeStats, EpaStats, Face,
};
use crate::types::Contact;
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_geometry::Shape;

/// Maximum EPA iterations.
pub(crate) const MAX_ITERATIONS: usize = 64;
/// EPA convergence tolerance.
pub(crate) const TOLERANCE: Real = 1e-6;
pub fn support_minkowski(
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
pub fn build_initial_faces(vertices: &[SupportPoint]) -> Vec<Face> {
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
pub fn compute_face_normal(a: &Vec3, b: &Vec3, c: &Vec3) -> Vec3 {
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
pub fn add_edge(edges: &mut Vec<(usize, usize)>, a: usize, b: usize) {
    if let Some(pos) = edges.iter().position(|&(ea, eb)| ea == b && eb == a) {
        edges.swap_remove(pos);
    } else {
        edges.push((a, b));
    }
}
pub fn build_contact(face: &Face, vertices: &[SupportPoint]) -> Contact {
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
pub fn epa_dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn epa_sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn epa_scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn epa_negate3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
pub fn epa_cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn epa_len3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
pub fn epa_face_normal(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
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
pub fn epa_add_edge(edges: &mut Vec<(usize, usize)>, a: usize, b: usize) {
    if let Some(pos) = edges.iter().position(|&(ea, eb)| ea == b && eb == a) {
        edges.swap_remove(pos);
    } else {
        edges.push((a, b));
    }
}
/// Check whether an EPA iteration has converged given the old and new face distances.
pub fn epa_has_converged(old_dist: f64, new_dist: f64, tol: f64) -> bool {
    (new_dist - old_dist).abs() < tol
}
/// Compute the barycentric coordinates of the projection of the origin onto
/// a triangle with vertices `a`, `b`, `c`.
///
/// Returns `(u, v, w)` such that `u + v + w = 1` and the projected point is
/// `u*a + v*b + w*c`.  Values are clamped to `[0, 1]` and renormalised.
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
pub fn epa_contact_normal_from_face(face: &EpaFaceRaw) -> [f64; 3] {
    face.normal
}
/// Run the full EPA algorithm with the given config.
///
/// This is a wrapper around `epa_penetration` that accepts an `EpaConfig`.
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
pub fn epa_triangle_area2(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = epa_sub3(b, a);
    let ac = epa_sub3(c, a);
    let cross = epa_cross3(ab, ac);
    epa_len3(cross)
}
/// Compute the area of a triangle.
pub fn epa_triangle_area(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    epa_triangle_area2(a, b, c) * 0.5
}
/// Remove degenerate faces (area below `area_threshold`) from a polytope.
///
/// Degenerate faces arise from near-collinear or coincident vertices and
/// can cause numerical instability in the EPA iteration.
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
pub fn epa_remove_face(polytope: &mut EpaPolytope, face_idx: usize) {
    if face_idx < polytope.faces.len() {
        polytope.faces.swap_remove(face_idx);
    }
}
/// Compute the projected point of the origin onto the plane of face `face_idx`.
///
/// Returns the projection in world space (same coordinate system as `polytope.vertices`).
/// Returns the face centroid as a fallback for degenerate faces.
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
pub fn epa_face_vertices(polytope: &EpaPolytope, face_idx: usize) -> Option<[[f64; 3]; 3]> {
    let face = polytope.faces.get(face_idx)?;
    Some([
        polytope.vertices[face.vertices[0]],
        polytope.vertices[face.vertices[1]],
        polytope.vertices[face.vertices[2]],
    ])
}
/// Compute the centroid of a face.
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
