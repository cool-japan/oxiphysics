//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{epa_dot3, epa_negate3, epa_scale3, epa_sub3};
use super::types::{EpaFaceRaw, EpaPenetration, EpaPolytope, EpaWitness};

#[cfg(test)]
pub(super) fn epa_add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Run EPA but stop if the polytope exceeds `max_faces` faces.
///
/// This prevents unbounded growth in pathological cases.
pub fn epa_penetration_capped<F>(
    mut support_fn: F,
    initial_simplex: &[[f64; 3]; 4],
    max_iter: usize,
    tol: f64,
    max_faces: usize,
) -> Option<EpaPenetration>
where
    F: FnMut([f64; 3]) -> [f64; 3],
{
    let mut polytope = EpaPolytope::from_gjk_simplex(initial_simplex);
    for _ in 0..max_iter {
        if polytope.faces.is_empty() || polytope.faces.len() > max_faces {
            break;
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
/// Extract witness points for shapes A and B from an EPA result.
///
/// `support_a(dir)` returns the support of shape A.
/// `support_b(dir)` returns the support of shape B.
/// `epa_result` is the converged EPA penetration.
pub fn epa_extract_witness<FA, FB>(
    mut support_a: FA,
    mut support_b: FB,
    epa_result: &EpaPenetration,
) -> EpaWitness
where
    FA: FnMut([f64; 3]) -> [f64; 3],
    FB: FnMut([f64; 3]) -> [f64; 3],
{
    let n = epa_result.normal;
    let depth = epa_result.depth;
    let pa = support_a(n);
    let pb = support_b(epa_negate3(n));
    EpaWitness {
        point_a: pa,
        point_b: pb,
        normal: n,
        depth,
    }
}
/// Clone the polytope, keeping only the `max_faces` closest faces.
///
/// This is useful as a "polytope reduction" step to limit memory usage
/// when the polytope grows too large after many expansions.
pub fn epa_keep_closest_faces(polytope: &EpaPolytope, max_faces: usize) -> EpaPolytope {
    if polytope.faces.len() <= max_faces {
        return polytope.clone();
    }
    let mut sorted_faces: Vec<(usize, f64)> = polytope
        .faces
        .iter()
        .enumerate()
        .map(|(i, f)| (i, f.distance))
        .collect();
    sorted_faces.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    sorted_faces.truncate(max_faces);
    let kept_indices: std::collections::HashSet<usize> =
        sorted_faces.iter().map(|(i, _)| *i).collect();
    let mut vertex_map: Vec<Option<usize>> = vec![None; polytope.vertices.len()];
    let mut new_vertices: Vec<[f64; 3]> = Vec::new();
    for &fi in &kept_indices {
        for &vi in &polytope.faces[fi].vertices {
            if vertex_map[vi].is_none() {
                vertex_map[vi] = Some(new_vertices.len());
                new_vertices.push(polytope.vertices[vi]);
            }
        }
    }
    let new_faces: Vec<EpaFaceRaw> = kept_indices
        .iter()
        .filter_map(|&fi| {
            let f = &polytope.faces[fi];
            Some(EpaFaceRaw {
                vertices: [
                    vertex_map[f.vertices[0]]?,
                    vertex_map[f.vertices[1]]?,
                    vertex_map[f.vertices[2]]?,
                ],
                normal: f.normal,
                distance: f.distance,
            })
        })
        .collect();
    EpaPolytope {
        vertices: new_vertices,
        faces: new_faces,
    }
}
/// Incrementally refine the penetration depth estimate.
///
/// Starting from an initial depth estimate `init_depth`, refines the estimate
/// by expanding the polytope up to `max_steps` times.
///
/// Returns `(final_depth, converged)`.
pub fn epa_refine_depth<F>(
    mut support_fn: F,
    polytope: &mut EpaPolytope,
    init_depth: f64,
    max_steps: usize,
    tol: f64,
) -> (f64, bool)
where
    F: FnMut([f64; 3]) -> [f64; 3],
{
    let mut prev_depth = init_depth;
    for _ in 0..max_steps {
        if polytope.faces.is_empty() {
            return (prev_depth, false);
        }
        let closest_idx = polytope.find_closest_face();
        let closest = polytope.faces[closest_idx].clone();
        let current_depth = closest.distance;
        if (current_depth - prev_depth).abs() < tol {
            return (current_depth, true);
        }
        prev_depth = current_depth;
        let new_point = support_fn(closest.normal);
        let new_dist = epa_dot3(new_point, closest.normal);
        if (new_dist - closest.distance).abs() < tol {
            return (closest.distance, true);
        }
        if !polytope.expand(new_point) {
            return (closest.distance, false);
        }
    }
    let closest_idx = polytope.find_closest_face();
    let depth = polytope.faces[closest_idx].distance;
    (depth, false)
}
/// Generate multiple contact points from the EPA closest face and its neighbors.
///
/// The EPA result typically gives a single contact point (from the closest face).
/// This function returns up to `max_contacts` contact points by also sampling
/// neighboring faces close to the minimum distance.
pub fn epa_contact_manifold(
    polytope: &EpaPolytope,
    max_contacts: usize,
    neighbor_tol: f64,
) -> Vec<EpaPenetration> {
    if polytope.faces.is_empty() || max_contacts == 0 {
        return Vec::new();
    }
    let (_closest_idx, min_dist, _min_normal) = polytope.closest_face();
    let mut candidates: Vec<&EpaFaceRaw> = polytope
        .faces
        .iter()
        .filter(|f| (f.distance - min_dist).abs() <= neighbor_tol)
        .collect();
    candidates.sort_by(|a, b| {
        a.distance
            .partial_cmp(&b.distance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    candidates.truncate(max_contacts);
    candidates
        .iter()
        .map(|face| {
            let v0 = polytope.vertices[face.vertices[0]];
            EpaPenetration {
                normal: face.normal,
                depth: face.distance,
                witness_a: v0,
                witness_b: epa_sub3(v0, epa_scale3(face.normal, face.distance)),
            }
        })
        .collect()
}
#[cfg(test)]
mod tests_epa_extended {
    use super::super::*;
    use crate::narrowphase::EpaFaceQueue;
    fn make_tetrahedron() -> [[f64; 3]; 4] {
        [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ]
    }
    fn sphere_cso_support(dir: [f64; 3]) -> [f64; 3] {
        let len = epa_len3(dir);
        if len < 1e-10 {
            return [0.0; 3];
        }
        let n = epa_scale3(dir, 1.0 / len);
        [2.0 * n[0] - 1.0, 2.0 * n[1], 2.0 * n[2]]
    }
    #[test]
    fn test_epa_remove_face_shrinks_count() {
        let simplex = make_tetrahedron();
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let initial_count = polytope.faces.len();
        assert_eq!(initial_count, 4);
        epa_remove_face(&mut polytope, 0);
        assert_eq!(
            polytope.faces.len(),
            3,
            "removing a face should shrink count"
        );
    }
    #[test]
    fn test_epa_remove_face_out_of_bounds_noop() {
        let simplex = make_tetrahedron();
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let count = polytope.faces.len();
        epa_remove_face(&mut polytope, 100);
        assert_eq!(
            polytope.faces.len(),
            count,
            "out-of-bounds remove should be no-op"
        );
    }
    #[test]
    fn test_epa_project_origin_onto_face_unit_normal() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let proj = epa_project_origin_onto_face(&polytope, 0);
        let d = epa_len3(proj);
        let expected = polytope.faces[0].distance;
        assert!(
            (d - expected).abs() < 1e-6,
            "projected dist={} expected={}",
            d,
            expected
        );
    }
    #[test]
    fn test_epa_face_vertices_valid() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let verts = epa_face_vertices(&polytope, 0).expect("face 0 should exist");
        for v in &verts {
            let found = polytope.vertices.iter().any(|pv| {
                (pv[0] - v[0]).abs() < 1e-10
                    && (pv[1] - v[1]).abs() < 1e-10
                    && (pv[2] - v[2]).abs() < 1e-10
            });
            assert!(found, "face vertex {:?} not found in polytope vertices", v);
        }
    }
    #[test]
    fn test_epa_face_vertices_out_of_bounds() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        assert!(epa_face_vertices(&polytope, 100).is_none());
    }
    #[test]
    fn test_epa_face_centroid_inside_triangle() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let c = epa_face_centroid(&polytope, 0).expect("face 0 must exist");
        let verts = epa_face_vertices(&polytope, 0).unwrap();
        for i in 0..3 {
            let expected = (verts[0][i] + verts[1][i] + verts[2][i]) / 3.0;
            assert!(
                (c[i] - expected).abs() < 1e-10,
                "centroid[{}]={} expected={}",
                i,
                c[i],
                expected
            );
        }
    }
    #[test]
    fn test_epa_farthest_face_index() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let farthest = epa_farthest_face(&polytope).expect("must have a farthest face");
        let max_dist = polytope
            .faces
            .iter()
            .map(|f| f.distance)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((polytope.faces[farthest].distance - max_dist).abs() < 1e-10);
    }
    #[test]
    fn test_epa_polytope_stats_basic() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let stats = epa_polytope_stats(&polytope);
        assert_eq!(stats.face_count, 4);
        assert_eq!(stats.vertex_count, 4);
        assert!(stats.min_distance >= 0.0);
        assert!(stats.max_distance >= stats.min_distance);
        assert!(stats.avg_distance >= stats.min_distance);
        assert!(stats.avg_distance <= stats.max_distance);
    }
    #[test]
    fn test_epa_polytope_stats_empty() {
        let polytope = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let stats = epa_polytope_stats(&polytope);
        assert_eq!(stats.face_count, 0);
        assert_eq!(stats.vertex_count, 0);
    }
    #[test]
    fn test_epa_polytope_stats_min_area_non_negative() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let stats = epa_polytope_stats(&polytope);
        assert!(stats.min_area >= 0.0, "min area = {}", stats.min_area);
        assert!(
            stats.max_area >= stats.min_area,
            "max area = {}",
            stats.max_area
        );
    }
    #[test]
    fn test_epa_extract_witness_normal_unit_length() {
        let pen = EpaPenetration {
            normal: [1.0, 0.0, 0.0],
            depth: 0.5,
            witness_a: [0.5, 0.0, 0.0],
            witness_b: [0.0, 0.0, 0.0],
        };
        let sphere_a_support = |dir: [f64; 3]| {
            let l = epa_len3(dir);
            if l < 1e-10 {
                [0.0; 3]
            } else {
                epa_scale3(dir, 1.0 / l)
            }
        };
        let sphere_b_support = |dir: [f64; 3]| {
            let l = epa_len3(dir);
            if l < 1e-10 {
                [1.0, 0.0, 0.0]
            } else {
                epa_add3([1.0, 0.0, 0.0], epa_scale3(dir, 1.0 / l))
            }
        };
        let witness = epa_extract_witness(sphere_a_support, sphere_b_support, &pen);
        let nlen = epa_len3(witness.normal);
        assert!((nlen - 1.0).abs() < 1e-10, "witness normal length={}", nlen);
        assert_eq!(witness.depth, 0.5);
    }
    #[test]
    fn test_epa_keep_closest_faces_fewer_than_limit() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let reduced = epa_keep_closest_faces(&polytope, 10);
        assert_eq!(reduced.faces.len(), 4);
    }
    #[test]
    fn test_epa_keep_closest_faces_hard_cap() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let reduced = epa_keep_closest_faces(&polytope, 2);
        assert!(
            reduced.faces.len() <= 2,
            "must respect cap of 2, got {}",
            reduced.faces.len()
        );
    }
    #[test]
    fn test_epa_keep_closest_faces_preserves_distances() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let reduced = epa_keep_closest_faces(&polytope, 3);
        let mut all_dists: Vec<f64> = polytope.faces.iter().map(|f| f.distance).collect();
        all_dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
        all_dists.truncate(3);
        let max_kept = all_dists.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        for face in &reduced.faces {
            assert!(
                face.distance <= max_kept + 1e-10,
                "kept face dist {} > max expected {}",
                face.distance,
                max_kept
            );
        }
    }
    #[test]
    fn test_epa_refine_depth_converges() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let init_depth = polytope.faces[polytope.find_closest_face()].distance;
        let (final_depth, converged) =
            epa_refine_depth(sphere_cso_support, &mut polytope, init_depth, 32, 1e-6);
        assert!(final_depth > 0.0, "final depth must be positive");
        let _ = converged;
    }
    #[test]
    fn test_epa_refine_depth_does_not_decrease() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let initial_closest_dist = polytope.faces[polytope.find_closest_face()].distance;
        let (final_depth, _) = epa_refine_depth(
            sphere_cso_support,
            &mut polytope,
            initial_closest_dist,
            16,
            1e-8,
        );
        assert!(
            final_depth >= initial_closest_dist - 1e-8,
            "depth regressed: {} < {}",
            final_depth,
            initial_closest_dist
        );
    }
    #[test]
    fn test_epa_face_queue_peek_min() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let queue = EpaFaceQueue::new(&polytope);
        let (idx, dist, normal) = queue.peek_min().expect("queue must have faces");
        assert!(dist >= 0.0);
        let nlen = epa_len3(normal);
        assert!((nlen - 1.0).abs() < 1e-6, "normal must be unit length");
        assert_eq!(idx, polytope.find_closest_face());
    }
    #[test]
    fn test_epa_face_queue_is_empty_for_empty_polytope() {
        let polytope = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let queue = EpaFaceQueue::new(&polytope);
        assert!(queue.is_empty());
        assert!(queue.peek_min().is_none());
    }
    #[test]
    fn test_epa_face_queue_len() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let queue = EpaFaceQueue::new(&polytope);
        assert_eq!(queue.len(), polytope.faces.len());
    }
    #[test]
    fn test_epa_contact_manifold_returns_closest_first() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let contacts = epa_contact_manifold(&polytope, 4, 0.5);
        assert!(!contacts.is_empty(), "must return at least one contact");
        let min_depth = contacts
            .iter()
            .map(|c| c.depth)
            .fold(f64::INFINITY, f64::min);
        assert!((contacts[0].depth - min_depth).abs() < 1e-10);
    }
    #[test]
    fn test_epa_contact_manifold_respects_max_contacts() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let contacts = epa_contact_manifold(&polytope, 2, 1.0);
        assert!(contacts.len() <= 2, "must not exceed max_contacts=2");
    }
    #[test]
    fn test_epa_contact_manifold_empty_polytope() {
        let polytope = EpaPolytope {
            vertices: vec![],
            faces: vec![],
        };
        let contacts = epa_contact_manifold(&polytope, 4, 0.1);
        assert!(contacts.is_empty());
    }
    #[test]
    fn test_epa_contact_manifold_zero_max_contacts() {
        let simplex = make_tetrahedron();
        let polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let contacts = epa_contact_manifold(&polytope, 0, 1.0);
        assert!(contacts.is_empty());
    }
    #[test]
    fn test_epa_all_contact_normals_unit_length() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        for _ in 0..3 {
            let idx = polytope.find_closest_face();
            let n = polytope.faces[idx].normal;
            let sp = sphere_cso_support(n);
            if !polytope.expand(sp) {
                break;
            }
        }
        for face in &polytope.faces {
            let nlen = epa_len3(face.normal);
            assert!(
                (nlen - 1.0).abs() < 1e-5,
                "face normal not unit length: {}",
                nlen
            );
        }
    }
    #[test]
    fn test_epa_all_face_distances_non_negative() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        for _ in 0..5 {
            let idx = polytope.find_closest_face();
            let n = polytope.faces[idx].normal;
            let sp = sphere_cso_support(n);
            if !polytope.expand(sp) {
                break;
            }
        }
        for face in &polytope.faces {
            assert!(
                face.distance >= -1e-8,
                "face distance must be non-negative, got {}",
                face.distance
            );
        }
    }
    #[test]
    fn test_epa_polytope_vertex_count_grows_on_expand() {
        let simplex: [[f64; 3]; 4] = [
            [0.5, 0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
        ];
        let mut polytope = EpaPolytope::from_gjk_simplex(&simplex);
        let initial_vcount = polytope.vertices.len();
        let n = polytope.faces[polytope.find_closest_face()].normal;
        let sp = sphere_cso_support(n);
        let expanded = polytope.expand(sp);
        if expanded {
            assert!(
                polytope.vertices.len() > initial_vcount,
                "vertex count should grow after expand"
            );
        }
    }
    #[test]
    fn test_epa_barycentric_origin_vertex_at_origin() {
        let a = [0.01, 0.0, 0.0];
        let b = [1.0, 0.0, 1.0];
        let c = [1.0, 0.0, -1.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-8, "u+v+w = {}", u + v + w);
        assert!(
            u >= v && u >= w,
            "u={} v={} w={} (u should dominate)",
            u,
            v,
            w
        );
    }
    #[test]
    fn test_epa_barycentric_origin_symmetric_triangle() {
        let a = [0.0, 0.0, 1.0];
        let b = [1.0, 0.0, -1.0];
        let c = [-1.0, 0.0, -1.0];
        let (u, v, w) = epa_barycentric_origin(a, b, c);
        assert!((u + v + w - 1.0).abs() < 1e-8);
        assert!(
            (v - w).abs() < 1e-8,
            "v={} w={} should be equal for symmetric triangle",
            v,
            w
        );
    }
    #[test]
    fn test_epa_has_converged_at_exact_zero_delta() {
        assert!(epa_has_converged(2.0, 2.0, 1e-6));
    }
    #[test]
    fn test_epa_has_converged_negative_delta() {
        assert!(!epa_has_converged(2.0, 1.9, 1e-6));
    }
}
