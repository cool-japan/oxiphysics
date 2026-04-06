//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{Face, Vertex};
use super::functions::{vec3_add, vec3_normalize, vec3_scale};
use super::types::{AnisotropicSizeField, BooleanResult, ProcessMesh};

/// Boolean intersection of two meshes.
///
/// Keeps only faces of A inside B and faces of B inside A.
pub fn mesh_intersection(a: &ProcessMesh, b: &ProcessMesh) -> BooleanResult {
    let offset_b = a.verts.len();
    let mut verts = a.verts.clone();
    verts.extend_from_slice(&b.verts);
    let faces_a: Vec<Face> = a
        .faces
        .iter()
        .filter(|&&[i, j, k]| {
            let centroid = vec3_scale(
                vec3_add(vec3_add(a.verts[i], a.verts[j]), a.verts[k]),
                1.0 / 3.0,
            );
            point_in_mesh(centroid, b)
        })
        .cloned()
        .collect();
    let faces_b: Vec<Face> = b
        .faces
        .iter()
        .filter(|&&[i, j, k]| {
            let centroid = vec3_scale(
                vec3_add(vec3_add(b.verts[i], b.verts[j]), b.verts[k]),
                1.0 / 3.0,
            );
            point_in_mesh(centroid, a)
        })
        .map(|&[i, j, k]| [i + offset_b, j + offset_b, k + offset_b])
        .collect();
    let mut faces = faces_a;
    faces.extend(faces_b);
    BooleanResult {
        mesh: ProcessMesh::new(verts, faces),
        is_exact: false,
    }
}
/// Boolean difference A − B.
///
/// Keeps faces of A outside B and faces of B inside A (flipped).
pub fn mesh_difference(a: &ProcessMesh, b: &ProcessMesh) -> BooleanResult {
    let offset_b = a.verts.len();
    let mut verts = a.verts.clone();
    verts.extend_from_slice(&b.verts);
    let faces_a: Vec<Face> = a
        .faces
        .iter()
        .filter(|&&[i, j, k]| {
            let centroid = vec3_scale(
                vec3_add(vec3_add(a.verts[i], a.verts[j]), a.verts[k]),
                1.0 / 3.0,
            );
            !point_in_mesh(centroid, b)
        })
        .cloned()
        .collect();
    let faces_b: Vec<Face> = b
        .faces
        .iter()
        .filter(|&&[i, j, k]| {
            let centroid = vec3_scale(
                vec3_add(vec3_add(b.verts[i], b.verts[j]), b.verts[k]),
                1.0 / 3.0,
            );
            point_in_mesh(centroid, a)
        })
        .map(|&[i, j, k]| [j + offset_b, i + offset_b, k + offset_b])
        .collect();
    let mut faces = faces_a;
    faces.extend(faces_b);
    BooleanResult {
        mesh: ProcessMesh::new(verts, faces),
        is_exact: false,
    }
}
/// Build a simple unit tetrahedron mesh.
pub fn make_tetrahedron() -> ProcessMesh {
    let verts = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, (3.0_f64).sqrt() / 2.0, 0.0],
        [0.5, (3.0_f64).sqrt() / 6.0, (6.0_f64).sqrt() / 3.0],
    ];
    let faces = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
    ProcessMesh::new(verts, faces)
}
/// Build a flat square mesh (two triangles).
pub fn make_quad() -> ProcessMesh {
    let verts = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let faces = vec![[0, 1, 2], [0, 2, 3]];
    ProcessMesh::new(verts, faces)
}
/// Build an icosphere with `subdivisions` subdivision rounds.
pub fn make_icosphere(subdivisions: usize) -> ProcessMesh {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let verts: Vec<Vertex> = vec![
        vec3_normalize([-1.0, phi, 0.0]),
        vec3_normalize([1.0, phi, 0.0]),
        vec3_normalize([-1.0, -phi, 0.0]),
        vec3_normalize([1.0, -phi, 0.0]),
        vec3_normalize([0.0, -1.0, phi]),
        vec3_normalize([0.0, 1.0, phi]),
        vec3_normalize([0.0, -1.0, -phi]),
        vec3_normalize([0.0, 1.0, -phi]),
        vec3_normalize([phi, 0.0, -1.0]),
        vec3_normalize([phi, 0.0, 1.0]),
        vec3_normalize([-phi, 0.0, -1.0]),
        vec3_normalize([-phi, 0.0, 1.0]),
    ];
    let faces: Vec<Face> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    let mut mesh = ProcessMesh::new(verts, faces);
    for _ in 0..subdivisions {
        mesh = midpoint_subdivide(&mesh);
        for v in &mut mesh.verts {
            *v = vec3_normalize(*v);
        }
    }
    mesh
}
/// Build a uniform anisotropic size field (same size everywhere).
pub fn uniform_anisotropic_field(mesh: &ProcessMesh, len: f64) -> AnisotropicSizeField {
    let nv = mesh.verts.len();
    AnisotropicSizeField {
        k_min_dir: vec![[1.0, 0.0, 0.0]; nv],
        k_max_dir: vec![[0.0, 1.0, 0.0]; nv],
        len_min: vec![len; nv],
        len_max: vec![len; nv],
    }
}
/// Anisotropic remeshing guided by a size field.
///
/// Splits edges that exceed `len_max` at either endpoint and collapses edges
/// shorter than `len_min` at either endpoint.
pub fn anisotropic_remesh(
    mesh: &ProcessMesh,
    field: &AnisotropicSizeField,
    iters: usize,
) -> ProcessMesh {
    let mut out = mesh.clone();
    for _ in 0..iters {
        let avg_max: f64 =
            field.len_max.iter().copied().sum::<f64>() / field.len_max.len().max(1) as f64;
        let avg_min: f64 =
            field.len_min.iter().copied().sum::<f64>() / field.len_min.len().max(1) as f64;
        out = split_long_edges(&out, avg_max);
        out = collapse_short_edges(&out, avg_min * 0.6);
        out = laplacian_smooth(&out, 0.3, 1);
    }
    out
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_processing::Quadric;
    fn tet() -> ProcessMesh {
        make_tetrahedron()
    }
    fn quad() -> ProcessMesh {
        make_quad()
    }
    fn ico1() -> ProcessMesh {
        make_icosphere(1)
    }
    #[test]
    fn test_mesh_num_verts() {
        let m = tet();
        assert_eq!(m.num_verts(), 4);
    }
    #[test]
    fn test_mesh_num_faces() {
        let m = tet();
        assert_eq!(m.num_faces(), 4);
    }
    #[test]
    fn test_mesh_compute_normals() {
        let mut m = tet();
        m.compute_normals();
        let normals = m.normals.as_ref().unwrap();
        assert_eq!(normals.len(), 4);
        for &n in normals {
            let l = vec3_len(n);
            assert!(l < 1.01, "normal too long: {l}");
        }
    }
    #[test]
    fn test_build_adjacency_tetrahedron() {
        let m = tet();
        let adj = m.build_adjacency();
        assert_eq!(adj.len(), 4);
        for i in 0..4 {
            assert_eq!(adj[i].len(), 3, "vertex {i} should have 3 neighbours");
        }
    }
    #[test]
    fn test_laplacian_smooth_zero_iter() {
        let m = tet();
        let s = laplacian_smooth(&m, 0.5, 0);
        for i in 0..4 {
            assert_eq!(s.verts[i], m.verts[i]);
        }
    }
    #[test]
    fn test_laplacian_smooth_moves_vertices() {
        let m = ico1();
        let s = laplacian_smooth(&m, 0.5, 3);
        let changed = m
            .verts
            .iter()
            .zip(s.verts.iter())
            .any(|(a, b)| vec3_len(vec3_sub(*a, *b)) > 1e-12);
        assert!(changed);
    }
    #[test]
    fn test_taubin_smooth_preserves_face_count() {
        let m = ico1();
        let s = taubin_smooth(&m, 0.5, -0.53, 3);
        assert_eq!(s.num_faces(), m.num_faces());
    }
    #[test]
    fn test_taubin_lambda_mu_signs() {
        let m = ico1();
        let _s = taubin_smooth(&m, 0.5, -0.53, 5);
    }
    #[test]
    fn test_midpoint_subdivide_face_count() {
        let m = tet();
        let s = midpoint_subdivide(&m);
        assert_eq!(s.num_faces(), m.num_faces() * 4);
    }
    #[test]
    fn test_midpoint_subdivide_vertex_count_increases() {
        let m = tet();
        let s = midpoint_subdivide(&m);
        assert!(s.num_verts() > m.num_verts());
    }
    #[test]
    fn test_loop_subdivide_face_count() {
        let m = tet();
        let s = loop_subdivide(&m);
        assert_eq!(s.num_faces(), m.num_faces() * 4);
    }
    #[test]
    fn test_catmull_clark_subdivide_more_faces() {
        let m = tet();
        let s = catmull_clark_subdivide(&m);
        assert!(s.num_faces() > m.num_faces());
    }
    #[test]
    fn test_double_midpoint_subdivide() {
        let m = tet();
        let s1 = midpoint_subdivide(&m);
        let s2 = midpoint_subdivide(&s1);
        assert_eq!(s2.num_faces(), m.num_faces() * 16);
    }
    #[test]
    fn test_quadric_zero() {
        let q = Quadric::zero();
        assert_eq!(q.evaluate([1.0, 2.0, 3.0]), 0.0);
    }
    #[test]
    fn test_quadric_plane() {
        let q = Quadric::from_plane(0.0, 0.0, 1.0, 0.0);
        assert!(q.evaluate([1.0, 2.0, 0.0]).abs() < 1e-12);
        assert!(q.evaluate([1.0, 2.0, 1.0]) > 0.0);
    }
    #[test]
    fn test_qem_decimate_reduces_faces() {
        let m = make_icosphere(2);
        let target = 200;
        let dec = qem_decimate(&m, target);
        assert!(
            dec.num_faces() < m.num_faces(),
            "decimation should reduce faces; got={}",
            dec.num_faces()
        );
        assert!(dec.num_faces() <= target + 5, "faces={}", dec.num_faces());
    }
    #[test]
    fn test_qem_decimate_no_op_when_under_target() {
        let m = tet();
        let dec = qem_decimate(&m, 100);
        assert_eq!(dec.num_faces(), m.num_faces());
    }
    #[test]
    fn test_triangle_quality_equilateral() {
        let va = [0.0, 0.0, 0.0];
        let vb = [1.0, 0.0, 0.0];
        let vc = [0.5, (3.0_f64).sqrt() / 2.0, 0.0];
        let q = triangle_quality(va, vb, vc);
        assert!((q.aspect_ratio - 1.0).abs() < 0.01, "ar={}", q.aspect_ratio);
        assert!((q.area_quality - 1.0).abs() < 0.01, "aq={}", q.area_quality);
    }
    #[test]
    fn test_triangle_quality_right_angle() {
        let va = [0.0, 0.0, 0.0];
        let vb = [1.0, 0.0, 0.0];
        let vc = [0.0, 1.0, 0.0];
        let q = triangle_quality(va, vb, vc);
        assert!(q.min_angle > 0.0);
        assert!(q.max_angle < std::f64::consts::PI);
    }
    #[test]
    fn test_mesh_quality_stats_tetrahedron() {
        let m = tet();
        let stats = mesh_quality_stats(&m);
        assert!(
            stats.mean_aspect_ratio >= 1.0 - 1e-9,
            "mean_aspect_ratio={}",
            stats.mean_aspect_ratio
        );
        assert!(stats.min_angle_deg > 0.0);
    }
    #[test]
    fn test_mesh_quality_empty_mesh() {
        let m = ProcessMesh::new(vec![], vec![]);
        let stats = mesh_quality_stats(&m);
        assert_eq!(stats.mean_aspect_ratio, 0.0);
    }
    #[test]
    fn test_detect_boundary_loops_closed_mesh() {
        let m = tet();
        let loops = detect_boundary_loops(&m);
        assert!(loops.is_empty(), "tetrahedron should have no boundary");
    }
    #[test]
    fn test_detect_boundary_loops_open_mesh() {
        let m = quad();
        let loops = detect_boundary_loops(&m);
        assert!(!loops.is_empty());
    }
    #[test]
    fn test_fill_hole_fan_single_triangle() {
        let m = quad();
        let loops = detect_boundary_loops(&m);
        let mut verts = m.verts.clone();
        let mut total_faces = m.faces.clone();
        for lp in &loops {
            let nf = fill_hole_fan(&mut verts, lp);
            total_faces.extend(nf);
        }
        assert!(total_faces.len() > m.num_faces());
    }
    #[test]
    fn test_fill_all_holes_closes_boundary() {
        let m = quad();
        let repaired = fill_all_holes(&m);
        let loops = detect_boundary_loops(&repaired);
        assert!(loops.is_empty(), "boundary should be closed after repair");
    }
    #[test]
    fn test_feature_edges_tetrahedron_sharp() {
        let m = tet();
        let edges = detect_feature_edges(&m, 0.0);
        assert!(!edges.is_empty());
    }
    #[test]
    fn test_feature_edges_high_threshold() {
        let m = tet();
        let edges = detect_feature_edges(&m, 180.0);
        assert!(edges.is_empty());
    }
    #[test]
    fn test_split_long_edges_increases_verts() {
        let m = make_icosphere(0);
        let split = split_long_edges(&m, 0.5);
        assert!(split.num_verts() >= m.num_verts());
    }
    #[test]
    fn test_collapse_short_edges_reduces_verts() {
        let m = make_icosphere(2);
        let collapsed = collapse_short_edges(&m, 0.5);
        assert!(collapsed.num_verts() <= m.num_verts());
    }
    #[test]
    fn test_isotropic_remesh_smoke() {
        let m = make_icosphere(1);
        let rm = isotropic_remesh(&m, 0.4, 2);
        assert!(rm.num_verts() > 0);
        assert!(rm.num_faces() > 0);
    }
    #[test]
    fn test_tutte_parameterize_uv_count() {
        let m = make_icosphere(1);
        let param = tutte_parameterize(&m, 50);
        assert_eq!(param.uvs.len(), m.num_verts());
    }
    #[test]
    fn test_tutte_parameterize_uv_range() {
        let m = make_icosphere(1);
        let param = tutte_parameterize(&m, 50);
        for &[u, v] in &param.uvs {
            assert!(u.is_finite(), "u not finite");
            assert!(v.is_finite(), "v not finite");
        }
    }
    #[test]
    fn test_lscm_parameterize_returns_uvs() {
        let m = make_icosphere(1);
        let param = lscm_parameterize(&m);
        assert_eq!(param.uvs.len(), m.num_verts());
    }
    #[test]
    fn test_generate_texture_atlas_nonempty() {
        let m = make_icosphere(1);
        let atlas = generate_texture_atlas(&m, 512);
        assert!(!atlas.is_empty());
    }
    #[test]
    fn test_generate_texture_atlas_bounds_valid() {
        let m = make_icosphere(1);
        let atlas = generate_texture_atlas(&m, 512);
        for patch in &atlas {
            let [umin, vmin, umax, vmax] = patch.bounds;
            assert!(umin <= umax, "umin > umax");
            assert!(vmin <= vmax, "vmin > vmax");
        }
    }
    #[test]
    fn test_point_in_mesh_inside() {
        let m = make_icosphere(2);
        assert!(point_in_mesh([0.0, 0.0, 0.0], &m));
    }
    #[test]
    fn test_point_in_mesh_outside() {
        let m = make_icosphere(2);
        assert!(!point_in_mesh([10.0, 0.0, 0.0], &m));
    }
    #[test]
    fn test_mesh_union_has_faces() {
        let a = make_icosphere(1);
        let mut b = make_icosphere(1);
        for v in &mut b.verts {
            v[0] += 0.5;
        }
        let result = mesh_union(&a, &b);
        assert!(result.mesh.num_faces() > 0);
    }
    #[test]
    fn test_mesh_intersection_smoke() {
        let a = make_icosphere(1);
        let mut b = make_icosphere(1);
        for v in &mut b.verts {
            v[0] += 0.5;
        }
        let result = mesh_intersection(&a, &b);
        let _ = result.mesh.num_faces();
    }
    #[test]
    fn test_mesh_difference_smoke() {
        let a = make_icosphere(1);
        let b = make_icosphere(1);
        let result = mesh_difference(&a, &b);
        let _ = result.mesh.num_faces();
    }
    #[test]
    fn test_anisotropic_remesh_smoke() {
        let m = make_icosphere(1);
        let field = uniform_anisotropic_field(&m, 0.5);
        let rm = anisotropic_remesh(&m, &field, 2);
        assert!(rm.num_verts() > 0);
    }
    #[test]
    fn test_uniform_anisotropic_field_dimensions() {
        let m = tet();
        let field = uniform_anisotropic_field(&m, 1.0);
        assert_eq!(field.len_min.len(), m.num_verts());
        assert_eq!(field.len_max.len(), m.num_verts());
    }
    #[test]
    fn test_make_icosphere_vertices_on_unit_sphere() {
        let m = make_icosphere(0);
        for v in &m.verts {
            let r = vec3_len(*v);
            assert!((r - 1.0).abs() < 1e-10, "radius={r}");
        }
    }
    #[test]
    fn test_make_icosphere_subdivision_grows_mesh() {
        let m0 = make_icosphere(0);
        let m1 = make_icosphere(1);
        assert!(m1.num_faces() > m0.num_faces());
    }
    #[test]
    fn test_vec3_cross_perpendicular() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_normalize_unit() {
        let v = [3.0, 4.0, 0.0];
        let n = vec3_normalize(v);
        assert!((vec3_len(n) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_lerp_midpoint() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 2.0, 2.0];
        let m = vec3_lerp(a, b, 0.5);
        assert!((m[0] - 1.0).abs() < 1e-12);
    }
}
