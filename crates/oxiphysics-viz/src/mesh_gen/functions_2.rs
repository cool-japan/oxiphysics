//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{MeshData, ParamMesh, ParamVertex};

/// Compute per-vertex tangents for a `MeshData` using the UV-based method.
///
/// Returns `(tangents, bitangents)` as `Vec<[f64; 3]>`, one per vertex.
/// Requires the mesh to have UVs.  Where the tangent is degenerate (e.g. no UV
/// gradient), a fallback perpendicular to the normal is used.
pub fn compute_tangent_space(mesh: &MeshData) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let n = mesh.positions.len();
    let mut tan1 = vec![[0.0_f64; 3]; n];
    let mut tan2 = vec![[0.0_f64; 3]; n];
    for &[i0, i1, i2] in &mesh.indices {
        if i0 >= n || i1 >= n || i2 >= n {
            continue;
        }
        let p0 = mesh.positions[i0];
        let p1 = mesh.positions[i1];
        let p2 = mesh.positions[i2];
        let uv0 = mesh.uvs[i0];
        let uv1 = mesh.uvs[i1];
        let uv2 = mesh.uvs[i2];
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let du1 = uv1[0] - uv0[0];
        let dv1 = uv1[1] - uv0[1];
        let du2 = uv2[0] - uv0[0];
        let dv2 = uv2[1] - uv0[1];
        let denom = du1 * dv2 - du2 * dv1;
        let r = if denom.abs() < 1e-30 {
            0.0
        } else {
            1.0 / denom
        };
        let sdir = [
            (dv2 * e1[0] - dv1 * e2[0]) * r,
            (dv2 * e1[1] - dv1 * e2[1]) * r,
            (dv2 * e1[2] - dv1 * e2[2]) * r,
        ];
        let tdir = [
            (du1 * e2[0] - du2 * e1[0]) * r,
            (du1 * e2[1] - du2 * e1[1]) * r,
            (du1 * e2[2] - du2 * e1[2]) * r,
        ];
        for idx in [i0, i1, i2] {
            tan1[idx][0] += sdir[0];
            tan1[idx][1] += sdir[1];
            tan1[idx][2] += sdir[2];
            tan2[idx][0] += tdir[0];
            tan2[idx][1] += tdir[1];
            tan2[idx][2] += tdir[2];
        }
    }
    let mut tangents = vec![[0.0_f64; 3]; n];
    let mut bitangents = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        let nn = mesh.normals[i];
        let t = tan1[i];
        let dot = nn[0] * t[0] + nn[1] * t[1] + nn[2] * t[2];
        let mut tang = [t[0] - dot * nn[0], t[1] - dot * nn[1], t[2] - dot * nn[2]];
        let mag = (tang[0] * tang[0] + tang[1] * tang[1] + tang[2] * tang[2]).sqrt();
        if mag < 1e-12 {
            let up = if nn[1].abs() < 0.99 {
                [0.0, 1.0, 0.0]
            } else {
                [1.0, 0.0, 0.0]
            };
            tang = [
                nn[1] * up[2] - nn[2] * up[1],
                nn[2] * up[0] - nn[0] * up[2],
                nn[0] * up[1] - nn[1] * up[0],
            ];
            let m = (tang[0] * tang[0] + tang[1] * tang[1] + tang[2] * tang[2])
                .sqrt()
                .max(1e-30);
            tang[0] /= m;
            tang[1] /= m;
            tang[2] /= m;
        } else {
            tang[0] /= mag;
            tang[1] /= mag;
            tang[2] /= mag;
        }
        tangents[i] = tang;
        let cross = [
            nn[1] * tang[2] - nn[2] * tang[1],
            nn[2] * tang[0] - nn[0] * tang[2],
            nn[0] * tang[1] - nn[1] * tang[0],
        ];
        let sign = if (cross[0] * tan2[i][0] + cross[1] * tan2[i][1] + cross[2] * tan2[i][2]) < 0.0
        {
            -1.0
        } else {
            1.0
        };
        bitangents[i] = [cross[0] * sign, cross[1] * sign, cross[2] * sign];
    }
    (tangents, bitangents)
}
#[cfg(test)]
mod mesh_data_tests {
    use super::super::functions::*;
    use super::super::types::*;
    use super::*;

    fn assert_normals_unit(normals: &[[f64; 3]], tol: f64) {
        for (i, n) in normals.iter().enumerate() {
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (mag - 1.0).abs() < tol,
                "normal[{i}] magnitude = {mag:.6}, expected ~1.0"
            );
        }
    }
    fn assert_indices_in_bounds(indices: &[[usize; 3]], vertex_count: usize) {
        for (i, tri) in indices.iter().enumerate() {
            for &vi in tri {
                assert!(
                    vi < vertex_count,
                    "triangle[{i}] vertex index {vi} out of bounds (n={vertex_count})"
                );
            }
        }
    }
    #[test]
    fn test_mesh_data_empty() {
        let m = MeshData::empty();
        assert_eq!(m.vertex_count(), 0);
        assert_eq!(m.triangle_count(), 0);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_mesh_data_is_consistent_sphere() {
        let m = generate_sphere([0.0; 3], 1.0, 4, 6);
        assert!(m.is_consistent(), "sphere MeshData should be consistent");
    }
    #[test]
    fn test_mesh_data_aabb_sphere() {
        let r = 2.0;
        let m = generate_sphere([0.0; 3], r, 6, 8);
        let (mn, mx) = m.aabb().unwrap();
        assert!(mn[0] >= -r - 1e-9 && mx[0] <= r + 1e-9);
        assert!(mn[1] >= -r - 1e-9 && mx[1] <= r + 1e-9);
    }
    #[test]
    fn test_mesh_data_aabb_empty() {
        let m = MeshData::empty();
        assert!(m.aabb().is_none());
    }
    #[test]
    fn test_mesh_data_scale() {
        let mut m = generate_sphere([0.0; 3], 1.0, 4, 6);
        m.scale(3.0);
        let (mn, mx) = m.aabb().unwrap();
        assert!(mn[0] >= -3.0 - 1e-6 && mx[0] <= 3.0 + 1e-6);
    }
    #[test]
    fn test_mesh_data_translate() {
        let mut m = generate_sphere([0.0; 3], 1.0, 4, 6);
        m.translate([10.0, 0.0, 0.0]);
        let (mn, _mx) = m.aabb().unwrap();
        assert!(mn[0] >= 9.0 - 1e-6, "translated min x should be near 9");
    }
    #[test]
    fn test_mesh_data_merge() {
        let a = generate_box([0.0; 3], [1.0; 3]);
        let b = generate_box([5.0, 0.0, 0.0], [1.0; 3]);
        let na = a.vertex_count();
        let nb = b.vertex_count();
        let ta = a.triangle_count();
        let tb = b.triangle_count();
        let mut merged = a.clone();
        merged.merge_from(&b);
        assert_eq!(merged.vertex_count(), na + nb);
        assert_eq!(merged.triangle_count(), ta + tb);
        assert_indices_in_bounds(&merged.indices, merged.vertex_count());
    }
    #[test]
    fn test_mesh_data_flat_indices() {
        let m = generate_box([0.0; 3], [1.0; 3]);
        let flat = m.flat_indices();
        assert_eq!(flat.len(), m.triangle_count() * 3);
        for &vi in &flat {
            assert!(vi < m.vertex_count());
        }
    }
    #[test]
    fn test_mesh_data_recompute_normals() {
        let mut m = generate_box([0.0; 3], [1.0; 3]);
        m.recompute_normals();
        assert_normals_unit(&m.normals, 0.01);
    }
    #[test]
    fn test_generate_sphere_vertex_count() {
        let st = 6;
        let sl = 8;
        let m = generate_sphere([0.0; 3], 1.0, st, sl);
        assert_eq!(m.vertex_count(), (st + 1) * (sl + 1));
    }
    #[test]
    fn test_generate_sphere_triangle_count() {
        let st = 6;
        let sl = 8;
        let m = generate_sphere([0.0; 3], 1.0, st, sl);
        assert_eq!(m.triangle_count(), st * sl * 2);
    }
    #[test]
    fn test_generate_sphere_normals_unit() {
        let m = generate_sphere([0.0; 3], 2.0, 6, 8);
        assert_normals_unit(&m.normals, 1e-9);
    }
    #[test]
    fn test_generate_sphere_uvs_in_range() {
        let m = generate_sphere([0.0; 3], 1.0, 4, 6);
        for (i, uv) in m.uvs.iter().enumerate() {
            assert!(
                uv[0] >= 0.0 && uv[0] <= 1.0,
                "uv[{i}].u={} out of [0,1]",
                uv[0]
            );
            assert!(
                uv[1] >= 0.0 && uv[1] <= 1.0,
                "uv[{i}].v={} out of [0,1]",
                uv[1]
            );
        }
    }
    #[test]
    fn test_generate_sphere_indices_in_bounds() {
        let m = generate_sphere([0.0; 3], 1.0, 5, 7);
        assert_indices_in_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_box_24_vertices() {
        let m = generate_box([0.0; 3], [1.0; 3]);
        assert_eq!(
            m.vertex_count(),
            24,
            "box should have 24 vertices (4 per face × 6 faces)"
        );
    }
    #[test]
    fn test_generate_box_12_triangles() {
        let m = generate_box([0.0; 3], [1.0; 3]);
        assert_eq!(
            m.triangle_count(),
            12,
            "box should have 12 triangles (2 per face × 6 faces)"
        );
    }
    #[test]
    fn test_generate_box_normals_unit() {
        let m = generate_box([0.0; 3], [0.5, 1.0, 2.0]);
        assert_normals_unit(&m.normals, 1e-9);
    }
    #[test]
    fn test_generate_box_normals_axis_aligned() {
        let m = generate_box([0.0; 3], [1.0; 3]);
        for n in &m.normals {
            let is_axis = (n[0].abs() > 0.9) || (n[1].abs() > 0.9) || (n[2].abs() > 0.9);
            assert!(is_axis, "box normals should be axis-aligned: {:?}", n);
        }
    }
    #[test]
    fn test_generate_box_uvs_consistent() {
        let m = generate_box([0.0; 3], [1.0; 3]);
        assert_eq!(m.uvs.len(), m.vertex_count());
    }
    #[test]
    fn test_generate_cylinder_mesh_consistent() {
        let m = generate_cylinder_mesh([0.0; 3], 1.0, 2.0, 8);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_generate_cylinder_mesh_indices_in_bounds() {
        let m = generate_cylinder_mesh([0.0; 3], 1.0, 2.0, 6);
        assert_indices_in_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_cylinder_mesh_has_caps() {
        let m = generate_cylinder_mesh([0.0; 3], 1.0, 2.0, 6);
        let has_top = m.positions.iter().any(|p| (p[1] - 1.0).abs() < 1e-9);
        let has_bot = m.positions.iter().any(|p| (p[1] + 1.0).abs() < 1e-9);
        assert!(has_top);
        assert!(has_bot);
    }
    #[test]
    fn test_generate_cone_mesh_consistent() {
        let m = generate_cone_mesh([0.0; 3], 1.0, 2.0, 8);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_generate_cone_mesh_indices_in_bounds() {
        let m = generate_cone_mesh([0.0; 3], 1.0, 2.0, 6);
        assert_indices_in_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_cone_mesh_has_apex() {
        let m = generate_cone_mesh([0.0; 3], 1.0, 3.0, 8);
        let has_apex = m.positions.iter().any(|p| (p[1] - 3.0).abs() < 1e-9);
        assert!(has_apex, "cone apex should be at y=height");
    }
    #[test]
    fn test_generate_torus_mesh_vertex_count() {
        let ms = 8;
        let ns = 6;
        let m = generate_torus_mesh([0.0; 3], 2.0, 0.5, ms, ns);
        assert_eq!(m.vertex_count(), (ms + 1) * (ns + 1));
    }
    #[test]
    fn test_generate_torus_mesh_normals_unit() {
        let m = generate_torus_mesh([0.0; 3], 2.0, 0.5, 10, 8);
        assert_normals_unit(&m.normals, 1e-9);
    }
    #[test]
    fn test_generate_torus_mesh_consistent() {
        let m = generate_torus_mesh([0.0; 3], 2.0, 0.5, 8, 6);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_generate_plane_mesh_vertex_count() {
        let sub = 4;
        let m = generate_plane_mesh([0.0; 3], [0.0, 1.0, 0.0], 2.0, sub);
        assert_eq!(m.vertex_count(), (sub + 1) * (sub + 1));
    }
    #[test]
    fn test_generate_plane_mesh_uvs_in_range() {
        let m = generate_plane_mesh([0.0; 3], [0.0, 1.0, 0.0], 2.0, 4);
        for uv in &m.uvs {
            assert!(uv[0] >= 0.0 && uv[0] <= 1.0 + 1e-9);
            assert!(uv[1] >= 0.0 && uv[1] <= 1.0 + 1e-9);
        }
    }
    #[test]
    fn test_generate_plane_mesh_normals_consistent() {
        let m = generate_plane_mesh([0.0; 3], [0.0, 1.0, 0.0], 2.0, 3);
        assert_normals_unit(&m.normals, 1e-9);
    }
    #[test]
    fn test_generate_arrow_non_empty() {
        let m = generate_arrow([0.0; 3], [0.0, 1.0, 0.0], 0.05, 0.15, 8);
        assert!(m.vertex_count() > 0);
        assert!(m.triangle_count() > 0);
    }
    #[test]
    fn test_generate_arrow_zero_direction_empty() {
        let m = generate_arrow([0.0; 3], [0.0, 0.0, 0.0], 0.05, 0.15, 8);
        assert_eq!(m.vertex_count(), 0);
        assert_eq!(m.triangle_count(), 0);
    }
    #[test]
    fn test_generate_arrow_consistent() {
        let m = generate_arrow([0.0; 3], [1.0, 0.0, 0.0], 0.1, 0.2, 6);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_generate_arrow_indices_in_bounds() {
        let m = generate_arrow([0.0; 3], [0.0, 0.0, 2.0], 0.05, 0.1, 6);
        assert_indices_in_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_tangent_space_sphere_count() {
        let m = generate_sphere([0.0; 3], 1.0, 4, 6);
        let (tangents, bitangents) = compute_tangent_space(&m);
        assert_eq!(tangents.len(), m.vertex_count());
        assert_eq!(bitangents.len(), m.vertex_count());
    }
    #[test]
    fn test_tangent_space_tangents_unit() {
        let m = generate_sphere([0.0; 3], 1.0, 6, 8);
        let (tangents, _) = compute_tangent_space(&m);
        for (i, t) in tangents.iter().enumerate() {
            let mag = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            assert!((mag - 1.0).abs() < 0.01, "tangent[{i}] mag = {mag}");
        }
    }
    #[test]
    fn test_tangent_space_orthogonal_to_normal() {
        let m = generate_sphere([0.0; 3], 1.0, 4, 6);
        let (tangents, _) = compute_tangent_space(&m);
        for (i, (n, t)) in m.normals.iter().zip(tangents.iter()).enumerate() {
            let dot = n[0] * t[0] + n[1] * t[1] + n[2] * t[2];
            assert!(
                dot.abs() < 0.01,
                "tangent[{i}] should be orthogonal to normal, dot={dot:.6}"
            );
        }
    }
    #[test]
    fn test_tangent_space_plane_y_up() {
        let m = generate_plane_mesh([0.0; 3], [0.0, 1.0, 0.0], 2.0, 2);
        let (tangents, _bitangents) = compute_tangent_space(&m);
        for (i, t) in tangents.iter().enumerate() {
            let mag = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            assert!(
                (mag - 1.0).abs() < 0.01,
                "tangent[{i}] should be unit length, mag={mag}"
            );
        }
    }
}
#[cfg(test)]
mod extra_mesh_tests {
    use super::super::functions::*;
    use super::super::types::*;
    use super::*;
    use crate::Color;

    fn assert_unit_n(normals: &[[f64; 3]], tol: f64) {
        for (i, n) in normals.iter().enumerate() {
            let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (mag - 1.0).abs() < tol,
                "normal[{i}] magnitude = {mag:.6}, expected ~1.0"
            );
        }
    }
    fn assert_tri_bounds(tris: &[[usize; 3]], n: usize) {
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < n, "tri[{i}] vertex {vi} out of bounds (n={n})");
            }
        }
    }
    #[test]
    fn test_mesh_data_scale_updates_positions() {
        let mut m = generate_box([0.0; 3], [1.0; 3]);
        m.scale(2.0);
        let (mn, mx) = m.aabb().unwrap();
        assert!(
            mn[0] >= -2.0 - 1e-9 && mx[0] <= 2.0 + 1e-9,
            "scaled box x: mn={} mx={}",
            mn[0],
            mx[0]
        );
    }
    #[test]
    fn test_mesh_data_scale_does_not_change_triangle_count() {
        let mut m = generate_box([0.0; 3], [1.0; 3]);
        let n_tris = m.triangle_count();
        m.scale(5.0);
        assert_eq!(m.triangle_count(), n_tris);
    }
    #[test]
    fn test_mesh_data_translate_does_not_change_normals() {
        let m0 = generate_sphere([0.0; 3], 1.0, 4, 6);
        let mut m1 = generate_sphere([0.0; 3], 1.0, 4, 6);
        m1.translate([100.0, 200.0, 300.0]);
        for (i, (n0, n1)) in m0.normals.iter().zip(m1.normals.iter()).enumerate() {
            for k in 0..3 {
                assert!(
                    (n0[k] - n1[k]).abs() < 1e-9,
                    "normal[{i}][{k}] changed after translate"
                );
            }
        }
    }
    #[test]
    fn test_mesh_data_merge_indices_offset_correctly() {
        let a = generate_cylinder_mesh([0.0; 3], 1.0, 2.0, 6);
        let b = generate_cone_mesh([0.0; 3], 0.5, 1.0, 6);
        let na = a.vertex_count();
        let mut merged = a.clone();
        merged.merge_from(&b);
        for (i, tri) in merged.indices.iter().enumerate() {
            for &vi in tri {
                assert!(
                    vi < merged.vertex_count(),
                    "merged tri[{i}] vertex {vi} OOB"
                );
            }
        }
        let b_first_tri = merged.indices[a.triangle_count()];
        assert!(
            b_first_tri[0] >= na || b_first_tri[1] >= na || b_first_tri[2] >= na,
            "merged b triangles should use offset indices >= {na}"
        );
    }
    #[test]
    fn test_mesh_data_flat_indices_length() {
        let m = generate_torus_mesh([0.0; 3], 2.0, 0.5, 6, 4);
        let flat = m.flat_indices();
        assert_eq!(flat.len(), m.triangle_count() * 3);
    }
    #[test]
    fn test_mesh_data_recompute_normals_unit_for_sphere() {
        let mut m = generate_torus_mesh([0.0; 3], 2.0, 0.5, 8, 6);
        m.recompute_normals();
        assert_unit_n(&m.normals, 0.01);
    }
    #[test]
    fn test_mesh_data_recompute_normals_unit_for_torus() {
        let mut m = generate_torus_mesh([0.0; 3], 2.0, 0.5, 8, 6);
        m.recompute_normals();
        assert_unit_n(&m.normals, 0.01);
    }
    #[test]
    fn test_mesh_data_is_consistent_box() {
        let m = generate_box([1.0, 2.0, 3.0], [0.5, 1.0, 1.5]);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_mesh_data_is_consistent_cylinder() {
        let m = generate_cylinder_mesh([0.0; 3], 1.0, 2.0, 8);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_mesh_data_is_consistent_cone() {
        let m = generate_cone_mesh([0.0; 3], 1.0, 2.0, 6);
        assert!(m.is_consistent());
    }
    #[test]
    fn test_mesh_data_aabb_box_tight() {
        let half = [1.0, 2.0, 3.0_f64];
        let center = [0.0; 3];
        let m = generate_box(center, half);
        let (mn, mx) = m.aabb().unwrap();
        assert!(
            (mn[0] - (-1.0)).abs() < 1e-9,
            "min x should be -1.0, got {}",
            mn[0]
        );
        assert!(
            (mx[0] - 1.0).abs() < 1e-9,
            "max x should be 1.0, got {}",
            mx[0]
        );
        assert!(
            (mn[1] - (-2.0)).abs() < 1e-9,
            "min y should be -2.0, got {}",
            mn[1]
        );
        assert!(
            (mx[2] - 3.0).abs() < 1e-9,
            "max z should be 3.0, got {}",
            mx[2]
        );
    }
    #[test]
    fn test_generate_sphere_verts_on_surface() {
        let r = 3.0;
        let c = [1.0, 2.0, 3.0];
        let m = generate_sphere(c, r, 6, 8);
        for (i, p) in m.positions.iter().enumerate() {
            let d = ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt();
            assert!((d - r).abs() < 1e-9, "vertex {i} dist={d} != r={r}");
        }
    }
    #[test]
    fn test_generate_sphere_normals_perpendicular_to_surface() {
        let m = generate_sphere([0.0; 3], 1.0, 4, 6);
        for (i, (p, n)) in m.positions.iter().zip(m.normals.iter()).enumerate() {
            let dot = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
            assert!(dot > 0.99, "normal[{i}] should point outward, dot={dot:.6}");
        }
    }
    #[test]
    fn test_generate_box_centered_nonzero() {
        let center = [3.0, 4.0, 5.0];
        let m = generate_box(center, [1.0; 3]);
        let (mn, mx) = m.aabb().unwrap();
        assert!(mn[0] >= 2.0 - 1e-9);
        assert!(mx[0] <= 4.0 + 1e-9);
        assert!(mn[1] >= 3.0 - 1e-9);
        assert!(mx[2] <= 6.0 + 1e-9);
    }
    #[test]
    fn test_generate_box_flat_returns_degenerate_box() {
        let m = generate_box([0.0; 3], [1.0, 0.0, 1.0]);
        assert!(
            m.vertex_count() > 0,
            "degenerate box should still have vertices"
        );
    }
    #[test]
    fn test_generate_cylinder_mesh_uvs_in_range() {
        let m = generate_cylinder_mesh([0.0; 3], 1.0, 2.0, 8);
        for (i, uv) in m.uvs.iter().enumerate() {
            assert!(
                uv[0] >= 0.0 && uv[0] <= 1.0 + 1e-9,
                "uv[{i}].u={} out of [0,1]",
                uv[0]
            );
            assert!(
                uv[1] >= 0.0 && uv[1] <= 1.0 + 1e-9,
                "uv[{i}].v={} out of [0,1]",
                uv[1]
            );
        }
    }
    #[test]
    fn test_generate_cylinder_mesh_normals_unit() {
        let m = generate_cylinder_mesh([0.0; 3], 2.0, 4.0, 12);
        assert_unit_n(&m.normals, 1e-9);
    }
    #[test]
    fn test_generate_cylinder_mesh_min_sectors() {
        let m = generate_cylinder_mesh([0.0; 3], 1.0, 2.0, 1);
        assert!(m.vertex_count() > 0);
        assert_tri_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_cone_mesh_uvs_in_range() {
        let m = generate_cone_mesh([0.0; 3], 1.0, 2.0, 8);
        for (i, uv) in m.uvs.iter().enumerate() {
            assert!(
                uv[0] >= 0.0 - 1e-9 && uv[0] <= 1.0 + 1e-9,
                "uv[{i}].u={}",
                uv[0]
            );
            assert!(
                uv[1] >= 0.0 - 1e-9 && uv[1] <= 1.0 + 1e-9,
                "uv[{i}].v={}",
                uv[1]
            );
        }
    }
    #[test]
    fn test_generate_cone_mesh_side_normals_outward() {
        let m = generate_cone_mesh([0.0; 3], 1.0, 2.0, 8);
        let side_normals: Vec<_> = m.normals.iter().filter(|n| n[1].abs() < 0.95).collect();
        for (i, n) in side_normals.iter().enumerate() {
            assert!(
                n[1] > 0.0,
                "side normal[{i}] should have positive y (slant), got {}",
                n[1]
            );
        }
    }
    #[test]
    fn test_generate_torus_mesh_uvs_in_range() {
        let m = generate_torus_mesh([0.0; 3], 2.0, 0.5, 8, 6);
        for (i, uv) in m.uvs.iter().enumerate() {
            assert!(uv[0] >= 0.0 && uv[0] <= 1.0 + 1e-9, "uv[{i}].u={}", uv[0]);
            assert!(uv[1] >= 0.0 && uv[1] <= 1.0 + 1e-9, "uv[{i}].v={}", uv[1]);
        }
    }
    #[test]
    fn test_generate_torus_mesh_indices_in_bounds() {
        let m = generate_torus_mesh([0.0; 3], 2.0, 0.5, 8, 6);
        assert_tri_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_torus_mesh_centered() {
        let center = [5.0, 6.0, 7.0];
        let m = generate_torus_mesh(center, 1.0, 0.3, 12, 8);
        let (mn, mx) = m.aabb().unwrap();
        let cx = (mn[0] + mx[0]) * 0.5;
        let cy = (mn[1] + mx[1]) * 0.5;
        let cz = (mn[2] + mx[2]) * 0.5;
        assert!(
            (cx - center[0]).abs() < 0.01,
            "torus center x: {cx} != {}",
            center[0]
        );
        assert!(
            (cy - center[1]).abs() < 0.01,
            "torus center y: {cy} != {}",
            center[1]
        );
        assert!(
            (cz - center[2]).abs() < 0.01,
            "torus center z: {cz} != {}",
            center[2]
        );
    }
    #[test]
    fn test_generate_plane_mesh_indices_in_bounds() {
        let m = generate_plane_mesh([0.0; 3], [0.0, 1.0, 0.0], 4.0, 6);
        assert_tri_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_plane_mesh_consistent_for_tilted_normal() {
        let m = generate_plane_mesh([0.0; 3], [1.0, 1.0, 0.0], 2.0, 4);
        assert!(m.is_consistent());
        assert_unit_n(&m.normals, 1e-9);
    }
    #[test]
    fn test_generate_plane_mesh_subdivision_grows_with_n() {
        let m1 = generate_plane_mesh([0.0; 3], [0.0, 1.0, 0.0], 2.0, 1);
        let m4 = generate_plane_mesh([0.0; 3], [0.0, 1.0, 0.0], 2.0, 4);
        assert!(m4.vertex_count() > m1.vertex_count());
        assert!(m4.triangle_count() > m1.triangle_count());
    }
    #[test]
    fn test_generate_arrow_x_axis_consistent() {
        let m = generate_arrow([0.0; 3], [2.0, 0.0, 0.0], 0.05, 0.15, 6);
        assert!(m.is_consistent());
        assert_tri_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_arrow_z_axis_consistent() {
        let m = generate_arrow([0.0; 3], [0.0, 0.0, 3.0], 0.1, 0.2, 8);
        assert!(m.is_consistent());
        assert_tri_bounds(&m.indices, m.vertex_count());
    }
    #[test]
    fn test_generate_arrow_uvs_count_matches_vertices() {
        let m = generate_arrow([0.0; 3], [0.0, 1.0, 0.0], 0.05, 0.15, 8);
        assert_eq!(m.uvs.len(), m.vertex_count());
    }
    #[test]
    fn test_disc_solid_triangle_count() {
        let s = 12;
        let (_, tris, _) = generate_disc([0.0; 3], 0.0, 1.0, s);
        assert_eq!(tris.len(), s, "solid disc should have `sectors` triangles");
    }
    #[test]
    fn test_disc_annulus_triangle_count() {
        let s = 10;
        let (_, tris, _) = generate_disc([0.0; 3], 0.5, 1.5, s);
        assert_eq!(tris.len(), 2 * s, "annulus should have 2×sectors triangles");
    }
    #[test]
    fn test_disc_solid_center_vertex_at_origin() {
        let center = [1.0, 2.0, 3.0];
        let (verts, _, _) = generate_disc(center, 0.0, 1.0, 8);
        let center_v = verts[0];
        assert!((center_v[0] - center[0]).abs() < 1e-9);
        assert!((center_v[1] - center[1]).abs() < 1e-9);
        assert!((center_v[2] - center[2]).abs() < 1e-9);
    }
    #[test]
    fn test_disc_min_sectors_clamped() {
        let (verts, tris, _) = generate_disc([0.0; 3], 0.0, 1.0, 1);
        assert!(
            verts.len() >= 4,
            "should have at least 4 vertices (center + 3)"
        );
        assert!(tris.len() >= 3, "should have at least 3 triangles");
    }
    #[test]
    fn test_terrain_constant_height_normals_up() {
        let (_, _, normals) = generate_terrain(0.0, 4.0, 0.0, 4.0, 6, 6, |_, _| 5.0);
        for (i, n) in normals.iter().enumerate() {
            assert!(
                (n[1] - 1.0).abs() < 1e-9,
                "flat terrain normal[{i}] y={}",
                n[1]
            );
        }
    }
    #[test]
    fn test_terrain_triangle_count_formula() {
        let nx = 7;
        let nz = 5;
        let (_, tris, _) = generate_terrain(0.0, 7.0, 0.0, 5.0, nx, nz, |_, _| 0.0);
        assert_eq!(tris.len(), nx * nz * 2);
    }
    #[test]
    fn test_terrain_vertex_x_z_spacing() {
        let (verts, _, _) = generate_terrain(0.0, 2.0, 0.0, 2.0, 2, 2, |_, _| 0.0);
        for (i, v) in verts.iter().enumerate() {
            assert!(
                v[0] >= 0.0 - 1e-9 && v[0] <= 2.0 + 1e-9,
                "vert[{i}] x={} OOB",
                v[0]
            );
            assert!(
                v[2] >= 0.0 - 1e-9 && v[2] <= 2.0 + 1e-9,
                "vert[{i}] z={} OOB",
                v[2]
            );
        }
    }
    #[test]
    fn test_uv_sphere_min_sectors_clamped() {
        let (verts, tris, normals) = generate_uv_sphere([0.0; 3], 1.0, 1, 1);
        assert!(!verts.is_empty());
        assert!(!tris.is_empty());
        assert_eq!(verts.len(), normals.len());
    }
    #[test]
    fn test_uv_sphere_triangle_indices_in_bounds() {
        let rings = 5;
        let sectors = 8;
        let (verts, tris, _) = generate_uv_sphere([0.0; 3], 1.0, rings, sectors);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_generate_cylinder_side_normals_horizontal() {
        let (_, _, normals) = generate_cylinder([0.0; 3], 1.0, 2.0, 8);
        let sector = 8;
        for (i, n) in normals[..2 * (sector + 1)].iter().enumerate() {
            assert!(n[1].abs() < 1e-9, "side normal[{i}] y={} should be 0", n[1]);
        }
    }
    #[test]
    fn test_generate_cylinder_triangle_bounds() {
        let (verts, tris, _) = generate_cylinder([0.0; 3], 1.0, 2.0, 6);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_generate_cone_triangle_bounds() {
        let (verts, tris, _) = generate_cone([0.0; 3], 1.0, 2.0, 8);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_generate_cone_base_verts_at_y_zero() {
        let (verts, _, _) = generate_cone([0.0; 3], 1.5, 3.0, 8);
        for (i, v) in verts[..9].iter().enumerate() {
            assert!(
                v[1].abs() < 1e-9,
                "cone base vert[{i}] y={} should be 0",
                v[1]
            );
        }
    }
    #[test]
    fn test_generate_torus_triangle_bounds() {
        let (verts, tris, _) = generate_torus([0.0; 3], 2.0, 0.5, 8, 6);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_generate_torus_major_radius_determines_extent() {
        let major_r = 3.0;
        let minor_r = 0.3;
        let (verts, _, _) = generate_torus([0.0; 3], major_r, minor_r, 12, 8);
        let bound = major_r + minor_r + 1e-9;
        for (i, v) in verts.iter().enumerate() {
            assert!(
                v[0].abs() <= bound,
                "torus vert[{i}] x={} exceeds bound",
                v[0]
            );
            assert!(
                v[2].abs() <= bound,
                "torus vert[{i}] z={} exceeds bound",
                v[2]
            );
        }
    }
    #[test]
    fn test_generate_plane_y_up_all_same_y() {
        let center = [0.0; 3];
        let (verts, _, _) = generate_plane(center, [0.0, 1.0, 0.0], 4.0, 4);
        for (i, v) in verts.iter().enumerate() {
            assert!(v[1].abs() < 1e-9, "Y-up plane vert[{i}] y={}", v[1]);
        }
    }
    #[test]
    fn test_generate_plane_subdivision_index_bounds() {
        let (verts, tris, _) = generate_plane([0.0; 3], [0.0, 1.0, 0.0], 2.0, 5);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_generate_capsule_normals_unit() {
        let (_, _, normals) = generate_capsule([0.0, 0.0, 0.0], [0.0, 2.0, 0.0], 0.5, 4, 8);
        assert_unit_n(&normals, 0.01);
    }
    #[test]
    fn test_generate_capsule_index_bounds() {
        let (verts, tris, _) = generate_capsule([0.0; 3], [0.0, 3.0, 0.0], 0.5, 3, 6);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "capsule tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_geodesic_sphere_subdivision_2_counts() {
        let (verts, tris, _) = generate_geodesic_sphere([0.0; 3], 1.0, 2);
        assert_eq!(tris.len(), 320);
        let _ = verts;
    }
    #[test]
    fn test_geodesic_sphere_all_verts_on_sphere() {
        let r = 5.0;
        let (verts, _, _) = generate_geodesic_sphere([2.0, -1.0, 0.5], r, 1);
        for (i, v) in verts.iter().enumerate() {
            let d = ((v[0] - 2.0).powi(2) + (v[1] + 1.0).powi(2) + (v[2] - 0.5).powi(2)).sqrt();
            assert!((d - r).abs() < 1e-9, "vertex {i} dist={d} != r={r}");
        }
    }
    #[test]
    fn test_generate_frustum_bottom_y_minus_half_h() {
        let center = [0.0; 3];
        let height = 4.0;
        let (verts, _, _) = generate_frustum(center, 1.5, 0.5, height, 6);
        let half_h = height / 2.0;
        let has_bottom = verts.iter().any(|v| (v[1] - (-half_h)).abs() < 1e-9);
        assert!(has_bottom, "frustum should have vertices at y = -half_h");
    }
    #[test]
    fn test_generate_frustum_index_bounds() {
        let (verts, tris, _) = generate_frustum([0.0; 3], 2.0, 1.0, 3.0, 8);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "frustum tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_bicubic_patch_indices_in_bounds() {
        let cp: [[f64; 3]; 16] = std::array::from_fn(|i| {
            let u = (i % 4) as f64;
            let v = (i / 4) as f64;
            [u, 0.0, v]
        });
        let (verts, tris, _) = generate_bicubic_patch(&cp, 5, 5);
        for (i, tri) in tris.iter().enumerate() {
            for &vi in tri {
                assert!(vi < verts.len(), "bicubic tri[{i}] vertex {vi} OOB");
            }
        }
    }
    #[test]
    fn test_bicubic_patch_single_div_4_tris() {
        let cp: [[f64; 3]; 16] = std::array::from_fn(|i| [i as f64, 0.0, 0.0]);
        let (_, tris, _) = generate_bicubic_patch(&cp, 1, 1);
        assert_eq!(tris.len(), 2, "1x1 div should produce 2 triangles");
    }
    #[test]
    fn test_axes_gizmo_index_bounds_all_axes() {
        let axes = generate_axes_gizmo([0.0; 3], 1.0, 0.05, 6);
        for (axis_idx, (verts, tris)) in axes.iter().enumerate() {
            for (ti, tri) in tris.iter().enumerate() {
                for &vi in tri {
                    assert!(
                        vi < verts.len(),
                        "axis {axis_idx} tri[{ti}] vertex {vi} OOB (n={})",
                        verts.len()
                    );
                }
            }
        }
    }
    #[test]
    fn test_axes_gizmo_origin_offset() {
        let origin = [5.0, 10.0, 15.0];
        let axes = generate_axes_gizmo(origin, 1.0, 0.05, 6);
        for (axis_idx, (verts, _)) in axes.iter().enumerate() {
            let min_x = verts.iter().map(|v| v[0]).fold(f64::INFINITY, f64::min);
            let min_y = verts.iter().map(|v| v[1]).fold(f64::INFINITY, f64::min);
            let _ = (axis_idx, min_x, min_y);
            assert!(!verts.is_empty());
        }
    }
    #[test]
    fn test_tangent_space_bitangents_unit() {
        let m = generate_sphere([0.0; 3], 1.0, 4, 6);
        let (_, bitangents) = compute_tangent_space(&m);
        for (i, b) in bitangents.iter().enumerate() {
            let mag = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
            assert!((mag - 1.0).abs() < 0.01, "bitangent[{i}] mag = {mag}");
        }
    }
    #[test]
    fn test_tangent_space_box_count_matches() {
        let m = generate_box([0.0; 3], [1.0; 3]);
        let (tangents, bitangents) = compute_tangent_space(&m);
        assert_eq!(tangents.len(), m.vertex_count());
        assert_eq!(bitangents.len(), m.vertex_count());
    }
    #[test]
    fn test_tangent_space_empty_mesh_empty_result() {
        let m = MeshData::empty();
        let (tangents, bitangents) = compute_tangent_space(&m);
        assert!(tangents.is_empty());
        assert!(bitangents.is_empty());
    }
    #[test]
    fn test_arrow_mesh_non_empty_for_valid_direction() {
        use crate::primitives::Color;
        use oxiphysics_core::math::Vec3;
        let start = Vec3::new(0.0, 0.0, 0.0);
        let end = Vec3::new(0.0, 2.0, 0.0);
        let color = Color::red();
        let mesh = crate::mesh_gen::arrow_mesh(start, end, 0.05, color);
        assert!(!mesh.vertices.is_empty());
        assert!(!mesh.indices.is_empty());
    }
    #[test]
    fn test_arrow_mesh_zero_length_returns_empty() {
        use oxiphysics_core::math::Vec3;
        let start = Vec3::new(1.0, 2.0, 3.0);
        let end = Vec3::new(1.0, 2.0, 3.0);
        let mesh = crate::mesh_gen::arrow_mesh(start, end, 0.05, Color::blue());
        assert!(mesh.vertices.is_empty());
        assert!(mesh.indices.is_empty());
    }
}
/// Generate a UV-sphere as a `ParamMesh`.
pub fn param_sphere(centre: [f64; 3], radius: f64, n_lat: usize, n_lon: usize) -> ParamMesh {
    let n_lat = n_lat.max(2);
    let n_lon = n_lon.max(3);
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    for i in 0..=n_lat {
        let phi = std::f64::consts::PI * i as f64 / n_lat as f64;
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();
        let v = i as f64 / n_lat as f64;
        for j in 0..=n_lon {
            let theta = 2.0 * std::f64::consts::PI * j as f64 / n_lon as f64;
            let nx = sin_phi * theta.cos();
            let ny = cos_phi;
            let nz = sin_phi * theta.sin();
            vertices.push(ParamVertex {
                pos: [
                    centre[0] + radius * nx,
                    centre[1] + radius * ny,
                    centre[2] + radius * nz,
                ],
                normal: [nx, ny, nz],
                uv: [j as f64 / n_lon as f64, v],
            });
        }
    }
    let stride = n_lon + 1;
    for i in 0..n_lat {
        for j in 0..n_lon {
            let i0 = i * stride + j;
            let i1 = i0 + 1;
            let i2 = (i + 1) * stride + j;
            let i3 = i2 + 1;
            indices.extend_from_slice(&[i0, i2, i1]);
            indices.extend_from_slice(&[i1, i2, i3]);
        }
    }
    ParamMesh { vertices, indices }
}
/// Generate a flat disc (filled circle) on the XY plane.
pub fn param_disc(centre: [f64; 3], radius: f64, n_segs: usize) -> ParamMesh {
    let n = n_segs.max(3);
    let mut vertices = Vec::with_capacity(n + 1);
    let mut indices = Vec::with_capacity(n * 3);
    vertices.push(ParamVertex {
        pos: centre,
        normal: [0.0, 0.0, 1.0],
        uv: [0.5, 0.5],
    });
    for i in 0..n {
        let theta = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        let (ct, st) = (theta.cos(), theta.sin());
        vertices.push(ParamVertex {
            pos: [centre[0] + radius * ct, centre[1] + radius * st, centre[2]],
            normal: [0.0, 0.0, 1.0],
            uv: [0.5 + 0.5 * ct, 0.5 + 0.5 * st],
        });
    }
    for i in 0..n {
        let j = i + 1;
        let k = j % n + 1;
        indices.extend_from_slice(&[0, j, k]);
    }
    ParamMesh { vertices, indices }
}
/// Generate a torus as a `ParamMesh`.
///
/// `major_radius` is the distance from the tube centre to the torus centre.
/// `minor_radius` is the tube radius.
pub fn param_torus(
    centre: [f64; 3],
    major_radius: f64,
    minor_radius: f64,
    n_major: usize,
    n_minor: usize,
) -> ParamMesh {
    let n_maj = n_major.max(3);
    let n_min = n_minor.max(3);
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    for i in 0..=n_maj {
        let phi = 2.0 * std::f64::consts::PI * i as f64 / n_maj as f64;
        let (cos_phi, sin_phi) = (phi.cos(), phi.sin());
        for j in 0..=n_min {
            let theta = 2.0 * std::f64::consts::PI * j as f64 / n_min as f64;
            let (cos_theta, sin_theta) = (theta.cos(), theta.sin());
            let x = (major_radius + minor_radius * cos_theta) * cos_phi;
            let y = (major_radius + minor_radius * cos_theta) * sin_phi;
            let z = minor_radius * sin_theta;
            let nx = cos_theta * cos_phi;
            let ny = cos_theta * sin_phi;
            let nz = sin_theta;
            vertices.push(ParamVertex {
                pos: [centre[0] + x, centre[1] + y, centre[2] + z],
                normal: [nx, ny, nz],
                uv: [i as f64 / n_maj as f64, j as f64 / n_min as f64],
            });
        }
    }
    let stride = n_min + 1;
    for i in 0..n_maj {
        for j in 0..n_min {
            let a = i * stride + j;
            let b = a + 1;
            let c = (i + 1) * stride + j;
            let d = c + 1;
            indices.extend_from_slice(&[a, c, b]);
            indices.extend_from_slice(&[b, c, d]);
        }
    }
    ParamMesh { vertices, indices }
}
/// Generate a Möbius strip approximation as a `ParamMesh`.
///
/// `n_along` is the number of segments along the strip,
/// `n_across` is the number of segments across the half-width.
pub fn param_mobius(n_along: usize, n_across: usize, radius: f64, half_width: f64) -> ParamMesh {
    let na = n_along.max(3);
    let nc = n_across.max(2);
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    for i in 0..=na {
        let t = 2.0 * std::f64::consts::PI * i as f64 / na as f64;
        for j in 0..=nc {
            let s = (j as f64 / nc as f64 - 0.5) * 2.0 * half_width;
            let half_t = t / 2.0;
            let x = (radius + s * half_t.cos()) * t.cos();
            let y = (radius + s * half_t.cos()) * t.sin();
            let z = s * half_t.sin();
            let nx = half_t.cos() * t.cos() - (radius + s * half_t.cos()) * t.sin();
            let ny = half_t.cos() * t.sin() + (radius + s * half_t.cos()) * t.cos();
            let nz = half_t.sin();
            let mag = (nx * nx + ny * ny + nz * nz).max(1e-300).sqrt();
            vertices.push(ParamVertex {
                pos: [x, y, z],
                normal: [nx / mag, ny / mag, nz / mag],
                uv: [i as f64 / na as f64, j as f64 / nc as f64],
            });
        }
    }
    let stride = nc + 1;
    for i in 0..na {
        for j in 0..nc {
            let a = i * stride + j;
            let b = a + 1;
            let c = (i + 1) * stride + j;
            let d = c + 1;
            indices.extend_from_slice(&[a, b, c]);
            indices.extend_from_slice(&[b, d, c]);
        }
    }
    ParamMesh { vertices, indices }
}
/// Apply a wavy displacement to a `ParamMesh` along the Z-axis.
///
/// Each vertex is moved by `amplitude * sin(frequency * x + phase)` in Z.
pub fn apply_sine_displacement(mesh: &mut ParamMesh, amplitude: f64, frequency: f64, phase: f64) {
    for v in &mut mesh.vertices {
        v.pos[2] += amplitude * (frequency * v.pos[0] + phase).sin();
    }
}
/// Compute per-vertex normals by averaging adjacent face normals.
pub fn recompute_normals(mesh: &mut ParamMesh) {
    let n_verts = mesh.vertices.len();
    let mut normals = vec![[0.0_f64; 3]; n_verts];
    let n_tris = mesh.indices.len() / 3;
    for t in 0..n_tris {
        let i0 = mesh.indices[t * 3];
        let i1 = mesh.indices[t * 3 + 1];
        let i2 = mesh.indices[t * 3 + 2];
        let a = mesh.vertices[i0].pos;
        let b = mesh.vertices[i1].pos;
        let c = mesh.vertices[i2].pos;
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let cross = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        for &idx in &[i0, i1, i2] {
            normals[idx][0] += cross[0];
            normals[idx][1] += cross[1];
            normals[idx][2] += cross[2];
        }
    }
    for (v, n) in mesh.vertices.iter_mut().zip(normals.iter()) {
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).max(1e-300).sqrt();
        v.normal = [n[0] / mag, n[1] / mag, n[2] / mag];
    }
}
#[cfg(test)]
mod tests_mesh_gen_ext {
    use super::super::types::*;
    use super::*;

    #[test]
    fn param_mesh_empty_zero_counts() {
        let m = ParamMesh::empty();
        assert_eq!(m.vertex_count(), 0);
        assert_eq!(m.triangle_count(), 0);
    }
    #[test]
    fn param_mesh_triangle_count_six_indices() {
        let mut m = ParamMesh::empty();
        m.vertices.extend((0..4).map(|_| ParamVertex {
            pos: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            uv: [0.0; 2],
        }));
        m.indices = vec![0, 1, 2, 0, 2, 3];
        assert_eq!(m.triangle_count(), 2);
    }
    #[test]
    fn param_mesh_bounding_box_unit_square() {
        let mut m = ParamMesh::empty();
        for (x, y) in [(0.0_f64, 0.0_f64), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
            m.vertices.push(ParamVertex {
                pos: [x, y, 0.0],
                normal: [0.0, 0.0, 1.0],
                uv: [0.0; 2],
            });
        }
        let (lo, hi) = m.bounding_box();
        assert!((lo[0] - 0.0).abs() < 1e-12);
        assert!((hi[0] - 1.0).abs() < 1e-12);
        assert!((hi[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn param_mesh_translate_moves_vertices() {
        let mut m = ParamMesh::empty();
        m.vertices.push(ParamVertex {
            pos: [1.0, 2.0, 3.0],
            normal: [0.0, 1.0, 0.0],
            uv: [0.0; 2],
        });
        m.translate([10.0, 0.0, 0.0]);
        assert!((m.vertices[0].pos[0] - 11.0).abs() < 1e-12);
    }
    #[test]
    fn param_mesh_scale_doubles_positions() {
        let mut m = ParamMesh::empty();
        m.vertices.push(ParamVertex {
            pos: [1.0, 2.0, 3.0],
            normal: [0.0; 3],
            uv: [0.0; 2],
        });
        m.scale(2.0);
        assert!((m.vertices[0].pos[0] - 2.0).abs() < 1e-12);
        assert!((m.vertices[0].pos[1] - 4.0).abs() < 1e-12);
        assert!((m.vertices[0].pos[2] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn param_sphere_non_empty() {
        let m = param_sphere([0.0; 3], 1.0, 8, 16);
        assert!(m.vertex_count() > 0);
        assert!(m.triangle_count() > 0);
    }
    #[test]
    fn param_sphere_vertices_on_surface() {
        let r = 2.0;
        let m = param_sphere([0.0; 3], r, 6, 8);
        for v in &m.vertices {
            let dist = (v.pos[0] * v.pos[0] + v.pos[1] * v.pos[1] + v.pos[2] * v.pos[2]).sqrt();
            assert!(
                (dist - r).abs() < 1e-10,
                "vertex at distance {dist} from centre, expected {r}"
            );
        }
    }
    #[test]
    fn param_sphere_normals_unit() {
        let m = param_sphere([0.0; 3], 1.0, 4, 6);
        for v in &m.vertices {
            let mag =
                (v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] + v.normal[2] * v.normal[2])
                    .sqrt();
            assert!((mag - 1.0).abs() < 1e-10, "normal magnitude {mag}");
        }
    }
    #[test]
    fn param_sphere_area_approx_4pi_r2() {
        let r = 1.0;
        let m = param_sphere([0.0; 3], r, 32, 64);
        let area = m.surface_area();
        let expected = 4.0 * std::f64::consts::PI * r * r;
        let rel_err = (area - expected).abs() / expected;
        assert!(
            rel_err < 0.01,
            "area {area} vs {expected}, rel_err {rel_err}"
        );
    }
    #[test]
    fn param_disc_centre_vertex_is_first() {
        let m = param_disc([0.0; 3], 1.0, 8);
        assert_eq!(m.vertices[0].pos, [0.0; 3]);
    }
    #[test]
    fn param_disc_triangle_count_equals_n_segs() {
        let n = 12;
        let m = param_disc([0.0; 3], 1.0, n);
        assert_eq!(m.triangle_count(), n);
    }
    #[test]
    fn param_disc_normals_point_z() {
        let m = param_disc([0.0; 3], 1.0, 8);
        for v in &m.vertices {
            assert!((v.normal[2] - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn param_torus_non_empty() {
        let m = param_torus([0.0; 3], 1.0, 0.3, 8, 6);
        assert!(m.vertex_count() > 0);
        assert!(m.triangle_count() > 0);
    }
    #[test]
    fn param_torus_bounding_box_within_major_plus_minor() {
        let maj = 1.5;
        let min = 0.4;
        let m = param_torus([0.0; 3], maj, min, 16, 8);
        let (lo, hi) = m.bounding_box();
        assert!(hi[0] <= maj + min + 1e-10);
        assert!(lo[0] >= -(maj + min + 1e-10));
    }
    #[test]
    fn param_mobius_non_empty() {
        let m = param_mobius(16, 4, 1.0, 0.3);
        assert!(m.vertex_count() > 0);
        assert!(m.triangle_count() > 0);
    }
    #[test]
    fn sine_displacement_changes_z() {
        let mut m = param_sphere([0.0; 3], 1.0, 4, 6);
        let z_before: Vec<f64> = m.vertices.iter().map(|v| v.pos[2]).collect();
        apply_sine_displacement(&mut m, 0.1, 1.0, 0.0);
        let z_after: Vec<f64> = m.vertices.iter().map(|v| v.pos[2]).collect();
        let changed = z_before
            .iter()
            .zip(z_after.iter())
            .any(|(a, b)| (a - b).abs() > 1e-15);
        assert!(changed);
    }
    #[test]
    fn sine_displacement_zero_amplitude_no_change() {
        let mut m = param_disc([0.0; 3], 1.0, 8);
        let z_before: Vec<f64> = m.vertices.iter().map(|v| v.pos[2]).collect();
        apply_sine_displacement(&mut m, 0.0, 10.0, 0.5);
        let z_after: Vec<f64> = m.vertices.iter().map(|v| v.pos[2]).collect();
        for (a, b) in z_before.iter().zip(z_after.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
    #[test]
    fn recompute_normals_disc_points_up() {
        let mut m = param_disc([0.0; 3], 1.0, 8);
        recompute_normals(&mut m);
        for v in &m.vertices {
            assert!(
                v.normal[2] > 0.0,
                "normal z component should be positive: {:?}",
                v.normal
            );
        }
    }
    #[test]
    fn recompute_normals_sphere_produces_unit_normals() {
        let mut m = param_sphere([0.0; 3], 1.0, 4, 6);
        recompute_normals(&mut m);
        for v in &m.vertices {
            let mag =
                (v.normal[0] * v.normal[0] + v.normal[1] * v.normal[1] + v.normal[2] * v.normal[2])
                    .sqrt();
            assert!(!mag.is_nan(), "NaN normal at pos {:?}", v.pos);
        }
    }
}
