//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::functions::{compute_normal, normalize3_f32, triangle_area};
use super::types::{StlMesh, StlQualityMetrics, TriangleMesh, WeldedMesh};

#[cfg(test)]
mod tests_ext {
    use super::*;
    use crate::stl::types::*;
    #[test]
    fn test_validate_single_triangle_not_watertight() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let report = validate_mesh(&m);
        assert!(!report.is_watertight, "single triangle has open edges");
        assert!(report.boundary_edge_count > 0);
    }
    #[test]
    fn test_validate_tetrahedron_watertight() {
        let p0 = [0.0f32, 0.0, 0.0];
        let p1 = [1.0, 0.0, 0.0];
        let p2 = [0.5, 1.0, 0.0];
        let p3 = [0.5, 0.33, 1.0];
        let mut m = StlMesh::new("tet");
        m.add_triangle(p0, p1, p2);
        m.add_triangle(p0, p1, p3);
        m.add_triangle(p0, p2, p3);
        m.add_triangle(p1, p2, p3);
        let report = validate_mesh(&m);
        assert!(report.is_watertight, "tetrahedron should be watertight");
    }
    #[test]
    fn test_validate_degenerate_triangle() {
        let mut m = StlMesh::new("d");
        m.triangles.push(StlTriangle {
            normal: [0.0, 0.0, 1.0],
            v0: [1.0, 1.0, 1.0],
            v1: [1.0, 1.0, 1.0],
            v2: [1.0, 1.0, 1.0],
        });
        let report = validate_mesh(&m);
        assert_eq!(report.degenerate_count, 1);
    }
    #[test]
    fn test_stl_color_encode_decode_roundtrip() {
        let c = StlColor::new(10, 20, 15);
        let attr = c.encode();
        let decoded = StlColor::decode(attr).expect("should decode");
        assert_eq!(decoded.r, 10);
        assert_eq!(decoded.g, 20);
        assert_eq!(decoded.b, 15);
    }
    #[test]
    fn test_stl_color_no_valid_bit() {
        assert!(StlColor::decode(0x0FFF).is_none());
    }
    #[test]
    fn test_stl_color_to_normalized() {
        let c = StlColor::new(31, 0, 0);
        let norm = c.to_normalized();
        assert!((norm[0] - 1.0).abs() < 1e-4, "R=31 → 1.0");
        assert!(norm[1].abs() < 1e-4);
    }
    #[test]
    fn test_binary_with_color_roundtrip() {
        let mut m = StlMesh::new("colored");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let colors = vec![StlColor::new(31, 10, 5)];
        let bytes = to_binary_bytes_with_color(&m, &colors);
        let attr = u16::from_le_bytes([bytes[132], bytes[133]]);
        let decoded = StlColor::decode(attr).expect("color should decode");
        assert_eq!(decoded.r, 31);
        assert_eq!(decoded.g, 10);
        assert_eq!(decoded.b, 5);
    }
    #[test]
    fn test_merge_meshes_count() {
        let mut a = StlMesh::new("a");
        a.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let mut b = StlMesh::new("b");
        b.add_triangle([2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]);
        b.add_triangle([3.0, 0.0, 0.0], [4.0, 0.0, 0.0], [3.0, 1.0, 0.0]);
        let merged = merge_meshes_named(&a, &b, "merged");
        assert_eq!(merged.triangles.len(), 3);
    }
    #[test]
    fn test_fix_normals() {
        let mut m = StlMesh::new("t");
        m.triangles.push(StlTriangle {
            normal: [1.0, 0.0, 0.0],
            v0: [0.0, 0.0, 0.0],
            v1: [1.0, 0.0, 0.0],
            v2: [0.0, 1.0, 0.0],
        });
        fix_normals(&mut m);
        let n = m.triangles[0].normal;
        assert!(
            n[2] > 0.9,
            "z-component of normal should be ~1, got {}",
            n[2]
        );
    }
    #[test]
    fn test_flip_winding_reverses_normal() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let n_before = m.triangles[0].normal[2];
        flip_winding(&mut m);
        let n_after = m.triangles[0].normal[2];
        assert!(n_before > 0.0, "before flip z>0");
        assert!(n_after < 0.0, "after flip z<0");
    }
    #[test]
    fn test_triangle_bounding_box_empty() {
        assert!(triangle_bounding_box(&[]).is_none());
    }
    #[test]
    fn test_triangle_bounding_box_values() {
        let tris = vec![StlTriangle {
            normal: [0.0, 0.0, 1.0],
            v0: [1.0, 2.0, 3.0],
            v1: [4.0, 5.0, 6.0],
            v2: [7.0, 8.0, 9.0],
        }];
        let (lo, hi) = triangle_bounding_box(&tris).unwrap();
        assert!((lo[0] - 1.0).abs() < 1e-5);
        assert!((hi[0] - 7.0).abs() < 1e-5);
    }
    #[test]
    fn test_stl_to_obj_roundtrip() {
        let mut m = StlMesh::new("cube_face");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let obj = stl_to_obj(&m);
        assert!(obj.contains("v 0 0 0"), "should contain vertex");
        assert!(obj.contains("f 1 2 3"), "should contain face");
    }
    #[test]
    fn test_obj_to_stl_roundtrip() {
        let mut m = StlMesh::new("round");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let obj = stl_to_obj(&m);
        let m2 = obj_to_stl(&obj);
        assert_eq!(m2.triangles.len(), 1);
        assert!((m2.triangles[0].v0[0]).abs() < 0.01);
    }
    #[test]
    fn test_repair_close_holes_adds_caps() {
        let mut m = StlMesh::new("open");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let before = m.triangles.len();
        repair_close_holes(&mut m);
        assert!(
            m.triangles.len() > before,
            "repair should add cap triangles"
        );
    }
}
/// Weld vertices of an STL mesh within `tolerance`.
///
/// Returns a [`WeldedMesh`] where vertices closer than `tolerance` are merged
/// into a single canonical vertex.
pub fn weld_vertices(mesh: &StlMesh, tolerance: f32) -> WeldedMesh {
    let mut unique: Vec<[f32; 3]> = Vec::new();
    let mut remap: Vec<usize> = Vec::new();
    let all_verts: Vec<[f32; 3]> = mesh
        .triangles
        .iter()
        .flat_map(|t| [t.v0, t.v1, t.v2])
        .collect();
    for v in &all_verts {
        let idx = unique.iter().position(|u| {
            let dx = u[0] - v[0];
            let dy = u[1] - v[1];
            let dz = u[2] - v[2];
            (dx * dx + dy * dy + dz * dz).sqrt() <= tolerance
        });
        if let Some(i) = idx {
            remap.push(i);
        } else {
            remap.push(unique.len());
            unique.push(*v);
        }
    }
    let triangles: Vec<[usize; 3]> = (0..mesh.triangles.len())
        .map(|i| [remap[i * 3], remap[i * 3 + 1], remap[i * 3 + 2]])
        .collect();
    WeldedMesh {
        vertices: unique,
        triangles,
    }
}
/// Count the number of unique vertices in an STL mesh within `tolerance`.
pub fn count_unique_vertices(mesh: &StlMesh, tolerance: f32) -> usize {
    weld_vertices(mesh, tolerance).vertices.len()
}
/// Compute quality metrics for an STL mesh.
pub fn compute_quality_metrics(mesh: &StlMesh) -> StlQualityMetrics {
    let n = mesh.triangles.len();
    if n == 0 {
        return StlQualityMetrics {
            max_aspect_ratio: 0.0,
            min_area: 0.0,
            max_area: 0.0,
            degenerate_fraction: 0.0,
            has_bad_geometry: false,
        };
    }
    let mut max_ar = 0.0_f32;
    let mut min_area = f32::MAX;
    let mut max_area = 0.0_f32;
    let mut degen = 0usize;
    let mut has_bad = false;
    for tri in &mesh.triangles {
        let verts = [tri.v0, tri.v1, tri.v2];
        if verts
            .iter()
            .any(|v| v.iter().any(|f| f.is_nan() || f.is_infinite()))
        {
            has_bad = true;
            continue;
        }
        let area = triangle_area(tri);
        if area < 1e-10 {
            degen += 1;
        }
        if area < min_area {
            min_area = area;
        }
        if area > max_area {
            max_area = area;
        }
        let edges = [
            edge_len(tri.v0, tri.v1),
            edge_len(tri.v1, tri.v2),
            edge_len(tri.v2, tri.v0),
        ];
        let max_edge = edges.iter().cloned().fold(f32::MIN, f32::max);
        if area > 1e-15 {
            let ar = max_edge * max_edge / (2.0 * area);
            if ar > max_ar {
                max_ar = ar;
            }
        }
    }
    StlQualityMetrics {
        max_aspect_ratio: max_ar,
        min_area,
        max_area,
        degenerate_fraction: degen as f32 / n as f32,
        has_bad_geometry: has_bad,
    }
}
pub(super) fn edge_len(a: [f32; 3], b: [f32; 3]) -> f32 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let dz = b[2] - a[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
/// Apply `iterations` of Laplacian smoothing to a welded mesh.
///
/// Each vertex is moved toward the average of its neighbours.
/// `factor` ∈ (0, 1] controls the strength (0.5 is typical).
pub fn laplacian_smooth(mesh: &mut WeldedMesh, iterations: usize, factor: f32) {
    let nv = mesh.vertices.len();
    if nv == 0 {
        return;
    }
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for tri in &mesh.triangles {
        let pairs = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])];
        for (a, b) in pairs {
            if !adj[a].contains(&b) {
                adj[a].push(b);
            }
            if !adj[b].contains(&a) {
                adj[b].push(a);
            }
        }
    }
    for _ in 0..iterations {
        let old = mesh.vertices.clone();
        for i in 0..nv {
            if adj[i].is_empty() {
                continue;
            }
            let n = adj[i].len() as f32;
            let avg = adj[i].iter().fold([0.0_f32; 3], |acc, &j| {
                [acc[0] + old[j][0], acc[1] + old[j][1], acc[2] + old[j][2]]
            });
            let avg = [avg[0] / n, avg[1] / n, avg[2] / n];
            mesh.vertices[i][0] = old[i][0] + factor * (avg[0] - old[i][0]);
            mesh.vertices[i][1] = old[i][1] + factor * (avg[1] - old[i][1]);
            mesh.vertices[i][2] = old[i][2] + factor * (avg[2] - old[i][2]);
        }
    }
}
/// Subdivide each triangle by inserting midpoints on each edge.
///
/// Each input triangle produces 4 output triangles.  This is the simplest
/// possible 1-to-4 refinement scheme.
pub fn subdivide_midpoint(mesh: &StlMesh) -> StlMesh {
    let mut out = StlMesh::new(&format!("{}_subdivided", mesh.name));
    for tri in &mesh.triangles {
        let m01 = midpoint(tri.v0, tri.v1);
        let m12 = midpoint(tri.v1, tri.v2);
        let m20 = midpoint(tri.v2, tri.v0);
        out.add_triangle(tri.v0, m01, m20);
        out.add_triangle(m01, tri.v1, m12);
        out.add_triangle(m20, m12, tri.v2);
        out.add_triangle(m01, m12, m20);
    }
    out
}
pub(super) fn midpoint(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        (a[0] + b[0]) / 2.0,
        (a[1] + b[1]) / 2.0,
        (a[2] + b[2]) / 2.0,
    ]
}
/// Convert an STL mesh into a proper indexed `TriangleMesh` with smooth normals.
///
/// Vertices are welded within `tolerance`; normals are averaged over all
/// triangles that share a vertex.
pub fn stl_to_triangle_mesh(mesh: &StlMesh, tolerance: f32) -> TriangleMesh {
    let welded = weld_vertices(mesh, tolerance);
    let nv = welded.vertices.len();
    let mut normal_sum: Vec<[f32; 3]> = vec![[0.0; 3]; nv];
    let mut normal_count: Vec<usize> = vec![0; nv];
    for tri in &welded.triangles {
        let v0 = welded.vertices[tri[0]];
        let v1 = welded.vertices[tri[1]];
        let v2 = welded.vertices[tri[2]];
        let n = compute_normal(v0, v1, v2);
        for &vi in tri {
            normal_sum[vi][0] += n[0];
            normal_sum[vi][1] += n[1];
            normal_sum[vi][2] += n[2];
            normal_count[vi] += 1;
        }
    }
    let normals: Vec<[f32; 3]> = (0..nv)
        .map(|i| {
            let c = normal_count[i].max(1) as f32;
            normalize3_f32([
                normal_sum[i][0] / c,
                normal_sum[i][1] / c,
                normal_sum[i][2] / c,
            ])
        })
        .collect();
    TriangleMesh {
        positions: welded.vertices,
        normals,
        indices: welded.triangles,
    }
}
/// Merge two STL meshes and weld shared vertices within `tolerance`.
pub fn merge_and_weld(a: &StlMesh, b: &StlMesh, tolerance: f32) -> WeldedMesh {
    let combined = merge_meshes(a, b);
    weld_vertices(&combined, tolerance)
}
#[cfg(test)]
mod extended_tests {
    use super::*;
    use crate::stl::StlColor;

    fn small_tet_mesh() -> StlMesh {
        let mut m = StlMesh::new("tet");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        m.add_triangle([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        m.add_triangle([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        m.add_triangle([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]);
        m
    }
    #[test]
    fn test_weld_vertices_reduces_count() {
        let m = small_tet_mesh();
        let w = weld_vertices(&m, 1e-4);
        assert!(
            w.vertices.len() <= 4,
            "welded tet should have ≤4 unique vertices, got {}",
            w.vertices.len()
        );
    }
    #[test]
    fn test_weld_vertices_preserves_triangle_count() {
        let m = small_tet_mesh();
        let w = weld_vertices(&m, 1e-4);
        assert_eq!(w.triangles.len(), 4);
    }
    #[test]
    fn test_weld_vertices_large_tolerance_collapses_all() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [0.001, 0.0, 0.0], [0.0, 0.001, 0.0]);
        let w = weld_vertices(&m, 0.01);
        assert!(
            w.vertices.len() == 1,
            "all verts should merge to 1, got {}",
            w.vertices.len()
        );
    }
    #[test]
    fn test_weld_vertices_no_tolerance_keeps_all() {
        let m = small_tet_mesh();
        let w = weld_vertices(&m, 0.0);
        assert!(w.vertices.len() <= 12);
    }
    #[test]
    fn test_count_unique_vertices() {
        let m = small_tet_mesh();
        let n = count_unique_vertices(&m, 1e-4);
        assert!(n <= 4, "got {n}");
    }
    #[test]
    fn test_quality_metrics_equilateral() {
        let mut m = StlMesh::new("equi");
        let h = (3.0_f32).sqrt() / 2.0;
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]);
        let q = compute_quality_metrics(&m);
        assert!(q.max_aspect_ratio > 0.0);
        assert!(q.min_area > 0.0);
        assert_eq!(q.degenerate_fraction, 0.0);
        assert!(!q.has_bad_geometry);
    }
    #[test]
    fn test_quality_metrics_degenerate() {
        let mut m = StlMesh::new("degen");
        m.add_triangle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let q = compute_quality_metrics(&m);
        assert!(q.degenerate_fraction > 0.9);
    }
    #[test]
    fn test_quality_metrics_empty() {
        let m = StlMesh::new("empty");
        let q = compute_quality_metrics(&m);
        assert_eq!(q.max_aspect_ratio, 0.0);
        assert!(!q.has_bad_geometry);
    }
    #[test]
    fn test_quality_metrics_two_triangles() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        m.add_triangle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let q = compute_quality_metrics(&m);
        assert_eq!(q.degenerate_fraction, 0.5);
    }
    #[test]
    fn test_laplacian_smooth_does_not_panic() {
        let m = small_tet_mesh();
        let mut w = weld_vertices(&m, 1e-4);
        laplacian_smooth(&mut w, 5, 0.5);
        for v in &w.vertices {
            assert!(v[0].is_finite() && v[1].is_finite() && v[2].is_finite());
        }
    }
    #[test]
    fn test_laplacian_smooth_zero_iterations_unchanged() {
        let m = small_tet_mesh();
        let mut w = weld_vertices(&m, 1e-4);
        let before = w.vertices.clone();
        laplacian_smooth(&mut w, 0, 0.5);
        assert_eq!(w.vertices, before);
    }
    #[test]
    fn test_laplacian_smooth_many_iterations_stays_finite() {
        let mut m = StlMesh::new("grid");
        for i in 0..2 {
            for j in 0..2 {
                let ix = i as f32;
                let jx = j as f32;
                m.add_triangle([ix, jx, 0.0], [ix + 1.0, jx, 0.0], [ix, jx + 1.0, 0.0]);
                m.add_triangle(
                    [ix + 1.0, jx, 0.0],
                    [ix + 1.0, jx + 1.0, 0.0],
                    [ix, jx + 1.0, 0.0],
                );
            }
        }
        let mut w = weld_vertices(&m, 1e-4);
        laplacian_smooth(&mut w, 50, 0.5);
        for v in &w.vertices {
            assert!(v[0].is_finite() && v[1].is_finite() && v[2].is_finite());
        }
    }
    #[test]
    fn test_subdivide_midpoint_count() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let sub = subdivide_midpoint(&m);
        assert_eq!(sub.triangles.len(), 4, "1 triangle → 4 after subdivision");
    }
    #[test]
    fn test_subdivide_midpoint_preserves_area() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        let area_before = m.surface_area();
        let sub = subdivide_midpoint(&m);
        let area_after = sub.surface_area();
        assert!(
            (area_before - area_after).abs() < 1e-4,
            "subdivision should preserve area: before={area_before}, after={area_after}"
        );
    }
    #[test]
    fn test_subdivide_twice_count() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let s1 = subdivide_midpoint(&m);
        let s2 = subdivide_midpoint(&s1);
        assert_eq!(s2.triangles.len(), 16, "2 rounds → 16 triangles");
    }
    #[test]
    fn test_subdivide_cube_total_area() {
        let mut m = StlMesh::new("cube");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]);
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        let area_before = m.surface_area();
        let sub = subdivide_midpoint(&m);
        let area_after = sub.surface_area();
        assert!((area_before - area_after).abs() < 1e-4);
    }
    #[test]
    fn test_stl_to_triangle_mesh_vertex_count() {
        let m = small_tet_mesh();
        let tm = stl_to_triangle_mesh(&m, 1e-4);
        assert!(
            tm.positions.len() <= 4,
            "welded tet should have ≤4 positions"
        );
        assert_eq!(tm.positions.len(), tm.normals.len());
    }
    #[test]
    fn test_stl_to_triangle_mesh_normals_unit_length() {
        let m = small_tet_mesh();
        let tm = stl_to_triangle_mesh(&m, 1e-4);
        for n in &tm.normals {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(
                (len - 1.0).abs() < 0.01 || len < 1e-10,
                "normal should be unit length, got {len}"
            );
        }
    }
    #[test]
    fn test_stl_to_triangle_mesh_index_validity() {
        let m = small_tet_mesh();
        let tm = stl_to_triangle_mesh(&m, 1e-4);
        let nv = tm.positions.len();
        for tri in &tm.indices {
            assert!(
                tri[0] < nv && tri[1] < nv && tri[2] < nv,
                "all triangle indices must be valid"
            );
        }
    }
    #[test]
    fn test_merge_and_weld_total_triangles() {
        let a = small_tet_mesh();
        let mut b = StlMesh::new("b");
        b.add_triangle([5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [5.0, 1.0, 0.0]);
        let w = merge_and_weld(&a, &b, 1e-4);
        assert_eq!(w.triangles.len(), 5);
    }
    #[test]
    fn test_merge_and_weld_deduplicates_shared_vertices() {
        let mut a = StlMesh::new("a");
        a.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let mut b = StlMesh::new("b");
        b.add_triangle([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let w = merge_and_weld(&a, &b, 1e-4);
        assert!(w.vertices.len() < 6);
    }
    #[test]
    fn test_edge_len_known_value() {
        let len = edge_len([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((len - 5.0).abs() < 1e-5, "3-4-5 triangle, got {len}");
    }
    #[test]
    fn test_edge_len_zero() {
        let len = edge_len([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]);
        assert!(len.abs() < 1e-6);
    }
    #[test]
    fn test_midpoint_known_value() {
        let m = midpoint([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        assert!((m[0] - 1.0).abs() < 1e-5);
        assert!((m[1] - 2.0).abs() < 1e-5);
        assert!((m[2] - 3.0).abs() < 1e-5);
    }
    #[test]
    fn test_validate_mesh_tet_watertight() {
        let m = small_tet_mesh();
        let r = validate_mesh(&m);
        assert!(r.is_watertight, "tetrahedron should be watertight");
        assert!(r.is_manifold, "tetrahedron should be manifold");
        assert_eq!(r.degenerate_count, 0);
    }
    #[test]
    fn test_validate_mesh_open_has_boundary() {
        let mut m = StlMesh::new("open");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]);
        let r = validate_mesh(&m);
        assert!(!r.is_watertight);
        assert!(r.boundary_edge_count > 0);
    }
    #[test]
    fn test_validate_mesh_degenerate_counted() {
        let mut m = StlMesh::new("degen");
        m.add_triangle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let r = validate_mesh(&m);
        assert_eq!(r.degenerate_count, 1);
    }
    #[test]
    fn test_stl_color_encode_decode_roundtrip() {
        let c = StlColor::new(15, 20, 5);
        let encoded = c.encode();
        let decoded = StlColor::decode(encoded).expect("should decode");
        assert_eq!(decoded.r, 15);
        assert_eq!(decoded.g, 20);
        assert_eq!(decoded.b, 5);
    }
    #[test]
    fn test_stl_color_decode_invalid() {
        let attr: u16 = 0x1234;
        assert!(StlColor::decode(attr).is_none());
    }
    #[test]
    fn test_stl_color_to_normalized() {
        let c = StlColor::new(31, 0, 0);
        let n = c.to_normalized();
        assert!((n[0] - 1.0).abs() < 1e-5);
        assert!(n[1].abs() < 1e-5);
    }
    #[test]
    fn test_stl_color_mask_overflow() {
        let c = StlColor::new(255, 255, 255);
        assert_eq!(c.r, 31);
        assert_eq!(c.g, 31);
        assert_eq!(c.b, 31);
    }
    #[test]
    fn test_binary_with_color_no_colors() {
        let mut m = StlMesh::new("t");
        m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let bytes = to_binary_bytes_with_color(&m, &[]);
        assert_eq!(bytes.len(), 84 + 50);
        let attr = u16::from_le_bytes([bytes[132], bytes[133]]);
        assert_eq!(attr, 0);
    }
    #[test]
    fn test_binary_with_color_correct_length() {
        let mut m = StlMesh::new("t");
        for _ in 0..5 {
            m.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        }
        let colors: Vec<StlColor> = (0..5).map(|i| StlColor::new(i as u8, 0, 0)).collect();
        let bytes = to_binary_bytes_with_color(&m, &colors);
        assert_eq!(bytes.len(), 84 + 5 * 50);
    }
    #[test]
    fn test_merge_meshes_named_empty_both() {
        let a = StlMesh::new("a");
        let b = StlMesh::new("b");
        let m = merge_meshes_named(&a, &b, "ab");
        assert_eq!(m.name, "ab");
        assert!(m.triangles.is_empty());
    }
    #[test]
    fn test_merge_meshes_named_preserves_order() {
        let mut a = StlMesh::new("a");
        a.add_triangle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let mut b = StlMesh::new("b");
        b.add_triangle([2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]);
        let m = merge_meshes_named(&a, &b, "result");
        assert_eq!(m.triangles[0].v0, [0.0, 0.0, 0.0]);
        assert_eq!(m.triangles[1].v0, [2.0, 0.0, 0.0]);
    }
}
