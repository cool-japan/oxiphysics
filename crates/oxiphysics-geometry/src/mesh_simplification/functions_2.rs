//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::mesh_simplification::HalfEdgeMesh;
    use crate::mesh_simplification::ProgressiveMesh;
    use crate::mesh_simplification::QemConfig;
    use crate::mesh_simplification::SymMat4;
    use crate::mesh_simplification::TextureSeamFlags;
    use crate::mesh_simplification::*;
    /// Create a simple unit cube mesh (8 vertices, 12 triangles).
    fn make_cube_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let tris = vec![
            [0, 1, 2],
            [0, 2, 3],
            [4, 6, 5],
            [4, 7, 6],
            [0, 5, 1],
            [0, 4, 5],
            [3, 2, 6],
            [3, 6, 7],
            [0, 3, 7],
            [0, 7, 4],
            [1, 5, 6],
            [1, 6, 2],
        ];
        (verts, tris)
    }
    /// Create a simple icosahedron-like mesh for testing.
    fn make_icosahedron() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let verts = vec![
            [-1.0, phi, 0.0],
            [1.0, phi, 0.0],
            [-1.0, -phi, 0.0],
            [1.0, -phi, 0.0],
            [0.0, -1.0, phi],
            [0.0, 1.0, phi],
            [0.0, -1.0, -phi],
            [0.0, 1.0, -phi],
            [phi, 0.0, -1.0],
            [phi, 0.0, 1.0],
            [-phi, 0.0, -1.0],
            [-phi, 0.0, 1.0],
        ];
        let tris = vec![
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
        (verts, tris)
    }
    /// Create a subdivided plane for testing.
    fn make_subdivided_plane(n: usize) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let mut verts = Vec::new();
        let mut tris = Vec::new();
        for iy in 0..=n {
            for ix in 0..=n {
                verts.push([ix as f64, iy as f64, 0.0]);
            }
        }
        let w = n + 1;
        for iy in 0..n {
            for ix in 0..n {
                let v00 = iy * w + ix;
                let v10 = iy * w + ix + 1;
                let v01 = (iy + 1) * w + ix;
                let v11 = (iy + 1) * w + ix + 1;
                tris.push([v00, v10, v11]);
                tris.push([v00, v11, v01]);
            }
        }
        (verts, tris)
    }
    #[test]
    fn test_symmat4_zero() {
        let m = SymMat4::zero();
        for &v in &m.m {
            assert!((v).abs() < 1e-15);
        }
    }
    #[test]
    fn test_symmat4_from_plane() {
        let m = SymMat4::from_plane(0.0, 0.0, 1.0, -1.0);
        assert!((m.m[7] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_symmat4_evaluate() {
        let m = SymMat4::from_plane(0.0, 0.0, 1.0, -2.0);
        let err = m.evaluate([5.0, 3.0, 2.0]);
        assert!(
            err.abs() < 1e-10,
            "Point on plane should have zero error, got {:.6e}",
            err
        );
    }
    #[test]
    fn test_symmat4_evaluate_off_plane() {
        let m = SymMat4::from_plane(0.0, 0.0, 1.0, 0.0);
        let err = m.evaluate([0.0, 0.0, 3.0]);
        assert!(
            (err - 9.0).abs() < 1e-10,
            "Error should be 9, got {:.6}",
            err
        );
    }
    #[test]
    fn test_symmat4_add() {
        let a = SymMat4::from_plane(1.0, 0.0, 0.0, 0.0);
        let b = SymMat4::from_plane(0.0, 1.0, 0.0, 0.0);
        let c = a.add(&b);
        assert!((c.m[0] - 1.0).abs() < 1e-10);
        assert!((c.m[4] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_symmat4_optimal_point() {
        let q1 = SymMat4::from_plane(1.0, 0.0, 0.0, -1.0);
        let q2 = SymMat4::from_plane(0.0, 1.0, 0.0, -2.0);
        let q3 = SymMat4::from_plane(0.0, 0.0, 1.0, -3.0);
        let q = q1.add(&q2).add(&q3);
        let opt = q.optimal_point();
        assert!(opt.is_some());
        let p = opt.unwrap();
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_solve_3x3_identity() {
        let a = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let b = [3.0, 5.0, 7.0];
        let x = solve_3x3(a, b).unwrap();
        assert!((x[0] - 3.0).abs() < 1e-10);
        assert!((x[1] - 5.0).abs() < 1e-10);
        assert!((x[2] - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_solve_3x3_singular() {
        let a = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let b = [1.0, 1.0, 1.0];
        assert!(solve_3x3(a, b).is_none());
    }
    #[test]
    fn test_half_edge_mesh_cube() {
        let (verts, tris) = make_cube_mesh();
        let mesh = HalfEdgeMesh::from_triangles(verts, &tris);
        assert_eq!(mesh.active_vertex_count(), 8);
        assert_eq!(mesh.active_face_count(), 12);
    }
    #[test]
    fn test_half_edge_boundary() {
        let (verts, tris) = make_subdivided_plane(2);
        let mesh = HalfEdgeMesh::from_triangles(verts, &tris);
        assert!(mesh.is_boundary_vertex(0));
    }
    #[test]
    fn test_half_edge_vertex_neighbors() {
        let (verts, tris) = make_cube_mesh();
        let mesh = HalfEdgeMesh::from_triangles(verts, &tris);
        let neighbors = mesh.vertex_neighbors(0);
        assert!(!neighbors.is_empty(), "Vertex 0 should have neighbors");
    }
    #[test]
    fn test_half_edge_face_normal() {
        let (verts, tris) = make_cube_mesh();
        let mesh = HalfEdgeMesh::from_triangles(verts, &tris);
        let n = mesh.face_normal(0);
        let l = len3(n);
        assert!(
            (l - 1.0).abs() < 1e-10,
            "Normal should be unit length, got {:.6}",
            l
        );
    }
    #[test]
    fn test_half_edge_face_area() {
        let (verts, tris) = make_cube_mesh();
        let mesh = HalfEdgeMesh::from_triangles(verts, &tris);
        let area = mesh.face_area(0);
        assert!(
            (area - 0.5).abs() < 1e-10,
            "Face area should be 0.5, got {:.6}",
            area
        );
    }
    #[test]
    fn test_triangle_aspect_ratio_equilateral() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, (3.0_f64).sqrt() / 2.0, 0.0];
        let ar = triangle_aspect_ratio(a, b, c);
        assert!(
            (ar - 2.0).abs() < 1e-6,
            "Equilateral AR should be 2, got {:.6}",
            ar
        );
    }
    #[test]
    fn test_triangle_min_angle_equilateral() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, (3.0_f64).sqrt() / 2.0, 0.0];
        let angle = triangle_min_angle(a, b, c);
        let expected = std::f64::consts::FRAC_PI_3;
        assert!(
            (angle - expected).abs() < 1e-6,
            "Min angle should be π/3, got {:.6}",
            angle
        );
    }
    #[test]
    fn test_triangle_quality_equilateral() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.5, (3.0_f64).sqrt() / 2.0, 0.0];
        let q = triangle_quality(a, b, c);
        assert!(
            (q - 1.0).abs() < 1e-6,
            "Equilateral quality should be 1.0, got {:.6}",
            q
        );
    }
    #[test]
    fn test_triangle_quality_degenerate() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [2.0, 0.0, 0.0];
        let q = triangle_quality(a, b, c);
        assert!(
            q < 0.01,
            "Degenerate quality should be near 0, got {:.6}",
            q
        );
    }
    #[test]
    fn test_qem_simplify_cube() {
        let (verts, tris) = make_cube_mesh();
        let config = QemConfig {
            target_triangles: 6,
            max_error: 0.0,
            preserve_boundary: false,
            max_normal_deviation: 0.0,
            preserve_texture_seams: false,
            boundary_weight: 0.0,
        };
        let (sv, st) = qem_simplify(&verts, &tris, &config);
        assert!(st.len() <= 12, "Should have fewer triangles");
        assert!(!sv.is_empty());
    }
    #[test]
    fn test_qem_simplify_icosahedron() {
        let (verts, tris) = make_icosahedron();
        let config = QemConfig {
            target_triangles: 10,
            max_error: 0.0,
            preserve_boundary: false,
            max_normal_deviation: 0.0,
            preserve_texture_seams: false,
            boundary_weight: 0.0,
        };
        let (sv, st) = qem_simplify(&verts, &tris, &config);
        assert!(st.len() <= 20);
        assert!(!sv.is_empty());
    }
    #[test]
    fn test_vertex_clustering() {
        let (verts, tris) = make_subdivided_plane(5);
        let (cv, ct) = vertex_clustering(&verts, &tris, 2.5);
        assert!(cv.len() < verts.len(), "Should reduce vertex count");
        assert!(ct.len() <= tris.len(), "Should not increase triangle count");
    }
    #[test]
    fn test_edge_length_decimation() {
        let (verts, tris) = make_subdivided_plane(4);
        let (dv, dt) = edge_length_decimation(&verts, &tris, 0.5);
        assert_eq!(dv.len(), verts.len());
        assert_eq!(dt.len(), tris.len());
    }
    #[test]
    fn test_edge_length_decimation_with_short_edges() {
        let mut verts = vec![
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
        ];
        let tris = vec![[0, 1, 3], [1, 2, 3]];
        let (dv, _dt) = edge_length_decimation(&verts, &tris, 0.5);
        assert!(dv.len() < verts.len(), "Should merge close vertices");
        let _ = &mut verts;
    }
    #[test]
    fn test_progressive_mesh_build() {
        let (verts, tris) = make_icosahedron();
        let pm = ProgressiveMesh::build(&verts, &tris, 10);
        let (bv, bt) = pm.base_mesh();
        assert!(!bv.is_empty());
        assert!(bt.len() <= 20);
    }
    #[test]
    fn test_point_to_triangle_dist_on_vertex() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let d = point_to_triangle_dist(a, a, b, c);
        assert!(d < 1e-10, "Distance to vertex should be 0, got {:.6e}", d);
    }
    #[test]
    fn test_point_to_triangle_dist_above() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let p = [0.25, 0.25, 1.0];
        let d = point_to_triangle_dist(p, a, b, c);
        assert!(
            (d - 1.0).abs() < 1e-10,
            "Distance should be 1.0, got {:.6}",
            d
        );
    }
    #[test]
    fn test_hausdorff_same_mesh() {
        let (verts, tris) = make_cube_mesh();
        let d = hausdorff_symmetric(&verts, &tris, &verts, &tris, 100);
        assert!(
            d < 1e-10,
            "Hausdorff of identical mesh should be 0, got {:.6e}",
            d
        );
    }
    #[test]
    fn test_validate_topology_cube() {
        let (verts, tris) = make_cube_mesh();
        let val = validate_topology(&verts, &tris);
        assert!(val.is_valid, "Cube should have valid topology");
        assert_eq!(val.non_manifold_edges, 0);
        assert_eq!(val.degenerate_triangles, 0);
    }
    #[test]
    fn test_validate_topology_degenerate() {
        let verts = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let tris = vec![[0, 0, 1]];
        let val = validate_topology(&verts, &tris);
        assert!(!val.is_valid);
        assert_eq!(val.degenerate_triangles, 1);
    }
    #[test]
    fn test_texture_seam_flags() {
        let mut flags = TextureSeamFlags::new();
        flags.mark_seam(0, 1);
        flags.mark_seam(3, 2);
        assert!(flags.is_seam(0, 1));
        assert!(flags.is_seam(1, 0));
        assert!(flags.is_seam(2, 3));
        assert!(!flags.is_seam(0, 2));
        assert_eq!(flags.count(), 2);
    }
    #[test]
    fn test_mesh_stats() {
        let (verts, tris) = make_cube_mesh();
        let stats = compute_mesh_stats(&verts, &tris);
        assert_eq!(stats.n_vertices, 8);
        assert_eq!(stats.n_triangles, 12);
        assert!(stats.avg_edge_length > 0.0);
        assert!(stats.min_quality > 0.0);
    }
    #[test]
    fn test_boundary_edges_cube() {
        let (_, tris) = make_cube_mesh();
        let boundary = find_boundary_edges(&tris);
        assert!(
            boundary.is_empty(),
            "Closed cube should have no boundary edges"
        );
    }
    #[test]
    fn test_boundary_edges_plane() {
        let (_, tris) = make_subdivided_plane(2);
        let boundary = find_boundary_edges(&tris);
        assert!(
            !boundary.is_empty(),
            "Open plane should have boundary edges"
        );
    }
    #[test]
    fn test_euler_characteristic_cube() {
        let (verts, tris) = make_cube_mesh();
        let chi = euler_characteristic(verts.len(), &tris);
        assert_eq!(chi, 2, "Cube Euler characteristic should be 2, got {}", chi);
    }
    #[test]
    fn test_simplify_to_ratio() {
        let (verts, tris) = make_subdivided_plane(4);
        let original_count = tris.len();
        let (_, st) = simplify_to_ratio(&verts, &tris, 0.5);
        assert!(
            st.len() <= original_count,
            "Should reduce triangles from {} to ≤ {}",
            original_count,
            original_count / 2
        );
    }
    #[test]
    fn test_mesh_surface_area_cube() {
        let (verts, tris) = make_cube_mesh();
        let area = mesh_surface_area(&verts, &tris);
        assert!(
            (area - 6.0).abs() < 1e-10,
            "Cube surface area should be 6, got {:.6}",
            area
        );
    }
    #[test]
    fn test_mesh_signed_volume_cube() {
        let (verts, tris) = make_cube_mesh();
        let vol = mesh_signed_volume(&verts, &tris);
        assert!(
            (vol.abs() - 1.0).abs() < 1e-10,
            "Cube volume should be 1.0, got {:.6}",
            vol.abs()
        );
    }
    #[test]
    fn test_connected_components_single() {
        let (verts, tris) = make_cube_mesh();
        let cc = count_connected_components(verts.len(), &tris);
        assert_eq!(cc, 1, "Cube should be one component");
    }
    #[test]
    fn test_connected_components_two() {
        let verts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [5.0, 5.0, 5.0],
            [6.0, 5.0, 5.0],
            [5.0, 6.0, 5.0],
        ];
        let tris = vec![[0, 1, 2], [3, 4, 5]];
        let cc = count_connected_components(verts.len(), &tris);
        assert_eq!(cc, 2, "Should have 2 components");
    }
    #[test]
    fn test_compact_mesh() {
        let verts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        let tris = vec![[0, 2, 3]];
        let (cv, ct) = compact_mesh(&verts, &tris);
        assert_eq!(cv.len(), 3, "Should have 3 vertices after compaction");
        assert_eq!(ct.len(), 1);
    }
    #[test]
    fn test_extract_triangles_from_halfedge() {
        let (verts, tris) = make_cube_mesh();
        let mesh = HalfEdgeMesh::from_triangles(verts, &tris);
        let extracted = mesh.extract_triangles();
        assert_eq!(extracted.len(), tris.len());
    }
    #[test]
    fn test_qem_config_default() {
        let config = QemConfig::default_config(100);
        assert_eq!(config.target_triangles, 50);
        assert!(config.preserve_boundary);
    }
    #[test]
    fn test_qem_config_error_threshold() {
        let config = QemConfig::with_error_threshold(0.01);
        assert_eq!(config.target_triangles, 0);
        assert!((config.max_error - 0.01).abs() < 1e-15);
    }
}
