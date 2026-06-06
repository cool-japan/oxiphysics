//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {
    use crate::geodesic::GeoMesh;
    use crate::geodesic::*;
    fn strip_mesh(n_tris: usize) -> GeoMesh {
        let mut verts = Vec::new();
        for i in 0..=n_tris {
            verts.push([i as f64, 0.0, 0.0]);
            verts.push([i as f64, 1.0, 0.0]);
        }
        let mut faces = Vec::new();
        for i in 0..n_tris {
            let b0 = 2 * i;
            let t0 = 2 * i + 1;
            let b1 = 2 * i + 2;
            let t1 = 2 * i + 3;
            faces.push([b0, b1, t0]);
            faces.push([b1, t1, t0]);
        }
        GeoMesh::new(verts, faces)
    }
    fn single_triangle() -> GeoMesh {
        let h = (3.0_f64).sqrt() / 2.0;
        GeoMesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0]],
            vec![[0, 1, 2]],
        )
    }
    fn tet_surface_mesh() -> GeoMesh {
        let v = vec![
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let f = vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
        GeoMesh::new(v, f)
    }
    #[test]
    fn test_geomesh_n_vertices() {
        let m = strip_mesh(3);
        assert_eq!(m.n_vertices(), 8);
    }
    #[test]
    fn test_geomesh_n_faces() {
        let m = strip_mesh(4);
        assert_eq!(m.n_faces(), 8);
    }
    #[test]
    fn test_geomesh_face_area_equilateral() {
        let m = single_triangle();
        let area = m.face_area(0);
        let expected = (3.0_f64).sqrt() / 4.0;
        assert!(
            (area - expected).abs() < 1e-10,
            "area={area}, expected={expected}"
        );
    }
    #[test]
    fn test_geomesh_face_area_right_triangle() {
        let m = GeoMesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let area = m.face_area(0);
        assert!((area - 0.5).abs() < 1e-10, "area={area}");
    }
    #[test]
    fn test_build_adjacency_single_triangle() {
        let m = single_triangle();
        let adj = m.build_adjacency();
        assert_eq!(adj.len(), 3);
        for (v, nbrs) in adj.iter().enumerate() {
            assert_eq!(nbrs.len(), 2, "vertex {v} should have 2 neighbours");
        }
    }
    #[test]
    fn test_build_adjacency_strip_mesh() {
        let m = strip_mesh(2);
        let adj = m.build_adjacency();
        for (v, nbrs) in adj.iter().enumerate() {
            assert!(!nbrs.is_empty(), "vertex {v} should have neighbours");
        }
    }
    #[test]
    fn test_dijkstra_source_zero() {
        let m = strip_mesh(3);
        let d = dijkstra_geodesic(&m, 0);
        assert_eq!(d[0], 0.0, "distance to source should be 0");
    }
    #[test]
    fn test_dijkstra_all_finite_connected() {
        let m = tet_surface_mesh();
        let d = dijkstra_geodesic(&m, 0);
        for (i, &di) in d.iter().enumerate() {
            assert!(di.is_finite(), "vertex {i} should be reachable");
        }
    }
    #[test]
    fn test_dijkstra_symmetric() {
        let m = tet_surface_mesh();
        let d0 = dijkstra_geodesic(&m, 0);
        let d1 = dijkstra_geodesic(&m, 1);
        assert!(
            (d0[1] - d1[0]).abs() < 1e-10,
            "geodesic should be symmetric: {} vs {}",
            d0[1],
            d1[0]
        );
    }
    #[test]
    fn test_dijkstra_non_negative() {
        let m = strip_mesh(4);
        let d = dijkstra_geodesic(&m, 0);
        for &di in &d {
            assert!(
                di >= 0.0 || !di.is_finite(),
                "distances must be non-negative"
            );
        }
    }
    #[test]
    fn test_dijkstra_triangle_distances() {
        let m = single_triangle();
        let d = dijkstra_geodesic(&m, 0);
        assert_eq!(d[0], 0.0);
        assert!(
            (d[1] - 1.0).abs() < 1e-10,
            "d[1] should be edge length 1.0, got {}",
            d[1]
        );
    }
    #[test]
    fn test_dijkstra_triangle_inequality() {
        let m = tet_surface_mesh();
        let d = dijkstra_geodesic(&m, 0);
        let adj = m.build_adjacency();
        for (u, nbrs) in adj.iter().enumerate() {
            for &(v, w) in nbrs {
                assert!(
                    d[v] <= d[u] + w + 1e-10,
                    "triangle inequality violated: d[{v}]={} > d[{u}]={} + w={w}",
                    d[v],
                    d[u]
                );
            }
        }
    }
    #[test]
    fn test_dijkstra_multi_source_bounded() {
        let m = strip_mesh(5);
        let d = dijkstra_geodesic_multi_source(&m, &[0, 10]);
        for &di in &d {
            assert!(
                di.is_finite(),
                "all vertices should be reachable from two sources"
            );
        }
    }
    #[test]
    fn test_dijkstra_multi_source_min_of_single() {
        let m = tet_surface_mesh();
        let d0 = dijkstra_geodesic(&m, 0);
        let d1 = dijkstra_geodesic(&m, 1);
        let dm = dijkstra_geodesic_multi_source(&m, &[0, 1]);
        for i in 0..m.n_vertices() {
            let expected = d0[i].min(d1[i]);
            assert!(
                (dm[i] - expected).abs() < 1e-10,
                "vertex {i}: multi_source={} vs min({}, {})={}",
                dm[i],
                d0[i],
                d1[i],
                expected
            );
        }
    }
    #[test]
    fn test_dijkstra_predecessor_source_pred_is_max() {
        let m = single_triangle();
        let (_, pred) = dijkstra_with_predecessors(&m, 0);
        assert_eq!(pred[0], usize::MAX, "source predecessor should be MAX");
    }
    #[test]
    fn test_extract_path_same_vertex() {
        let pred = vec![usize::MAX, 0, 1];
        let path = extract_geodesic_path(0, 0, &pred);
        assert_eq!(path, Some(vec![0]));
    }
    #[test]
    fn test_extract_path_direct_edge() {
        let m = single_triangle();
        let (_, pred) = dijkstra_with_predecessors(&m, 0);
        let path = extract_geodesic_path(0, 2, &pred);
        assert!(path.is_some(), "path from 0 to 2 should exist");
        let path = path.unwrap();
        assert_eq!(path[0], 0, "path should start at 0");
        assert_eq!(*path.last().unwrap(), 2, "path should end at 2");
    }
    #[test]
    fn test_extract_path_unreachable() {
        let pred = vec![usize::MAX, usize::MAX];
        let path = extract_geodesic_path(0, 1, &pred);
        assert!(path.is_none(), "unreachable vertex should return None");
    }
    #[test]
    fn test_path_to_polyline_length() {
        let m = single_triangle();
        let (_, pred) = dijkstra_with_predecessors(&m, 0);
        let path = extract_geodesic_path(0, 1, &pred).unwrap();
        let poly = path_to_polyline(&m, &path);
        assert_eq!(poly.len(), path.len());
    }
    #[test]
    fn test_polyline_length_single_edge() {
        let poly = vec![[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]];
        let len = polyline_length(&poly);
        assert!((len - 5.0).abs() < 1e-10, "length should be 5.0, got {len}");
    }
    #[test]
    fn test_polyline_length_empty() {
        let poly: Vec<[f64; 3]> = Vec::new();
        assert_eq!(polyline_length(&poly), 0.0);
    }
    #[test]
    fn test_heat_geodesic_source_distance_zero() {
        let m = tet_surface_mesh();
        let result = heat_geodesic(&m, 0, 0.0, 3);
        assert_eq!(result.distances[0], 0.0, "source distance should be 0");
    }
    #[test]
    fn test_heat_geodesic_all_finite() {
        let m = tet_surface_mesh();
        let result = heat_geodesic(&m, 0, 0.0, 5);
        for (i, &d) in result.distances.iter().enumerate() {
            assert!(d.is_finite(), "distance to vertex {i} should be finite");
        }
    }
    #[test]
    fn test_heat_geodesic_non_negative() {
        let m = strip_mesh(3);
        let result = heat_geodesic(&m, 0, 0.1, 3);
        for &d in &result.distances {
            assert!(d >= 0.0, "distances should be non-negative");
        }
    }
    #[test]
    fn test_heat_geodesic_gradient_size() {
        let m = strip_mesh(3);
        let result = heat_geodesic(&m, 0, 0.1, 2);
        assert_eq!(
            result.gradient.len(),
            m.n_faces(),
            "gradient should have one entry per face"
        );
    }
    #[test]
    fn test_extract_isolines_single_triangle() {
        let m = single_triangle();
        let distances = vec![0.0, 1.0, 0.866_025];
        let segs = extract_isolines(&m, &distances, 0.5);
        assert!(!segs.is_empty(), "expected at least one isoline segment");
    }
    #[test]
    fn test_extract_isolines_outside_range_empty() {
        let m = single_triangle();
        let distances = vec![0.0, 0.3, 0.2];
        let segs = extract_isolines(&m, &distances, 1.0);
        assert!(
            segs.is_empty(),
            "no crossings expected when level is above all distances"
        );
    }
    #[test]
    fn test_voronoi_regions_single_source() {
        let m = tet_surface_mesh();
        let regions = geodesic_voronoi_regions(&m, &[0]);
        assert!(
            regions.iter().all(|&r| r == 0),
            "all regions should be 0 with one source"
        );
    }
    #[test]
    fn test_voronoi_regions_two_sources() {
        let m = strip_mesh(4);
        let regions = geodesic_voronoi_regions(&m, &[0, 8]);
        assert_eq!(regions[0], 0);
        assert_eq!(regions[8], 1);
    }
    #[test]
    fn test_voronoi_regions_length() {
        let m = tet_surface_mesh();
        let regions = geodesic_voronoi_regions(&m, &[0, 2]);
        assert_eq!(regions.len(), m.n_vertices());
    }
    #[test]
    fn test_geodesic_diameter_positive() {
        let m = strip_mesh(4);
        let (diam, u, v) = geodesic_diameter(&m);
        assert!(diam > 0.0, "diameter should be positive");
        assert_ne!(u, v, "endpoints of diameter should be different");
    }
    #[test]
    fn test_geodesic_diameter_empty_mesh() {
        let m = GeoMesh::new(Vec::new(), Vec::new());
        let (diam, _, _) = geodesic_diameter(&m);
        assert_eq!(diam, 0.0);
    }
    #[test]
    fn test_geodesic_diameter_single_vertex() {
        let m = GeoMesh::new(vec![[0.0, 0.0, 0.0]], Vec::new());
        let (diam, u, v) = geodesic_diameter(&m);
        assert_eq!(u, v);
        assert_eq!(diam, 0.0);
    }
    #[test]
    fn test_mean_geodesic_distance_non_negative() {
        let m = single_triangle();
        let mgd = mean_geodesic_distance(&m);
        for (i, &d) in mgd.iter().enumerate() {
            assert!(
                d >= 0.0,
                "mean geodesic distance for vertex {i} should be non-negative"
            );
        }
    }
    #[test]
    fn test_mean_geodesic_distance_equilateral_equal() {
        let m = single_triangle();
        let mgd = mean_geodesic_distance(&m);
        assert!(
            (mgd[0] - mgd[1]).abs() < 1e-10,
            "equilateral: all mean distances should be equal"
        );
        assert!((mgd[1] - mgd[2]).abs() < 1e-10);
    }
    #[test]
    fn test_eccentricity_non_negative() {
        let m = tet_surface_mesh();
        let ecc = geodesic_eccentricity(&m);
        for (i, &e) in ecc.iter().enumerate() {
            assert!(e >= 0.0, "eccentricity at {i} should be non-negative");
        }
    }
    #[test]
    fn test_eccentricity_length() {
        let m = strip_mesh(3);
        let ecc = geodesic_eccentricity(&m);
        assert_eq!(ecc.len(), m.n_vertices());
    }
    #[test]
    fn test_all_pairs_diagonal_zero() {
        let m = single_triangle();
        let d = all_pairs_geodesic(&m);
        let n = m.n_vertices();
        for i in 0..n {
            assert_eq!(d[i * n + i], 0.0, "diagonal should be 0");
        }
    }
    #[test]
    fn test_all_pairs_symmetric() {
        let m = single_triangle();
        let d = all_pairs_geodesic(&m);
        let n = m.n_vertices();
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (d[i * n + j] - d[j * n + i]).abs() < 1e-10,
                    "not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_all_pairs_triangle_inequality() {
        let m = single_triangle();
        let d = all_pairs_geodesic(&m);
        let n = m.n_vertices();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    assert!(
                        d[i * n + j] <= d[i * n + k] + d[k * n + j] + 1e-10,
                        "triangle inequality failed"
                    );
                }
            }
        }
    }
    #[test]
    fn test_geodesic_centroid_equilateral() {
        let m = single_triangle();
        let (centroid, _) = geodesic_centroid(&m);
        assert!(centroid < m.n_vertices());
    }
    #[test]
    fn test_geodesic_centroid_non_negative_sum() {
        let m = tet_surface_mesh();
        let (_, sum) = geodesic_centroid(&m);
        assert!(
            sum >= 0.0,
            "sum of squared distances should be non-negative"
        );
    }
    #[test]
    fn test_smooth_path_preserves_endpoints() {
        let path = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.5, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.3, 0.0],
        ];
        let smoothed = smooth_geodesic_path(&path, 5);
        assert!(
            (smoothed[0][0] - path[0][0]).abs() < 1e-10,
            "start should be fixed"
        );
        assert!(
            (smoothed.last().unwrap()[0] - path.last().unwrap()[0]).abs() < 1e-10,
            "end should be fixed"
        );
    }
    #[test]
    fn test_smooth_path_empty() {
        let path: Vec<[f64; 3]> = Vec::new();
        let smoothed = smooth_geodesic_path(&path, 3);
        assert!(smoothed.is_empty());
    }
    #[test]
    fn test_smooth_path_single() {
        let path = vec![[1.0, 2.0, 3.0]];
        let smoothed = smooth_geodesic_path(&path, 3);
        assert_eq!(smoothed.len(), 1);
    }
    #[test]
    fn test_smooth_path_finite_values() {
        let path = vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]];
        let smoothed = smooth_geodesic_path(&path, 10);
        for pt in &smoothed {
            for &c in pt {
                assert!(c.is_finite());
            }
        }
    }
    #[test]
    fn test_angular_defect_length() {
        let m = tet_surface_mesh();
        let defect = angular_defect(&m);
        assert_eq!(defect.len(), m.n_vertices());
    }
    #[test]
    fn test_angular_defect_finite() {
        let m = tet_surface_mesh();
        let defect = angular_defect(&m);
        for (i, &d) in defect.iter().enumerate() {
            assert!(
                d.is_finite(),
                "angular defect at vertex {i} should be finite"
            );
        }
    }
    #[test]
    fn test_total_gaussian_curvature_sphere_like() {
        let m = tet_surface_mesh();
        let total = total_gaussian_curvature(&m);
        let expected = 4.0 * std::f64::consts::PI;
        assert!(
            (total - expected).abs() < 0.5,
            "total Gaussian curvature should be ≈4π={:.4}, got {:.4}",
            expected,
            total
        );
    }
    #[test]
    fn test_angular_defect_single_triangle_boundary() {
        let m = single_triangle();
        let defect = angular_defect(&m);
        assert_eq!(defect.len(), 3);
        for &d in &defect {
            assert!(d.is_finite());
        }
    }
    #[test]
    fn test_angular_defect_sum_gauss_bonnet() {
        let m = tet_surface_mesh();
        let total = total_gaussian_curvature(&m);
        assert!(
            total > 0.0,
            "total curvature for closed convex mesh should be positive"
        );
    }
    #[test]
    fn test_cotangent_laplacian_length() {
        let m = tet_surface_mesh();
        let lap = build_cotangent_laplacian(&m);
        assert_eq!(lap.len(), m.n_vertices());
    }
    #[test]
    fn test_cotangent_laplacian_non_empty() {
        let m = tet_surface_mesh();
        let lap = build_cotangent_laplacian(&m);
        for (i, row) in lap.iter().enumerate() {
            assert!(!row.is_empty(), "vertex {i} should have cotangent weights");
        }
    }
    #[test]
    fn test_cotangent_weight_right_angle() {
        let cot = cotangent_weight(1.0, 1.0, 2.0_f64.sqrt());
        assert!(cot.abs() < 1e-9, "cot(π/2) should be 0, got {cot}");
    }
    #[test]
    fn test_cotangent_weight_equilateral() {
        let s = 1.0_f64;
        let cot = cotangent_weight(s, s, s);
        let expected = 1.0 / 3.0_f64.sqrt();
        assert!(
            (cot - expected).abs() < 1e-9,
            "cot(60°) should be ≈{expected}, got {cot}"
        );
    }
    #[test]
    fn test_intrinsic_delaunay_face_count_preserved() {
        let m = tet_surface_mesh();
        let id = intrinsic_delaunay(&m);
        assert_eq!(
            id.faces.len(),
            m.n_faces(),
            "face count should be preserved"
        );
    }
    #[test]
    fn test_intrinsic_delaunay_vertex_count_preserved() {
        let m = tet_surface_mesh();
        let id = intrinsic_delaunay(&m);
        assert_eq!(
            id.vertices.len(),
            m.n_vertices(),
            "vertex count should be preserved"
        );
    }
    #[test]
    fn test_intrinsic_delaunay_edge_lengths_positive() {
        let m = tet_surface_mesh();
        let id = intrinsic_delaunay(&m);
        for (fi, lengths) in id.edge_lengths.iter().enumerate() {
            for (k, &l) in lengths.iter().enumerate() {
                assert!(
                    l > 0.0,
                    "edge length [{fi}][{k}] should be positive, got {l}"
                );
            }
        }
    }
    #[test]
    fn test_intrinsic_delaunay_single_triangle() {
        let m = single_triangle();
        let id = intrinsic_delaunay(&m);
        assert_eq!(id.faces.len(), 1);
        assert_eq!(id.vertices.len(), 3);
    }
    #[test]
    fn test_dijkstra_intrinsic_source_zero() {
        let m = tet_surface_mesh();
        let id = intrinsic_delaunay(&m);
        let d = dijkstra_intrinsic(&id, 0);
        assert_eq!(d[0], 0.0, "source distance should be 0");
    }
    #[test]
    fn test_dijkstra_intrinsic_all_finite() {
        let m = tet_surface_mesh();
        let id = intrinsic_delaunay(&m);
        let d = dijkstra_intrinsic(&id, 0);
        for (i, &di) in d.iter().enumerate() {
            assert!(di.is_finite(), "vertex {i} should be reachable");
        }
    }
    #[test]
    fn test_dijkstra_intrinsic_non_negative() {
        let m = strip_mesh(4);
        let id = intrinsic_delaunay(&m);
        let d = dijkstra_intrinsic(&id, 0);
        for &di in &d {
            assert!(
                di >= 0.0 || !di.is_finite(),
                "distances must be non-negative"
            );
        }
    }
    #[test]
    fn test_parallel_transport_empty_path() {
        let m = single_triangle();
        let v = [1.0, 0.0, 0.0];
        let result = parallel_transport(&m, &[], v);
        assert_eq!(result, v, "empty path should return unchanged vector");
    }
    #[test]
    fn test_parallel_transport_single_vertex() {
        let m = single_triangle();
        let v = [1.0, 0.0, 0.0];
        let result = parallel_transport(&m, &[0], v);
        assert_eq!(
            result, v,
            "single vertex path should return unchanged vector"
        );
    }
    #[test]
    fn test_parallel_transport_finite_output() {
        let m = tet_surface_mesh();
        let (_, pred) = dijkstra_with_predecessors(&m, 0);
        let path = extract_geodesic_path(0, 3, &pred).unwrap();
        let v = [1.0, 0.0, 0.0];
        let result = parallel_transport(&m, &path, v);
        for &c in &result {
            assert!(c.is_finite(), "parallel transport output should be finite");
        }
    }
    #[test]
    fn test_varadhan_distances_source_is_zero() {
        let m = tet_surface_mesh();
        let d = varadhan_distances(&m, 0, 0.01);
        assert!(d[0] < 1e-6, "source distance should be ≈0, got {}", d[0]);
    }
    #[test]
    fn test_varadhan_distances_non_negative() {
        let m = strip_mesh(3);
        let d = varadhan_distances(&m, 0, 0.05);
        for &di in &d {
            assert!(di >= 0.0, "distances should be non-negative");
        }
    }
    #[test]
    fn test_varadhan_distances_length() {
        let m = tet_surface_mesh();
        let d = varadhan_distances(&m, 0, 0.1);
        assert_eq!(d.len(), m.n_vertices());
    }
    #[test]
    fn test_geo_voronoi_cells_count() {
        let m = tet_surface_mesh();
        let cells = geodesic_voronoi_cells(&m, &[0, 2]);
        assert_eq!(cells.len(), 2, "should return one cell per source");
    }
    #[test]
    fn test_geo_voronoi_cells_area_positive() {
        let m = strip_mesh(4);
        let cells = geodesic_voronoi_cells(&m, &[0, 8]);
        for (i, cell) in cells.iter().enumerate() {
            assert!(cell.area >= 0.0, "cell {i} area should be non-negative");
        }
    }
    #[test]
    fn test_geo_voronoi_cells_total_area() {
        let m = tet_surface_mesh();
        let cells = geodesic_voronoi_cells(&m, &[0, 1, 2, 3]);
        let total: f64 = cells.iter().map(|c| c.area).sum();
        let mesh_area: f64 = (0..m.n_faces()).map(|fi| m.face_area(fi)).sum();
        assert!(
            (total - mesh_area).abs() < 1e-9,
            "total cell area {total} should equal mesh area {mesh_area}"
        );
    }
    #[test]
    fn test_geo_voronoi_cells_vertices_coverage() {
        let m = tet_surface_mesh();
        let sources = vec![0, 1, 2, 3];
        let cells = geodesic_voronoi_cells(&m, &sources);
        let mut covered = vec![false; m.n_vertices()];
        for cell in &cells {
            for &v in &cell.vertices {
                covered[v] = true;
            }
        }
        assert!(covered.iter().all(|&c| c), "all vertices should be covered");
    }
    #[test]
    fn test_geo_voronoi_cells_empty_sources() {
        let m = tet_surface_mesh();
        let cells = geodesic_voronoi_cells(&m, &[]);
        assert!(cells.is_empty(), "empty sources should produce no cells");
    }
    #[test]
    fn test_geo_voronoi_cells_source_in_own_cell() {
        let m = tet_surface_mesh();
        let cells = geodesic_voronoi_cells(&m, &[0, 3]);
        assert!(
            cells[0].vertices.contains(&0),
            "source 0 should be in cell 0"
        );
        assert!(
            cells[1].vertices.contains(&3),
            "source 3 should be in cell 1"
        );
    }
}
