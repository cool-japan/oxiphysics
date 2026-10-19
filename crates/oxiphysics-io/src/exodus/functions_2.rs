//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_expanded {
    use crate::exodus::*;
    #[test]
    fn test_tet4_quality_unit_tet() {
        let mesh = build_unit_tet_mesh("q");
        let qual = tet4_quality(&mesh);
        assert_eq!(qual.len(), 1);
        let q = &qual[0];
        assert!(q.aspect_ratio > 0.0 && q.aspect_ratio <= 1.0);
        assert!(q.min_edge > 0.0);
        assert!(q.max_edge >= q.min_edge);
    }
    #[test]
    fn test_hex8_quality_unit_hex() {
        let mesh = build_unit_hex_mesh("q");
        let qual = hex8_quality(&mesh);
        assert_eq!(qual.len(), 1);
        let q = &qual[0];
        assert!((q.aspect_ratio - 1.0).abs() < 1e-10);
        assert!((q.scaled_jacobian - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_hex8_quality_skipped_for_non_hex() {
        let mesh = build_unit_tet_mesh("q");
        let qual = hex8_quality(&mesh);
        assert!(qual.is_empty(), "Should skip non-HEX8 blocks");
    }
    #[test]
    fn test_tet4_quality_skipped_for_non_tet() {
        let mesh = build_unit_hex_mesh("q");
        let qual = tet4_quality(&mesh);
        assert!(qual.is_empty(), "Should skip non-TET4 blocks");
    }
    #[test]
    fn test_node_distance() {
        let d = node_distance([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-12);
        let d2 = node_distance([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!((d2 - f64::sqrt(3.0)).abs() < 1e-12);
    }
    #[test]
    fn test_idw_exact_node() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        let field = vec![10.0, 20.0];
        let v = idw_interpolate(&mesh, &field, [0.0, 0.0, 0.0], 2.0);
        assert!((v - 10.0).abs() < 1e-6, "Got {}", v);
    }
    #[test]
    fn test_idw_midpoint() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(2.0, 0.0, 0.0);
        let field = vec![0.0, 1.0];
        let v = idw_interpolate(&mesh, &field, [1.0, 0.0, 0.0], 2.0);
        assert!((v - 0.5).abs() < 1e-6, "Got {}", v);
    }
    #[test]
    fn test_find_nearest_node() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(10.0, 0.0, 0.0);
        mesh.add_node(5.0, 0.0, 0.0);
        let idx = find_nearest_node(&mesh, [4.9, 0.0, 0.0]).unwrap();
        assert_eq!(idx, 2);
    }
    #[test]
    fn test_find_nearest_node_empty() {
        let mesh = ExodusMesh::new("empty");
        assert!(find_nearest_node(&mesh, [0.0, 0.0, 0.0]).is_none());
    }
    #[test]
    fn test_nodes_in_sphere() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        mesh.add_node(10.0, 0.0, 0.0);
        let inside = nodes_in_sphere(&mesh, [0.0, 0.0, 0.0], 1.5);
        assert!(inside.contains(&0));
        assert!(inside.contains(&1));
        assert!(!inside.contains(&2));
    }
    #[test]
    fn test_extract_block_unit_tet() {
        let mesh = build_unit_tet_mesh("orig");
        let sub = extract_block(&mesh, 1).unwrap();
        assert_eq!(sub.node_count(), 4);
        assert_eq!(sub.element_count(), 1);
    }
    #[test]
    fn test_extract_block_not_found() {
        let mesh = build_unit_tet_mesh("orig");
        assert!(extract_block(&mesh, 99).is_none());
    }
    #[test]
    fn test_extract_node_set_coords() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        mesh.add_node(2.0, 0.0, 0.0);
        mesh.add_node_set(1, "left", vec![1, 2]);
        let coords = extract_node_set_coords(&mesh, 1).unwrap();
        assert_eq!(coords.len(), 2);
        assert_eq!(coords[0], [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_refine_tet4_increases_elements() {
        let mesh = build_unit_tet_mesh("orig");
        let refined = refine_tet4_mesh(&mesh);
        assert_eq!(refined.element_count(), 8);
    }
    #[test]
    fn test_refine_tet4_adds_midpoint_nodes() {
        let mesh = build_unit_tet_mesh("orig");
        let refined = refine_tet4_mesh(&mesh);
        assert_eq!(refined.node_count(), 10);
    }
    #[test]
    fn test_all_boundary_nodes() {
        let mut mesh = ExodusMesh::new("t");
        for _ in 0..5 {
            mesh.add_node(0.0, 0.0, 0.0);
        }
        mesh.add_node_set(1, "left", vec![1, 2]);
        mesh.add_node_set(2, "right", vec![4, 5]);
        let bn = all_boundary_nodes(&mesh);
        assert_eq!(bn, vec![1, 2, 4, 5]);
    }
    #[test]
    fn test_node_valence() {
        let mut mesh = ExodusMesh::new("t");
        for _ in 0..4 {
            mesh.add_node(0.0, 0.0, 0.0);
        }
        let bi = mesh.add_block(1, "b", "TET4", 4);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        let valence = node_valence(&mesh);
        assert_eq!(valence.len(), 4);
        assert_eq!(valence[0], 1);
        assert_eq!(valence[3], 1);
    }
    #[test]
    fn test_write_read_result_text_roundtrip() {
        let mut mesh = ExodusMesh::new("rt_test");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        let mut result = ExodusResult {
            mesh,
            time_steps: Vec::new(),
        };
        add_node_variable_step(&mut result, 0.5, "T", vec![300.0, 310.0]);
        add_node_variable_step(&mut result, 1.0, "T", vec![305.0, 315.0]);
        let text = write_result_text(&result);
        let parsed = parse_result_text(&text).expect("parse failed");
        assert_eq!(parsed.time_steps.len(), 2);
        assert!((parsed.time_steps[0].time - 0.5).abs() < 1e-9);
        assert!((parsed.time_steps[1].time - 1.0).abs() < 1e-9);
        assert_eq!(parsed.time_steps[0].variables[0].name, "T");
        assert!((parsed.time_steps[0].variables[0].values[0] - 300.0).abs() < 1e-6);
    }
    #[test]
    fn test_write_result_text_no_timesteps() {
        let mesh = ExodusMesh::new("empty");
        let result = ExodusResult {
            mesh,
            time_steps: Vec::new(),
        };
        let text = write_result_text(&result);
        let parsed = parse_result_text(&text).expect("parse failed");
        assert!(parsed.time_steps.is_empty());
    }
    #[test]
    fn test_partition_round_robin_count() {
        let mut mesh = ExodusMesh::new("p");
        for _ in 0..4 {
            mesh.add_node(0.0, 0.0, 0.0);
        }
        let bi = mesh.add_block(1, "b", "TET4", 4);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        let parts = partition_mesh_round_robin(&mesh, 2);
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].mesh.element_count(), 1);
        assert_eq!(parts[1].mesh.element_count(), 1);
    }
    #[test]
    fn test_partition_zero_parts_empty() {
        let mesh = ExodusMesh::new("p");
        let parts = partition_mesh_round_robin(&mesh, 0);
        assert!(parts.is_empty());
    }
    #[test]
    fn test_field_stats_uniform() {
        let v = vec![3.0_f64; 5];
        let s = field_stats(&v);
        assert!((s.min - 3.0).abs() < 1e-12);
        assert!((s.max - 3.0).abs() < 1e-12);
        assert!((s.mean - 3.0).abs() < 1e-12);
        assert!(s.std_dev.abs() < 1e-12);
        assert_eq!(s.count, 5);
    }
    #[test]
    fn test_field_stats_varied() {
        let v = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let s = field_stats(&v);
        assert!((s.min - 1.0).abs() < 1e-12);
        assert!((s.max - 5.0).abs() < 1e-12);
        assert!((s.mean - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_field_stats_empty() {
        let s = field_stats(&[]);
        assert_eq!(s.count, 0);
    }
    #[test]
    fn test_normalize_field() {
        let v = vec![0.0_f64, 5.0, 10.0];
        let n = normalize_field(&v);
        assert!((n[0] - 0.0).abs() < 1e-12);
        assert!((n[1] - 0.5).abs() < 1e-12);
        assert!((n[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalize_field_uniform() {
        let v = vec![7.0_f64; 4];
        let n = normalize_field(&v);
        for &x in &n {
            assert!((x - 0.5).abs() < 1e-12);
        }
    }
    #[test]
    fn test_field_l2_norm() {
        let v = vec![3.0_f64, 4.0];
        let norm = field_l2_norm(&v);
        assert!((norm - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_export_nodes_csv_header() {
        let mesh = build_unit_tet_mesh("t");
        let csv = export_nodes_csv(&mesh);
        assert!(csv.starts_with("node_id,x,y,z\n"));
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 5);
    }
    #[test]
    fn test_export_field_csv() {
        let mut mesh = ExodusMesh::new("t");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        let field = vec![100.0, 200.0];
        let csv = export_field_csv(&mesh, "temp", &field);
        assert!(csv.contains("temp"));
        assert!(csv.contains("100"));
        assert!(csv.contains("200"));
    }
    #[test]
    fn test_export_connectivity_csv() {
        let mesh = build_unit_tet_mesh("t");
        let csv = export_connectivity_csv(&mesh);
        assert!(csv.contains("block_id"));
        assert!(csv.contains("element_id"));
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 2);
    }
}
#[cfg(test)]
mod tests_exodus_writer {
    use crate::exodus::*;
    #[test]
    fn test_nodal_var_step_contains_name() {
        let s = ExodusWriter::write_nodal_variable_step("temperature", 0, 0.0, &[300.0, 310.0]);
        assert!(s.contains("temperature"), "Should contain var name");
    }
    #[test]
    fn test_nodal_var_step_contains_step_and_time() {
        let s = ExodusWriter::write_nodal_variable_step("pressure", 3, 1.5e-3, &[1.0]);
        assert!(s.contains("STEP 3"), "Should contain step index");
        assert!(s.contains("TIME"), "Should contain TIME keyword");
    }
    #[test]
    fn test_nodal_var_step_values_in_output() {
        let vals = vec![1.0_f64, 2.0, 3.0];
        let s = ExodusWriter::write_nodal_variable_step("vel_x", 0, 0.0, &vals);
        assert!(
            s.contains("1.00000000e0") || s.contains("1.00000000e") || s.contains("1.0"),
            "Values should appear in output"
        );
    }
    #[test]
    fn test_nodal_var_step_ends_with_end_marker() {
        let s = ExodusWriter::write_nodal_variable_step("disp", 1, 0.01, &[0.1]);
        assert!(
            s.trim_end().ends_with("END_NODAL_VAR"),
            "Should end with END_NODAL_VAR"
        );
    }
    #[test]
    fn test_global_var_contains_name() {
        let v = ExodusGlobalVariable {
            name: "kinetic_energy".to_string(),
            time: 0.5,
            value: 1234.5,
        };
        let s = ExodusWriter::write_global_variable(&v);
        assert!(s.contains("kinetic_energy"), "Should contain var name");
    }
    #[test]
    fn test_global_var_contains_time_and_value_keywords() {
        let v = ExodusGlobalVariable {
            name: "total_mass".to_string(),
            time: 0.1,
            value: 99.9,
        };
        let s = ExodusWriter::write_global_variable(&v);
        assert!(s.contains("TIME"), "Should contain TIME");
        assert!(s.contains("VALUE"), "Should contain VALUE");
    }
    #[test]
    fn test_global_var_single_line() {
        let v = ExodusGlobalVariable {
            name: "x".to_string(),
            time: 0.0,
            value: 0.0,
        };
        let s = ExodusWriter::write_global_variable(&v);
        assert_eq!(s.lines().count(), 1, "Global var should be a single line");
    }
    #[test]
    fn test_side_set_contains_name_and_id() {
        let data = ExodusSideSetData {
            id: 10,
            name: "inlet".to_string(),
            element_ids: vec![1, 2],
            side_ids: vec![3, 4],
            dist_factors: vec![],
        };
        let s = ExodusWriter::write_side_set(&data);
        assert!(s.contains("inlet"), "Should contain set name");
        assert!(s.contains("10"), "Should contain set id");
    }
    #[test]
    fn test_side_set_element_side_pairs_in_output() {
        let data = ExodusSideSetData {
            id: 1,
            name: "wall".to_string(),
            element_ids: vec![5, 6],
            side_ids: vec![2, 3],
            dist_factors: vec![],
        };
        let s = ExodusWriter::write_side_set(&data);
        assert!(s.contains("5 2"), "Should contain element/side pair 5 2");
        assert!(s.contains("6 3"), "Should contain element/side pair 6 3");
    }
    #[test]
    fn test_side_set_dist_factors_written_when_present() {
        let data = ExodusSideSetData {
            id: 2,
            name: "outlet".to_string(),
            element_ids: vec![1],
            side_ids: vec![1],
            dist_factors: vec![0.5],
        };
        let s = ExodusWriter::write_side_set(&data);
        assert!(
            s.contains("DIST_FACTORS"),
            "Should contain DIST_FACTORS section"
        );
    }
}
#[cfg(test)]
mod tests_refine_tet4 {
    use crate::exodus::*;

    fn tet_mesh(vertices: &[[f64; 3]; 4]) -> ExodusMesh {
        let mut mesh = ExodusMesh::new("refine");
        for v in vertices {
            mesh.add_node(v[0], v[1], v[2]);
        }
        let bi = mesh.add_block(1, "tet_block", "TET4", 4);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        mesh
    }

    /// Signed volume of the tetrahedron with 1-based node ids `ids`.
    fn signed_volume(mesh: &ExodusMesh, ids: &[usize]) -> f64 {
        let p = |k: usize| mesh.coordinates[ids[k] - 1];
        let (p0, p1, p2, p3) = (p(0), p(1), p(2), p(3));
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let e3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
        (e1[0] * (e2[1] * e3[2] - e2[2] * e3[1]) - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
            + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]))
            / 6.0
    }

    fn element_volumes(mesh: &ExodusMesh) -> Vec<f64> {
        let mut out = Vec::new();
        for block in mesh.blocks.iter().filter(|b| b.element_type == "TET4") {
            for ids in block.connectivity.chunks_exact(4) {
                out.push(signed_volume(mesh, ids));
            }
        }
        out
    }

    /// 1-based id of the node of `mesh` located at `target`.
    fn node_at(mesh: &ExodusMesh, target: [f64; 3]) -> usize {
        let found = mesh.coordinates.iter().position(|c| {
            (c[0] - target[0]).abs() < 1e-12
                && (c[1] - target[1]).abs() < 1e-12
                && (c[2] - target[2]).abs() < 1e-12
        });
        found.map(|i| i + 1).unwrap_or(0)
    }

    fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            0.5 * (a[0] + b[0]),
            0.5 * (a[1] + b[1]),
            0.5 * (a[2] + b[2]),
        ]
    }

    fn distance_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
        (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
    }

    /// Number of elements of `mesh` that contain both node ids.
    fn elements_containing_both(mesh: &ExodusMesh, a: usize, b: usize) -> usize {
        mesh.blocks
            .iter()
            .flat_map(|blk| blk.connectivity.chunks_exact(4))
            .filter(|ids| ids.contains(&a) && ids.contains(&b))
            .count()
    }

    const IRREGULAR_TET: [[f64; 3]; 4] = [
        [0.1, -0.2, 0.3],
        [2.3, 0.4, -0.1],
        [0.5, 1.9, 0.2],
        [-0.3, 0.6, 2.2],
    ];

    #[test]
    fn unit_tet_children_all_have_positive_volume() {
        let mesh = build_unit_tet_mesh("unit");
        let parent = signed_volume(&mesh, &[1, 2, 3, 4]);
        assert!((parent - 1.0 / 6.0).abs() < 1e-15);
        let refined = refine_tet4_mesh(&mesh);
        let volumes = element_volumes(&refined);
        assert_eq!(volumes.len(), 8);
        for (i, v) in volumes.iter().enumerate() {
            assert!(*v > 0.0, "child {i} is inverted: signed volume {v}");
            assert!(
                (v - parent / 8.0).abs() < 1e-15,
                "child {i} volume {v} differs from parent/8 = {}",
                parent / 8.0
            );
        }
        for q in tet4_quality(&refined) {
            assert!(
                q.scaled_jacobian > 0.0,
                "element {} has non-positive scaled Jacobian {}",
                q.element_index,
                q.scaled_jacobian
            );
        }
    }

    #[test]
    fn irregular_tet_children_are_positive_and_sum_to_parent() {
        let mesh = tet_mesh(&IRREGULAR_TET);
        let parent = signed_volume(&mesh, &[1, 2, 3, 4]);
        assert!(parent > 0.0, "fixture must be positively oriented");
        let refined = refine_tet4_mesh(&mesh);
        let volumes = element_volumes(&refined);
        assert_eq!(volumes.len(), 8);
        for (i, v) in volumes.iter().enumerate() {
            assert!(*v > 0.0, "child {i} is inverted: signed volume {v}");
            assert!(
                (v - parent / 8.0).abs() < 1e-12,
                "child {i} volume {v} differs from parent/8 = {}",
                parent / 8.0
            );
        }
        let total: f64 = volumes.iter().sum();
        assert!(
            (total - parent).abs() < 1e-12,
            "children sum to {total}, parent is {parent}"
        );
        for q in tet4_quality(&refined) {
            assert!(q.scaled_jacobian > 0.0);
        }
    }

    #[test]
    fn inverted_parent_yields_inverted_children() {
        let mut swapped = IRREGULAR_TET;
        swapped.swap(1, 2);
        let mesh = tet_mesh(&swapped);
        let parent = signed_volume(&mesh, &[1, 2, 3, 4]);
        assert!(parent < 0.0, "fixture must be negatively oriented");
        let refined = refine_tet4_mesh(&mesh);
        let volumes = element_volumes(&refined);
        assert_eq!(volumes.len(), 8);
        for (i, v) in volumes.iter().enumerate() {
            assert!(*v < 0.0, "child {i} flipped orientation: volume {v}");
            assert!((v - parent / 8.0).abs() < 1e-12);
        }
    }

    #[test]
    fn interior_split_follows_the_shortest_diagonal() {
        // Each fixture makes a different octahedron diagonal the shortest one:
        // the pairs of opposite edges are (v0v2, v1v3), (v0v1, v2v3) and
        // (v0v3, v1v2), so a diagonal is short when the sums of the two vertex
        // pairs nearly coincide.
        let fixtures: [([[f64; 3]; 4], usize); 3] = [
            (
                [
                    [0.0, 0.0, 0.0],
                    [1.0, -1.0, -0.5],
                    [2.0, 0.0, 0.0],
                    [1.0, 1.0, 1.0],
                ],
                0,
            ),
            (
                [
                    [0.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [1.0, 1.0, 1.0],
                    [1.0, -1.0, -0.5],
                ],
                1,
            ),
            (
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 1.0, 1.0],
                    [1.0, -1.0, -0.5],
                    [2.0, 0.0, 0.0],
                ],
                2,
            ),
        ];
        for (vertices, expected) in fixtures {
            let mesh = tet_mesh(&vertices);
            let parent = signed_volume(&mesh, &[1, 2, 3, 4]);
            assert!(
                parent > 0.0,
                "fixture {expected} must be positively oriented"
            );
            let m = |a: usize, b: usize| midpoint(vertices[a], vertices[b]);
            let diagonals = [(m(0, 2), m(1, 3)), (m(0, 1), m(2, 3)), (m(0, 3), m(1, 2))];
            let lengths: Vec<f64> = diagonals.iter().map(|&(a, b)| distance_sq(a, b)).collect();
            for (idx, len) in lengths.iter().enumerate() {
                if idx != expected {
                    assert!(
                        lengths[expected] < *len - 1e-9,
                        "fixture {expected}: diagonal {idx} ({len}) is not longer than {}",
                        lengths[expected]
                    );
                }
            }
            let refined = refine_tet4_mesh(&mesh);
            for (idx, &(a, b)) in diagonals.iter().enumerate() {
                let (ia, ib) = (node_at(&refined, a), node_at(&refined, b));
                assert!(ia > 0 && ib > 0, "midpoint node missing in refined mesh");
                let shared = elements_containing_both(&refined, ia, ib);
                if idx == expected {
                    assert_eq!(
                        shared, 4,
                        "shortest diagonal {idx} must be an edge of 4 children"
                    );
                } else {
                    assert_eq!(shared, 0, "longer diagonal {idx} must not be an edge");
                }
            }
            let volumes = element_volumes(&refined);
            assert_eq!(volumes.len(), 8);
            for v in &volumes {
                assert!(
                    (v - parent / 8.0).abs() < 1e-12,
                    "fixture {expected}: child volume {v}, parent/8 = {}",
                    parent / 8.0
                );
            }
        }
    }

    #[test]
    fn neighbouring_tets_share_midpoints_and_stay_positive() {
        let mut mesh = ExodusMesh::new("pair");
        mesh.add_node(0.0, 0.0, 0.0);
        mesh.add_node(1.0, 0.0, 0.0);
        mesh.add_node(0.0, 1.0, 0.0);
        mesh.add_node(0.0, 0.0, 1.0);
        mesh.add_node(0.0, 0.0, -1.0);
        let bi = mesh.add_block(1, "tets", "TET4", 4);
        mesh.add_element_to_block(bi, &[1, 2, 3, 4]);
        mesh.add_element_to_block(bi, &[1, 3, 2, 5]);
        let before = element_volumes(&mesh);
        assert!(before.iter().all(|v| *v > 0.0));
        let refined = refine_tet4_mesh(&mesh);
        assert_eq!(refined.element_count(), 16);
        // 5 original nodes + (6 + 6 - 3 shared) edge midpoints.
        assert_eq!(refined.node_count(), 14);
        for i in 0..refined.coordinates.len() {
            for j in (i + 1)..refined.coordinates.len() {
                assert!(
                    distance_sq(refined.coordinates[i], refined.coordinates[j]) > 1e-18,
                    "nodes {i} and {j} coincide"
                );
            }
        }
        let after = element_volumes(&refined);
        assert!(after.iter().all(|v| *v > 0.0));
        let total_before: f64 = before.iter().sum();
        let total_after: f64 = after.iter().sum();
        assert!((total_before - total_after).abs() < 1e-12);
    }
}
