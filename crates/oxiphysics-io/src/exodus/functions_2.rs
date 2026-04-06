//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use crate::exodus::types::*;
#[cfg(test)]
mod tests_expanded {
    use super::*;
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
    use super::*;
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
