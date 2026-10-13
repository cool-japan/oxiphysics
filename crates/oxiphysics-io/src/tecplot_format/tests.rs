// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the Tecplot ASCII reader and writer.

use super::*;

// ── TecplotZoneType ─────────────────────────────────────────────────────

#[test]
fn test_zone_type_as_str_ordered() {
    assert_eq!(TecplotZoneType::Ordered.as_tecplot_str(), "ORDERED");
}

#[test]
fn test_zone_type_as_str_fe_triangle() {
    assert_eq!(TecplotZoneType::FETriangle.as_tecplot_str(), "FETRIANGLE");
}

#[test]
fn test_zone_type_as_str_fe_quad() {
    assert_eq!(TecplotZoneType::FEQuad.as_tecplot_str(), "FEQUADRILATERAL");
}

#[test]
fn test_zone_type_as_str_fe_tetra() {
    assert_eq!(TecplotZoneType::FETetra.as_tecplot_str(), "FETETRAHEDRON");
}

#[test]
fn test_zone_type_as_str_fe_brick() {
    assert_eq!(TecplotZoneType::FEBrick.as_tecplot_str(), "FEBRICK");
}

#[test]
fn test_zone_type_from_str_ordered() {
    assert_eq!(
        TecplotZoneType::from_keyword("ORDERED"),
        TecplotZoneType::Ordered
    );
}

#[test]
fn test_zone_type_from_str_fe_triangle() {
    assert_eq!(
        TecplotZoneType::from_keyword("FETRIANGLE"),
        TecplotZoneType::FETriangle
    );
}

#[test]
fn test_zone_type_from_str_fe_tetra() {
    assert_eq!(
        TecplotZoneType::from_keyword("FETETRAHEDRON"),
        TecplotZoneType::FETetra
    );
}

#[test]
fn test_zone_type_from_str_case_insensitive() {
    assert_eq!(
        TecplotZoneType::from_keyword("ordered"),
        TecplotZoneType::Ordered
    );
}

#[test]
fn test_zone_type_from_str_unknown_defaults_to_ordered() {
    assert_eq!(
        TecplotZoneType::from_keyword("UNKNOWN"),
        TecplotZoneType::Ordered
    );
}

// ── TecplotVariable ─────────────────────────────────────────────────────

#[test]
fn test_variable_new() {
    let v = TecplotVariable::new("Pressure", vec![1.0, 2.0, 3.0]);
    assert_eq!(v.name, "Pressure");
    assert_eq!(v.len(), 3);
    assert!(!v.is_empty());
}

#[test]
fn test_variable_empty() {
    let v = TecplotVariable::new("X", vec![]);
    assert!(v.is_empty());
}

#[test]
fn test_variable_data_values() {
    let v = TecplotVariable::new("T", vec![100.0, 200.0]);
    assert!((v.data[0] - 100.0).abs() < 1e-12);
    assert!((v.data[1] - 200.0).abs() < 1e-12);
}

// ── TecplotZone ─────────────────────────────────────────────────────────

#[test]
fn test_zone_new_ordered() {
    let z = TecplotZone::new_ordered("Zone 1", 5, 3, 2);
    assert_eq!(z.zone_type, TecplotZoneType::Ordered);
    assert_eq!(z.i_dim, 5);
    assert_eq!(z.j_dim, 3);
    assert_eq!(z.k_dim, 2);
    assert_eq!(z.point_count(), 30);
}

#[test]
fn test_zone_new_fe() {
    let connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7];
    let z = TecplotZone::new_fe(
        "FE Zone",
        TecplotZoneType::FETetra,
        8,
        2,
        connectivity.clone(),
    );
    assert_eq!(z.zone_type, TecplotZoneType::FETetra);
    assert_eq!(z.n_nodes, 8);
    assert_eq!(z.n_elements, 2);
    assert_eq!(z.point_count(), 8);
    assert_eq!(z.connectivity, connectivity);
}

#[test]
fn test_zone_type_nodes_per_element() {
    assert_eq!(TecplotZoneType::Ordered.nodes_per_element(), 0);
    assert_eq!(TecplotZoneType::FETriangle.nodes_per_element(), 3);
    assert_eq!(TecplotZoneType::FEQuad.nodes_per_element(), 4);
    assert_eq!(TecplotZoneType::FETetra.nodes_per_element(), 4);
    assert_eq!(TecplotZoneType::FEBrick.nodes_per_element(), 8);
}

#[test]
fn test_zone_type_from_element_type() {
    assert_eq!(
        TecplotZoneType::from_element_type("Triangle"),
        Some(TecplotZoneType::FETriangle)
    );
    assert_eq!(
        TecplotZoneType::from_element_type("QUADRILATERAL"),
        Some(TecplotZoneType::FEQuad)
    );
    assert_eq!(
        TecplotZoneType::from_element_type("tetrahedron"),
        Some(TecplotZoneType::FETetra)
    );
    assert_eq!(
        TecplotZoneType::from_element_type("BRICK"),
        Some(TecplotZoneType::FEBrick)
    );
    assert_eq!(TecplotZoneType::from_element_type("POLYGON"), None);
}

#[test]
fn test_zone_element_nodes_accessor() {
    let z = TecplotZone::new_fe(
        "Z",
        TecplotZoneType::FETriangle,
        4,
        2,
        vec![0, 1, 2, 1, 3, 2],
    );
    assert_eq!(z.element_nodes(0), Some(&[0usize, 1, 2][..]));
    assert_eq!(z.element_nodes(1), Some(&[1usize, 3, 2][..]));
    assert_eq!(z.element_nodes(2), None);
    let ordered = TecplotZone::new_ordered("O", 2, 1, 1);
    assert_eq!(ordered.element_nodes(0), None);
    // A connectivity array shorter than the element count is not indexed
    // out of range.
    let short = TecplotZone::new_fe("S", TecplotZoneType::FEQuad, 4, 2, vec![0, 1, 2, 3]);
    assert_eq!(short.element_nodes(0), Some(&[0usize, 1, 2, 3][..]));
    assert_eq!(short.element_nodes(1), None);
}

#[test]
fn test_zone_point_count_ordered_1d() {
    let z = TecplotZone::new_ordered("1D", 10, 1, 1);
    assert_eq!(z.point_count(), 10);
}

#[test]
fn test_zone_stores_variables() {
    let mut z = TecplotZone::new_ordered("Z", 2, 1, 1);
    z.variables.push(TecplotVariable::new("X", vec![0.0, 1.0]));
    assert_eq!(z.variables.len(), 1);
    assert_eq!(z.variables[0].name, "X");
}

// ── TecplotDataset ──────────────────────────────────────────────────────

#[test]
fn test_dataset_new() {
    let ds = TecplotDataset::new("Test");
    assert_eq!(ds.title, "Test");
    assert!(ds.zones.is_empty());
    assert!(ds.variables.is_empty());
}

#[test]
fn test_dataset_default() {
    let ds = TecplotDataset::default();
    assert!(ds.title.is_empty());
}

// ── TecplotWriter ───────────────────────────────────────────────────────

fn make_1d_dataset(path: &str, n: usize) -> TecplotDataset {
    let _ = path; // unused in builder
    let mut ds = TecplotDataset::new("Test Dataset");
    ds.variables = vec!["X".to_string(), "P".to_string()];
    let mut zone = TecplotZone::new_ordered("Zone 1", n, 1, 1);
    let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let p: Vec<f64> = (0..n).map(|i| (i as f64) * 2.0).collect();
    zone.variables.push(TecplotVariable::new("X", x));
    zone.variables.push(TecplotVariable::new("P", p));
    ds.zones.push(zone);
    ds
}

#[test]
fn test_write_creates_file() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_write_test.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 5);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .expect("write failed");
    assert!(path.exists());
}

#[test]
fn test_write_contains_title() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_title.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("Test Dataset"), "no title in output");
}

#[test]
fn test_write_contains_variables() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_vars.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("VARIABLES"), "no VARIABLES header");
    assert!(content.contains('"'), "variables not quoted");
}

#[test]
fn test_write_contains_zone_keyword() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_zone_kw.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("ZONE"), "no ZONE header");
}

#[test]
fn test_write_ordered_has_ijk() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_ijk.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("I="), "no I= in ZONE header");
    assert!(content.contains("J="), "no J= in ZONE header");
    assert!(content.contains("K="), "no K= in ZONE header");
}

#[test]
fn test_write_data_lines_count() {
    // A 1D ordered zone with I=4 should produce exactly 4 data rows
    let path = std::env::temp_dir().join("oxiphysics_tecplot_datarows.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 4);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    // Count lines that start with a digit or '-' (data lines)
    let data_lines = content
        .lines()
        .filter(|l| {
            let t = l.trim();
            t.starts_with(|c: char| c.is_ascii_digit() || c == '-')
        })
        .count();
    assert_eq!(data_lines, 4, "expected 4 data rows, got {data_lines}");
}

#[test]
fn test_write_fe_zone() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_fe_zone.dat");
    let mut ds = TecplotDataset::new("FE Dataset");
    ds.variables = vec!["X".to_string()];
    let mut zone = TecplotZone::new_fe("FE Zone", TecplotZoneType::FETetra, 4, 1, vec![0, 1, 2, 3]);
    zone.variables
        .push(TecplotVariable::new("X", vec![0.0, 1.0, 0.0, 0.0]));
    ds.zones.push(zone);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("FETETRAHEDRON"));
    fs::remove_file(&path).ok();
}

// ── Write/Read roundtrip ────────────────────────────────────────────────

#[test]
fn test_roundtrip_title() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_title.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(parsed.title, "Test Dataset");
}

#[test]
fn test_roundtrip_variable_names() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_varnames.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(parsed.variables, vec!["X", "P"]);
}

#[test]
fn test_roundtrip_zone_count() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_nzones.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(parsed.zones.len(), 1, "expected 1 zone");
}

#[test]
fn test_roundtrip_zone_type_ordered() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_ztype.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(parsed.zones[0].zone_type, TecplotZoneType::Ordered);
}

#[test]
fn test_roundtrip_zone_i_dim() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_idim.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 5);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(parsed.zones[0].i_dim, 5);
}

#[test]
fn test_roundtrip_data_values() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_data.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    let zone = &parsed.zones[0];
    // X variable should be [0,1,2]
    let x_var = zone.variables.iter().find(|v| v.name == "X").unwrap();
    assert!((x_var.data[0]).abs() < 1e-6);
    assert!((x_var.data[1] - 1.0).abs() < 1e-6);
    assert!((x_var.data[2] - 2.0).abs() < 1e-6);
}

#[test]
fn test_roundtrip_second_variable() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_p.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 3);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    let zone = &parsed.zones[0];
    let p_var = zone.variables.iter().find(|v| v.name == "P").unwrap();
    assert!((p_var.data[0]).abs() < 1e-6, "P[0]={}", p_var.data[0]);
    assert!((p_var.data[1] - 2.0).abs() < 1e-6, "P[1]={}", p_var.data[1]);
    assert!((p_var.data[2] - 4.0).abs() < 1e-6, "P[2]={}", p_var.data[2]);
}

#[test]
fn test_roundtrip_multi_zone() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_rt_mz.dat");
    let mut ds = TecplotDataset::new("Multi");
    ds.variables = vec!["X".to_string()];
    for i in 0..3usize {
        let mut z = TecplotZone::new_ordered(format!("Zone {i}"), 2, 1, 1);
        z.variables
            .push(TecplotVariable::new("X", vec![i as f64, i as f64 + 1.0]));
        ds.zones.push(z);
    }
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(parsed.zones.len(), 3, "expected 3 zones");
}

#[test]
fn test_read_missing_file_returns_error() {
    let path = std::env::temp_dir().join("no_such_file_oxiphysics_tec.dat");
    let result = TecplotDataset::read(path.to_str().unwrap_or(""));
    assert!(result.is_err());
}

#[test]
fn test_roundtrip_large_ordered_zone() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_large.dat");
    let n = 50;
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), n);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    let zone = &parsed.zones[0];
    assert_eq!(zone.i_dim, n);
    let x_var = zone.variables.iter().find(|v| v.name == "X").unwrap();
    assert_eq!(x_var.data.len(), n);
}

#[test]
fn test_roundtrip_3d_zone() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_3d.dat");
    let mut ds = TecplotDataset::new("3D");
    ds.variables = vec!["X".to_string()];
    let n = 2usize;
    let total = n * n * n;
    let mut zone = TecplotZone::new_ordered("3D Zone", n, n, n);
    let x: Vec<f64> = (0..total).map(|i| i as f64).collect();
    zone.variables.push(TecplotVariable::new("X", x));
    ds.zones.push(zone);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    assert_eq!(parsed.zones[0].i_dim, n);
    assert_eq!(parsed.zones[0].j_dim, n);
    assert_eq!(parsed.zones[0].k_dim, n);
}

// ── TecplotWriter helpers ────────────────────────────────────────────────

#[test]
fn test_writer_datapacking_point_keyword() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_dpkw.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 2);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("DATAPACKING=POINT"));
}

#[test]
fn test_writer_zone_title_in_output() {
    let path = std::env::temp_dir().join("oxiphysics_tecplot_ztitle.dat");
    let ds = make_1d_dataset(path.to_str().unwrap_or(""), 2);
    TecplotWriter::new()
        .write(&ds, path.to_str().unwrap_or(""))
        .unwrap();
    let content = fs::read_to_string(&path).unwrap();
    assert!(content.contains("Zone 1"), "zone title missing");
}

// ── Variable line parser ────────────────────────────────────────────────

#[test]
fn test_parse_variables_line_three_vars() {
    let line = r#"VARIABLES = "X" "Y" "Z""#;
    let vars = TecplotReader::parse_variables_line(line);
    assert_eq!(vars, vec!["X", "Y", "Z"]);
}

#[test]
fn test_parse_variables_line_with_commas() {
    let line = r#"VARIABLES = "X", "P""#;
    let vars = TecplotReader::parse_variables_line(line);
    assert_eq!(vars, vec!["X", "P"]);
}

#[test]
fn test_parse_variables_line_single() {
    let line = r#"VARIABLES = "T""#;
    let vars = TecplotReader::parse_variables_line(line);
    assert_eq!(vars, vec!["T"]);
}

// ── FE zone connectivity (F2496) ────────────────────────────────────────

/// Dataset with one FE zone of the given type: `n_nodes` nodes on a line
/// carrying variables X and P, `connectivity` as supplied.
fn fe_dataset(
    zone_type: TecplotZoneType,
    n_nodes: usize,
    connectivity: Vec<usize>,
) -> TecplotDataset {
    let npe = zone_type.nodes_per_element();
    let n_elements = connectivity.len() / npe;
    let mut ds = TecplotDataset::new("FE");
    ds.variables = vec!["X".to_string(), "P".to_string()];
    let mut zone = TecplotZone::new_fe("Mesh", zone_type, n_nodes, n_elements, connectivity);
    let x: Vec<f64> = (0..n_nodes).map(|i| i as f64 * 0.5).collect();
    let p: Vec<f64> = (0..n_nodes).map(|i| 100.0 + i as f64).collect();
    zone.variables.push(TecplotVariable::new("X", x));
    zone.variables.push(TecplotVariable::new("P", p));
    ds.zones.push(zone);
    ds
}

fn temp_path(name: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!("oxiphysics_tecplot_{name}.dat"))
}

#[test]
fn test_fe_two_triangles_roundtrip_connectivity() {
    // 4 nodes, 2 triangles sharing the edge (1, 2).
    let path = temp_path("fe_two_triangles");
    let ds = fe_dataset(TecplotZoneType::FETriangle, 4, vec![0, 1, 2, 1, 3, 2]);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    fs::remove_file(&path).ok();

    assert_eq!(parsed.zones.len(), 1);
    let zone = &parsed.zones[0];
    assert_eq!(zone.zone_type, TecplotZoneType::FETriangle);
    assert_eq!(zone.n_nodes, 4);
    assert_eq!(zone.n_elements, 2);
    assert_eq!(zone.connectivity, vec![0, 1, 2, 1, 3, 2]);
    assert_eq!(zone.element_nodes(1), Some(&[1usize, 3, 2][..]));
    // Nodal data is not disturbed by the connectivity block.
    let x = zone.variables.iter().find(|v| v.name == "X").unwrap();
    let p = zone.variables.iter().find(|v| v.name == "P").unwrap();
    assert_eq!(x.data.len(), 4);
    assert_eq!(p.data.len(), 4);
    assert!((x.data[3] - 1.5).abs() < 1e-12);
    assert!((p.data[2] - 102.0).abs() < 1e-12);
}

#[test]
fn test_fe_written_file_has_one_based_connectivity_after_node_rows() {
    let path = temp_path("fe_layout");
    let ds = fe_dataset(TecplotZoneType::FETriangle, 4, vec![0, 1, 2, 1, 3, 2]);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let content = fs::read_to_string(&path).unwrap();
    fs::remove_file(&path).ok();

    let lines: Vec<&str> = content.lines().collect();
    let zone_line = lines.iter().position(|l| l.starts_with("ZONE")).unwrap();
    assert!(lines[zone_line].contains("N=4"));
    assert!(lines[zone_line].contains("E=2"));
    assert!(lines[zone_line].contains("ZONETYPE=FETRIANGLE"));
    // 4 node rows with 2 variables each, then 2 element rows.
    let body = &lines[zone_line + 1..];
    assert_eq!(body.len(), 6, "body: {body:?}");
    for row in &body[..4] {
        assert_eq!(row.split_whitespace().count(), 2, "node row: {row}");
    }
    assert_eq!(body[4], "1 2 3");
    assert_eq!(body[5], "2 4 3");
}

#[test]
fn test_fe_all_zone_types_roundtrip() {
    let cases: Vec<(TecplotZoneType, usize, Vec<usize>)> = vec![
        (TecplotZoneType::FEQuad, 6, vec![0, 1, 4, 3, 1, 2, 5, 4]),
        (TecplotZoneType::FETetra, 5, vec![0, 1, 2, 3, 1, 2, 3, 4]),
        (
            TecplotZoneType::FEBrick,
            12,
            vec![0, 1, 3, 2, 6, 7, 9, 8, 4, 5, 11, 10, 0, 1, 2, 3],
        ),
    ];
    for (zone_type, n_nodes, connectivity) in cases {
        let name = format!("fe_roundtrip_{}", zone_type.as_tecplot_str());
        let path = temp_path(&name);
        let ds = fe_dataset(zone_type.clone(), n_nodes, connectivity.clone());
        ds.write(path.to_str().unwrap_or("")).unwrap();
        let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
        fs::remove_file(&path).ok();
        let zone = &parsed.zones[0];
        assert_eq!(zone.zone_type, zone_type);
        assert_eq!(zone.n_nodes, n_nodes);
        assert_eq!(
            zone.n_elements,
            connectivity.len() / zone_type.nodes_per_element()
        );
        assert_eq!(zone.connectivity, connectivity, "{zone_type:?}");
        assert_eq!(zone.variables[0].data.len(), n_nodes);
    }
}

#[test]
fn test_fe_mixed_ordered_and_fe_zones_roundtrip() {
    let path = temp_path("fe_mixed_zones");
    let mut ds = fe_dataset(TecplotZoneType::FEQuad, 4, vec![0, 1, 2, 3]);
    let mut ordered = TecplotZone::new_ordered("Line", 3, 1, 1);
    ordered
        .variables
        .push(TecplotVariable::new("X", vec![0.0, 1.0, 2.0]));
    ordered
        .variables
        .push(TecplotVariable::new("P", vec![5.0, 6.0, 7.0]));
    ds.zones.push(ordered);
    let second_fe = fe_dataset(TecplotZoneType::FETriangle, 3, vec![2, 1, 0]);
    ds.zones.push(second_fe.zones[0].clone());
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    fs::remove_file(&path).ok();

    assert_eq!(parsed.zones.len(), 3);
    assert_eq!(parsed.zones[0].connectivity, vec![0, 1, 2, 3]);
    assert!(parsed.zones[1].connectivity.is_empty());
    assert_eq!(parsed.zones[1].variables[1].data, vec![5.0, 6.0, 7.0]);
    assert_eq!(parsed.zones[2].connectivity, vec![2, 1, 0]);
    assert_eq!(parsed.zones[2].zone_type, TecplotZoneType::FETriangle);
}

#[test]
fn test_fe_write_rejects_connectivity_of_wrong_length() {
    let path = temp_path("fe_bad_length");
    fs::remove_file(&path).ok();
    // 2 triangles need 6 entries, only 5 supplied.
    let mut ds = fe_dataset(TecplotZoneType::FETriangle, 4, vec![0, 1, 2, 1, 3, 2]);
    ds.zones[0].connectivity.pop();
    let err = ds.write(path.to_str().unwrap_or("")).unwrap_err();
    assert!(matches!(err, IoError::General(_)), "{err:?}");
    assert!(err.to_string().contains("connectivity"), "{err}");
    // An invalid zone must not leave a half-written file behind.
    assert!(!path.exists());
}

#[test]
fn test_fe_write_rejects_missing_connectivity() {
    // The old constructor-less state: elements declared, no connectivity.
    let path = temp_path("fe_no_connectivity");
    fs::remove_file(&path).ok();
    let mut ds = fe_dataset(TecplotZoneType::FETetra, 4, vec![0, 1, 2, 3]);
    ds.zones[0].connectivity.clear();
    assert!(ds.write(path.to_str().unwrap_or("")).is_err());
    assert!(!path.exists());
}

#[test]
fn test_fe_write_rejects_out_of_range_node() {
    let path = temp_path("fe_bad_node");
    fs::remove_file(&path).ok();
    let ds = fe_dataset(TecplotZoneType::FETriangle, 4, vec![0, 1, 4]);
    let err = ds.write(path.to_str().unwrap_or("")).unwrap_err();
    assert!(matches!(err, IoError::General(_)), "{err:?}");
    assert!(!path.exists());
}

#[test]
fn test_fe_write_accepts_zero_elements() {
    let path = temp_path("fe_zero_elements");
    let ds = fe_dataset(TecplotZoneType::FETriangle, 3, Vec::new());
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    fs::remove_file(&path).ok();
    assert_eq!(parsed.zones[0].n_elements, 0);
    assert!(parsed.zones[0].connectivity.is_empty());
    assert_eq!(parsed.zones[0].variables[0].data.len(), 3);
}

/// Write `content` to a temp file, parse it and delete it again.
fn read_text(name: &str, content: &str) -> Result<TecplotDataset, IoError> {
    let path = temp_path(name);
    fs::write(&path, content).unwrap();
    let result = TecplotDataset::read(path.to_str().unwrap_or(""));
    fs::remove_file(&path).ok();
    result
}

const TRIANGLE_FILE_HEAD: &str = "TITLE = \"t\"\nVARIABLES = \"X\" \"Y\"\n";

#[test]
fn test_read_fe_connectivity_may_wrap_across_lines() {
    // Tecplot is free-format: element rows can be split anywhere.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\", N=4, E=2, ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n\
         0 0\n1 0\n0 1\n1 1\n1 2\n3\n2 4 3\n"
    );
    let ds = read_text("fe_wrapped", &text).unwrap();
    assert_eq!(ds.zones[0].connectivity, vec![0, 1, 2, 1, 3, 2]);
    assert_eq!(ds.zones[0].variables[1].data, vec![0.0, 0.0, 1.0, 1.0]);
}

#[test]
fn test_read_fe_truncated_connectivity_is_error() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\", N=4, E=2, ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n\
         0 0\n1 0\n0 1\n1 1\n1 2 3\n"
    );
    let err = read_text("fe_truncated", &text).unwrap_err();
    assert!(matches!(err, IoError::Parse(_)), "{err:?}");
    assert!(err.to_string().contains("connectivity"), "{err}");
}

#[test]
fn test_read_fe_file_without_connectivity_is_error() {
    // What the writer used to produce: node rows only.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n\
         0 0\n1 0\n0 1\n"
    );
    assert!(read_text("fe_none", &text).is_err());
}

#[test]
fn test_read_fe_node_number_above_range_is_error() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n\
         0 0\n1 0\n0 1\n1 2 4\n"
    );
    let err = read_text("fe_high", &text).unwrap_err();
    assert!(matches!(err, IoError::Parse(_)), "{err:?}");
}

#[test]
fn test_read_fe_node_number_zero_is_error() {
    // Tecplot numbers nodes from 1; 0 would become index -1.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n\
         0 0\n1 0\n0 1\n0 1 2\n"
    );
    assert!(read_text("fe_zero", &text).is_err());
}

#[test]
fn test_read_fe_non_integer_connectivity_is_error() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n\
         0 0\n1 0\n0 1\n1 2 3.5\n"
    );
    assert!(read_text("fe_float", &text).is_err());
}

#[test]
fn test_read_header_with_zonetype_before_counts() {
    // `E=` must not be matched inside `ZONETYPE=`.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE ZONETYPE=FETriangle, N=3, E=1, DATAPACKING=POINT\n\
         0 0\n1 0\n0 1\n1 2 3\n"
    );
    let ds = read_text("fe_zonetype_first", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.zone_type, TecplotZoneType::FETriangle);
    assert_eq!(zone.n_nodes, 3);
    assert_eq!(zone.n_elements, 1);
    assert_eq!(zone.connectivity, vec![0, 1, 2]);
}

#[test]
fn test_read_header_title_containing_keywords() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"N=99, E=77 ZONETYPE=FEBRICK\", N=3, E=1, \
         ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n0 0\n1 0\n0 1\n3 2 1\n"
    );
    let ds = read_text("fe_title_keywords", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.title, "N=99, E=77 ZONETYPE=FEBRICK");
    assert_eq!(zone.zone_type, TecplotZoneType::FETriangle);
    assert_eq!(zone.n_nodes, 3);
    assert_eq!(zone.n_elements, 1);
    assert_eq!(zone.connectivity, vec![2, 1, 0]);
}

#[test]
fn test_read_header_long_node_and_element_keywords() {
    // Tecplot 360 writes Nodes= / Elements= and a trailing DT=(...) list.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\" STRANDID=0, SOLUTIONTIME=0 Nodes=4, Elements=2, \
         ZONETYPE=FEQuadrilateral DATAPACKING=POINT DT=(SINGLE SINGLE)\n\
         0 0\n1 0\n1 1\n0 1\n1 2 3 4\n4 3 2 1\n"
    );
    let ds = read_text("fe_long_keys", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.zone_type, TecplotZoneType::FEQuad);
    assert_eq!(zone.n_nodes, 4);
    assert_eq!(zone.n_elements, 2);
    assert_eq!(zone.connectivity, vec![0, 1, 2, 3, 3, 2, 1, 0]);
}

#[test]
fn test_read_legacy_fepoint_header() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"old\", N=4, E=1, F=FEPOINT, ET=TETRAHEDRON\n\
         0 0\n1 0\n0 1\n1 1\n1 2 3 4\n"
    );
    let ds = read_text("fe_legacy", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.zone_type, TecplotZoneType::FETetra);
    assert_eq!(zone.connectivity, vec![0, 1, 2, 3]);
    assert_eq!(zone.variables[0].data, vec![0.0, 1.0, 0.0, 1.0]);
}

#[test]
fn test_read_legacy_fepoint_header_defaults_to_triangle() {
    let text =
        format!("{TRIANGLE_FILE_HEAD}ZONE T=\"old\", N=3, E=1, F=FEPOINT\n0 0\n1 0\n0 1\n1 2 3\n");
    let ds = read_text("fe_legacy_default", &text).unwrap();
    assert_eq!(ds.zones[0].zone_type, TecplotZoneType::FETriangle);
    assert_eq!(ds.zones[0].connectivity, vec![0, 1, 2]);
}

#[test]
fn test_read_ordered_zone_with_point_format_keyword_is_unchanged() {
    let text = format!("{TRIANGLE_FILE_HEAD}ZONE T=\"o\", I=3, J=1, F=POINT\n0 1\n2 3\n4 5\n");
    let ds = read_text("ordered_f_point", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.zone_type, TecplotZoneType::Ordered);
    assert_eq!(zone.i_dim, 3);
    assert_eq!(zone.variables[0].data, vec![0.0, 2.0, 4.0]);
    assert_eq!(zone.variables[1].data, vec![1.0, 3.0, 5.0]);
    assert!(zone.connectivity.is_empty());
}

#[test]
fn test_read_skips_keyword_lines_between_node_data_and_connectivity() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"z\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=POINT\n\
         AUXDATA Source=\"unit test\"\n0 0\n1 0\n0 1\nAUXDATA Note=\"mid\"\n1 2 3\n"
    );
    let ds = read_text("fe_auxdata", &text).unwrap();
    assert_eq!(ds.zones[0].connectivity, vec![0, 1, 2]);
    assert_eq!(ds.zones[0].variables[0].data.len(), 3);
}

#[test]
fn test_parse_header_pairs_handles_quotes_parentheses_and_spacing() {
    let pairs = TecplotReader::parse_header_pairs(
        "ZONE T = \"a, b=c\" , VARLOCATION=([3-4]=CELLCENTERED), I=2,J=3 K= 4",
    );
    let get = |k: &str| {
        pairs
            .iter()
            .find(|(key, _)| key == k)
            .map(|(_, v)| v.as_str())
    };
    assert_eq!(get("T"), Some("a, b=c"));
    assert_eq!(get("VARLOCATION"), Some("([3-4]=CELLCENTERED)"));
    assert_eq!(get("I"), Some("2"));
    assert_eq!(get("J"), Some("3"));
    assert_eq!(get("K"), Some("4"));
    assert_eq!(pairs.len(), 5);
}

// ── BLOCK data packing ──────────────────────────────────────────────────

#[test]
fn test_read_block_packed_ordered_zone() {
    // BLOCK: all X values first, then all P values.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"b\", I=3, J=1, K=1, DATAPACKING=BLOCK\n\
         0 1 2\n10 11 12\n"
    );
    let ds = read_text("block_ordered", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.variables[0].data, vec![0.0, 1.0, 2.0]);
    assert_eq!(zone.variables[1].data, vec![10.0, 11.0, 12.0]);
}

#[test]
fn test_read_block_packed_fe_zone_with_connectivity() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"b\", N=4, E=2, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK\n\
         0 1 0 1\n5 6 7 8\n1 2 3\n2 4 3\n"
    );
    let ds = read_text("block_fe", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.variables[0].data, vec![0.0, 1.0, 0.0, 1.0]);
    assert_eq!(zone.variables[1].data, vec![5.0, 6.0, 7.0, 8.0]);
    assert_eq!(zone.connectivity, vec![0, 1, 2, 1, 3, 2]);
}

#[test]
fn test_read_legacy_block_headers() {
    let ordered =
        format!("{TRIANGLE_FILE_HEAD}ZONE T=\"o\", I=2, J=2, F=BLOCK\n1 2 3 4\n5 6 7 8\n");
    let ds = read_text("legacy_block", &ordered).unwrap();
    assert_eq!(ds.zones[0].variables[0].data, vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(ds.zones[0].variables[1].data, vec![5.0, 6.0, 7.0, 8.0]);

    let fe = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"f\", N=4, E=1, F=FEBLOCK, ET=QUADRILATERAL\n\
         0 1 1 0\n0 0 1 1\n1 2 3 4\n"
    );
    let ds = read_text("legacy_feblock", &fe).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.zone_type, TecplotZoneType::FEQuad);
    assert_eq!(zone.variables[0].data, vec![0.0, 1.0, 1.0, 0.0]);
    assert_eq!(zone.variables[1].data, vec![0.0, 0.0, 1.0, 1.0]);
    assert_eq!(zone.connectivity, vec![0, 1, 2, 3]);
}

#[test]
fn test_read_point_packing_keyword_is_explicit_or_default() {
    let explicit =
        format!("{TRIANGLE_FILE_HEAD}ZONE T=\"p\", I=2, J=1, K=1, DATAPACKING=POINT\n1 5\n2 6\n");
    let ds = read_text("point_explicit", &explicit).unwrap();
    assert_eq!(ds.zones[0].variables[0].data, vec![1.0, 2.0]);
    assert_eq!(ds.zones[0].variables[1].data, vec![5.0, 6.0]);
    let implicit = format!("{TRIANGLE_FILE_HEAD}ZONE T=\"p\", I=2, J=1, K=1\n1 5\n2 6\n");
    let ds = read_text("point_default", &implicit).unwrap();
    assert_eq!(ds.zones[0].variables[0].data, vec![1.0, 2.0]);
}

#[test]
fn test_read_cell_centered_variables_are_rejected() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"c\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK, \
         VARLOCATION=([2]=CELLCENTERED)\n0 1 0\n7\n1 2 3\n"
    );
    let err = read_text("cell_centered", &text).unwrap_err();
    assert!(matches!(err, IoError::Parse(_)), "{err:?}");
    assert!(err.to_string().contains("VARLOCATION"), "{err}");
}

#[test]
fn test_read_nodal_varlocation_is_accepted() {
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"n\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK, \
         VARLOCATION=([1-2]=NODAL)\n0 1 0\n5 6 7\n1 2 3\n"
    );
    let ds = read_text("nodal_varlocation", &text).unwrap();
    assert_eq!(ds.zones[0].variables[1].data, vec![5.0, 6.0, 7.0]);
    assert_eq!(ds.zones[0].connectivity, vec![0, 1, 2]);
}

// ── hostile headers ─────────────────────────────────────────────────────

#[test]
fn test_read_absurd_ordered_dimensions_do_not_overflow() {
    // 2^32 * 2^32 * 2^32 does not fit in a usize: the header must be
    // accepted (the data is simply shorter than declared) without an
    // arithmetic-overflow panic.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"huge\", I=4294967296, J=4294967296, K=4294967296\n\
         1 2\n3 4\n"
    );
    let ds = read_text("huge_ordered", &text).unwrap();
    let zone = &ds.zones[0];
    assert_eq!(zone.zone_type, TecplotZoneType::Ordered);
    assert_eq!(zone.point_count(), usize::MAX);
    assert_eq!(zone.variables[0].data, vec![1.0, 3.0]);
    assert_eq!(zone.variables[1].data, vec![2.0, 4.0]);
}

#[test]
fn test_read_absurd_element_count_is_a_parse_error_not_an_allocation() {
    // E claims far more connectivity than the file holds.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"huge\", N=3, E=18446744073709551615, \
         ZONETYPE=FEBRICK\n0 0\n1 0\n0 1\n1 2 3 4 5 6 7 8\n"
    );
    let err = read_text("huge_elements", &text).unwrap_err();
    assert!(matches!(err, IoError::Parse(_)), "{err:?}");
}

#[test]
fn test_new_ordered_with_absurd_dimensions_saturates() {
    let zone = TecplotZone::new_ordered("big", usize::MAX, 2, 3);
    assert_eq!(zone.n_nodes, usize::MAX);
    assert_eq!(zone.point_count(), usize::MAX);
}

// ── writer: one value per point, exact layout ───────────────────────────

#[test]
fn test_write_rows_use_fifteen_digit_scientific_notation() {
    let path = temp_path("row_format");
    let mut ds = TecplotDataset::new("fmt");
    ds.variables = vec!["X".to_string(), "P".to_string()];
    let mut zone = TecplotZone::new_ordered("Z", 2, 1, 1);
    zone.variables
        .push(TecplotVariable::new("X", vec![1.0, -2.5]));
    zone.variables
        .push(TecplotVariable::new("P", vec![0.0, 1.0e-3]));
    ds.zones.push(zone);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let content = fs::read_to_string(&path).unwrap();
    fs::remove_file(&path).ok();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(
        lines,
        vec![
            "TITLE = \"fmt\"",
            "VARIABLES = \"X\", \"P\"",
            "ZONE T=\"Z\", I=2, J=1, K=1, DATAPACKING=POINT",
            "1.000000000000000e0 0.000000000000000e0",
            "-2.500000000000000e0 1.000000000000000e-3",
        ]
    );
}

#[test]
fn test_write_non_finite_values_round_trip() {
    let path = temp_path("non_finite");
    let mut ds = TecplotDataset::new("nf");
    ds.variables = vec!["X".to_string()];
    let mut zone = TecplotZone::new_ordered("Z", 3, 1, 1);
    zone.variables.push(TecplotVariable::new(
        "X",
        vec![f64::NAN, f64::INFINITY, f64::NEG_INFINITY],
    ));
    ds.zones.push(zone);
    ds.write(path.to_str().unwrap_or("")).unwrap();
    let parsed = TecplotDataset::read(path.to_str().unwrap_or("")).unwrap();
    fs::remove_file(&path).ok();
    let data = &parsed.zones[0].variables[0].data;
    assert!(data[0].is_nan());
    assert_eq!(data[1], f64::INFINITY);
    assert_eq!(data[2], f64::NEG_INFINITY);
}

#[test]
fn test_write_rejects_variable_shorter_than_point_count() {
    let path = temp_path("short_variable");
    fs::remove_file(&path).ok();
    let mut ds = make_1d_dataset("", 4);
    ds.zones[0].variables[1].data.pop();
    let err = ds.write(path.to_str().unwrap_or("")).unwrap_err();
    assert!(matches!(err, IoError::General(_)), "{err:?}");
    let message = err.to_string();
    assert!(message.contains("\"P\""), "{message}");
    assert!(message.contains("3 values"), "{message}");
    assert!(message.contains("4 points"), "{message}");
    // Missing values are never padded with zeros, and nothing is written.
    assert!(!path.exists());
}

#[test]
fn test_write_rejects_variable_longer_than_point_count() {
    let path = temp_path("long_variable");
    fs::remove_file(&path).ok();
    let mut ds = make_1d_dataset("", 3);
    ds.zones[0].variables[0].data.push(9.0);
    let err = ds.write(path.to_str().unwrap_or("")).unwrap_err();
    assert!(matches!(err, IoError::General(_)), "{err:?}");
    assert!(!path.exists());
}

#[test]
fn test_write_rejects_fe_zone_with_wrong_variable_length() {
    let path = temp_path("fe_short_variable");
    fs::remove_file(&path).ok();
    let mut ds = fe_dataset(TecplotZoneType::FETriangle, 4, vec![0, 1, 2, 1, 3, 2]);
    ds.zones[0].variables[0].data.truncate(2);
    let err = ds.write(path.to_str().unwrap_or("")).unwrap_err();
    assert!(matches!(err, IoError::General(_)), "{err:?}");
    assert!(err.to_string().contains("\"X\""), "{err}");
    assert!(!path.exists());
}

#[test]
fn test_write_does_not_expand_absurd_declared_sizes() {
    // Zones whose declared size is far larger than their data must be
    // rejected at once instead of being padded row by row until memory
    // runs out.
    let path = temp_path("absurd_size");
    fs::remove_file(&path).ok();
    let mut ds = TecplotDataset::new("absurd");
    ds.variables = vec!["X".to_string()];
    let mut ordered = TecplotZone::new_ordered("big", usize::MAX, 2, 3);
    ordered
        .variables
        .push(TecplotVariable::new("X", vec![1.0, 2.0]));
    ds.zones.push(ordered);
    assert!(ds.write(path.to_str().unwrap_or("")).is_err());
    assert!(!path.exists());

    let mut fe = TecplotZone::new_fe("fe", TecplotZoneType::FETriangle, usize::MAX, 0, Vec::new());
    fe.variables.push(TecplotVariable::new("X", vec![1.0]));
    ds.zones = vec![fe];
    assert!(ds.write(path.to_str().unwrap_or("")).is_err());
    assert!(!path.exists());
}

#[test]
fn test_read_hostile_header_then_write_fails_fast() {
    // The reader accepts a header that declares 2^96 points over two data
    // rows (see `test_read_absurd_ordered_dimensions_do_not_overflow`);
    // writing that zone back must report the inconsistency.
    let text = format!(
        "{TRIANGLE_FILE_HEAD}ZONE T=\"huge\", I=4294967296, J=4294967296, K=4294967296\n\
         1 2\n3 4\n"
    );
    let ds = read_text("hostile_roundtrip_in", &text).unwrap();
    let path = temp_path("hostile_roundtrip_out");
    fs::remove_file(&path).ok();
    let err = ds.write(path.to_str().unwrap_or("")).unwrap_err();
    assert!(matches!(err, IoError::General(_)), "{err:?}");
    assert!(!path.exists());
}

#[test]
fn test_read_keyword_detection_ignores_case_and_leading_space() {
    let text = "title = \"lower\"\n  Variables = \"X\"\n\
                zone t=\"z\", i=2, j=1, k=1\n  1.5\n2.5\n";
    let ds = read_text("keyword_case", text).unwrap();
    assert_eq!(ds.title, "lower");
    assert_eq!(ds.variables, vec!["X"]);
    assert_eq!(ds.zones.len(), 1);
    assert_eq!(ds.zones[0].title, "z");
    assert_eq!(ds.zones[0].variables[0].data, vec![1.5, 2.5]);
}
