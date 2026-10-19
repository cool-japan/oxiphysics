// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The Fluent reader against the layout Fluent itself writes: the data of a
//! section starts on the header line, header numbers are hexadecimal, faces of
//! a fixed-type zone carry no type prefix and two-dimensional meshes have two
//! coordinates per node.  The meshes below are typed out in that layout instead
//! of being produced by [`FluentWriter`], so the reader is checked against the
//! format and not only against its own writer.

use super::*;

fn parse(text: &str) -> crate::Result<FluentMesh> {
    FluentReader::parse(text.as_bytes())
}

fn parsed(text: &str) -> FluentMesh {
    parse(text).expect("the mesh should parse")
}

fn error_text(text: &str) -> String {
    parse(text)
        .expect_err("the mesh should be rejected")
        .to_string()
}

/// A square cut into two triangles, as a 2-D mesh from Fluent looks like.
const SQUARE_2D: &str = "\
(0 \"GAMBIT to Fluent File\")
(0 \"\")
(1 \"2-D mesh of a square (two triangles)\")
(0 \"Dimension:\")
(2 2)
(0 \"Node Section:\")
(10 (0 1 4 0 2))
(10 (1 1 4 1 2)(
0.00000000e+00 0.00000000e+00
1.00000000e+00 0.00000000e+00
1.00000000e+00 1.00000000e+00
0.00000000e+00 1.00000000e+00
))
(0 \"Cell Section (cells):\")
(12 (0 1 2 0))
(12 (3 1 2 1 1))
(0 \"Face Section:\")
(13 (0 1 5 0))
(13 (4 1 1 2 2)(
1 3 1 2
))
(13 (5 2 5 3 2)(
1 2 1 0
2 3 1 0
3 4 2 0
4 1 2 0
))
(0 \"Zone Sections\")
(45 (3 fluid fluid-3)())
(45 (4 interior interior-4)())
(45 (5 wall wall-5)())
";

/// Two tetrahedra glued along a face: two node zones (the second one starting
/// at node 4), an interior face zone of fixed type and a boundary face zone of
/// mixed type.
const TETS_3D: &str = "\
(0 \"Dimension:\")
(2 3)
(10 (0 1 5 0 3))
(10 (1 1 3 2 3)(
0 0 0
1 0 0
0 1 0
))
(10 (2 4 5 1 3)(
0 0 1
1.5e+00 -1.25e-01 1
))
(12 (0 1 2 0))
(12 (1 1 2 1 2))
(13 (0 1 7 0))
(13 (2 1 1 2 3)(
2 3 4 1 2
))
(13 (3 2 7 3 0)(
3 1 2 3 1 0
3 1 2 4 1 0
3 1 3 4 1 0
3 2 3 5 2 0
3 2 4 5 2 0
3 3 4 5 2 0
))
(45 (1 fluid fluid-1)())
(45 (2 interior interior-2)())
(45 (3 wall wall-3)())
";

#[test]
fn two_dimensional_mesh_is_read_with_two_coordinates_per_node() {
    let mesh = parsed(SQUARE_2D);
    assert_eq!(mesh.node_count(), 4);
    let coordinates: Vec<[f64; 3]> = mesh.nodes.iter().map(|n| n.coordinates).collect();
    assert_eq!(
        coordinates,
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0]
        ]
    );
    let ids: Vec<usize> = mesh.nodes.iter().map(|n| n.id).collect();
    assert_eq!(ids, vec![1, 2, 3, 4]);
}

#[test]
fn two_dimensional_cells_faces_and_zones() {
    let mesh = parsed(SQUARE_2D);
    assert_eq!(
        mesh.cells,
        vec![
            FluentCell::new(1, FluentCellType::Tri, 3),
            FluentCell::new(2, FluentCellType::Tri, 3),
        ]
    );
    assert_eq!(mesh.face_count(), 5);
    assert_eq!(
        mesh.faces[0],
        FluentFace::new(1, FluentFaceType::Line, vec![1, 3], 1, 2)
    );
    assert!(!mesh.faces[0].is_boundary());
    assert_eq!(
        mesh.faces[3],
        FluentFace::new(4, FluentFaceType::Line, vec![3, 4], 2, 0)
    );
    assert!(mesh.faces[1..].iter().all(FluentFace::is_boundary));
    let face_ids: Vec<usize> = mesh.faces.iter().map(|f| f.id).collect();
    assert_eq!(face_ids, vec![1, 2, 3, 4, 5]);
    assert_eq!(
        mesh.zones,
        vec![
            (3, FluentZoneType::Fluid),
            (4, FluentZoneType::Interior),
            (5, FluentZoneType::Wall)
        ]
    );
}

#[test]
fn three_dimensional_mesh_with_two_node_zones() {
    let mesh = parsed(TETS_3D);
    let ids: Vec<usize> = mesh.nodes.iter().map(|n| n.id).collect();
    assert_eq!(ids, vec![1, 2, 3, 4, 5]);
    assert_eq!(mesh.nodes[2].coordinates, [0.0, 1.0, 0.0]);
    assert_eq!(mesh.nodes[3].coordinates, [0.0, 0.0, 1.0]);
    assert_eq!(mesh.nodes[4].coordinates, [1.5, -0.125, 1.0]);
}

#[test]
fn fixed_type_and_mixed_face_zones_in_one_mesh() {
    let mesh = parsed(TETS_3D);
    assert_eq!(mesh.face_count(), 7);
    let ids: Vec<usize> = mesh.faces.iter().map(|f| f.id).collect();
    assert_eq!(ids, vec![1, 2, 3, 4, 5, 6, 7]);
    // Fixed-type zone: no prefix, three nodes, the two cells at the end.
    assert_eq!(
        mesh.faces[0],
        FluentFace::new(1, FluentFaceType::Triangle, vec![2, 3, 4], 1, 2)
    );
    // Mixed zone: every row starts with its type code.
    assert_eq!(
        mesh.faces[1],
        FluentFace::new(2, FluentFaceType::Triangle, vec![1, 2, 3], 1, 0)
    );
    assert_eq!(
        mesh.faces[6],
        FluentFace::new(7, FluentFaceType::Triangle, vec![3, 4, 5], 2, 0)
    );
    assert_eq!(mesh.cell_count(), 2);
    assert!(
        mesh.cells
            .iter()
            .all(|c| c.cell_type == FluentCellType::Tet)
    );
    assert_eq!(mesh.zone_type(1), Some(&FluentZoneType::Fluid));
    assert_eq!(mesh.zone_type(3), Some(&FluentZoneType::Wall));
}

#[test]
fn mixed_and_polygonal_zones_read_lines_triangles_and_quads() {
    let text = "\
(13 (1 1 3 3 0)(
2 1 2 1 0
3 1 2 3 1 0
4 1 2 3 4 2 0
))
(13 (2 4 5 3 5)(
4 5 6 7 8 1 0
3 5 6 7 1 0
))
";
    let mesh = parsed(text);
    let types: Vec<FluentFaceType> = mesh.faces.iter().map(|f| f.face_type.clone()).collect();
    assert_eq!(
        types,
        vec![
            FluentFaceType::Line,
            FluentFaceType::Triangle,
            FluentFaceType::Quad,
            FluentFaceType::Quad,
            FluentFaceType::Triangle
        ]
    );
    assert_eq!(mesh.faces[2].node_ids, vec![1, 2, 3, 4]);
    assert_eq!(mesh.faces[2].left_cell, 2);
    assert_eq!(mesh.faces[3].id, 4);
    assert_eq!(mesh.faces[3].node_ids, vec![5, 6, 7, 8]);
    assert_eq!(mesh.faces[4].node_ids, vec![5, 6, 7]);
}

#[test]
fn quadrilateral_zone_reads_four_nodes_without_prefix() {
    let mesh = parsed("(13 (1 1 2 3 4)(\n1 2 3 4 1 0\n4 3 2 1 2 1\n))\n");
    assert_eq!(mesh.faces.len(), 2);
    assert_eq!(mesh.faces[0].face_type, FluentFaceType::Quad);
    assert_eq!(mesh.faces[0].node_ids, vec![1, 2, 3, 4]);
    assert_eq!(mesh.faces[1].node_ids, vec![4, 3, 2, 1]);
    assert_eq!((mesh.faces[1].left_cell, mesh.faces[1].right_cell), (2, 1));
}

#[test]
fn header_numbers_are_hexadecimal() {
    let text = "\
(10 (A a C 1 3)(
0 0 0
1 0 0
2 0 0
))
(13 (b 8 a 3 3)(
a b c 1 0
b C a 2 1
1F 20 21 0 3
))
(12 (c 1f 21 1 4))
";
    let mesh = parsed(text);
    let node_ids: Vec<usize> = mesh.nodes.iter().map(|n| n.id).collect();
    assert_eq!(node_ids, vec![10, 11, 12]);
    let face_ids: Vec<usize> = mesh.faces.iter().map(|f| f.id).collect();
    assert_eq!(face_ids, vec![8, 9, 10]);
    assert_eq!(mesh.faces[0].node_ids, vec![10, 11, 12]);
    assert_eq!(mesh.faces[1].node_ids, vec![11, 12, 10]);
    assert_eq!(mesh.faces[2].node_ids, vec![31, 32, 33]);
    assert_eq!((mesh.faces[2].left_cell, mesh.faces[2].right_cell), (0, 3));
    assert_eq!(mesh.cell_count(), 3);
    assert_eq!(mesh.cells[0].id, 31);
    assert_eq!(mesh.cells[2].id, 33);
    assert_eq!(mesh.cells[0].zone_id, 12);
    assert_eq!(mesh.cells[0].cell_type, FluentCellType::Hex);
}

#[test]
fn line_breaks_inside_the_data_do_not_matter() {
    let text = "(2 3)\r\n(10 (1 1 2 1 3)(\r\n0 0\r\n0 1 0\r\n0\r\n))\r\n(13 (1 1 1 2 3)( 1 2\r\n3 1\r\n0 ))\r\n";
    let mesh = parsed(text);
    assert_eq!(mesh.nodes[0].coordinates, [0.0, 0.0, 0.0]);
    assert_eq!(mesh.nodes[1].coordinates, [1.0, 0.0, 0.0]);
    assert_eq!(mesh.faces[0].node_ids, vec![1, 2, 3]);
    assert_eq!(mesh.faces[0].right_cell, 0);
}

#[test]
fn data_may_start_on_the_header_line_or_on_the_next() {
    let same_line = "(10 (1 1 2 1 3)(\n0 0 0\n1 0 0\n))\n";
    let own_line = "(10 (1 1 2 1 3)\n(\n0 0 0\n1 0 0\n))\n";
    let squeezed = "(10 (1 1 2 1 3)(0 0 0 1 0 0))";
    let a = parsed(same_line);
    assert_eq!(a.nodes.len(), 2);
    assert_eq!(parsed(own_line).nodes, a.nodes);
    assert_eq!(parsed(squeezed).nodes, a.nodes);
}

#[test]
fn node_dimension_falls_back_to_the_dimension_section_then_to_the_data() {
    // No nd in the header: the (2 2) section decides.
    let flat = parsed("(2 2)\n(10 (1 1 3 1)(\n0 0\n1 0\n0 1\n))\n");
    assert_eq!(flat.nodes[2].coordinates, [0.0, 1.0, 0.0]);
    // Neither: the amount of data decides when it divides evenly.
    let two = parsed("(10 (1 1 3 1)(0 0 1 0 0 1))");
    assert_eq!(two.nodes[1].coordinates, [1.0, 0.0, 0.0]);
    let three = parsed("(10 (1 1 2 1)(0 0 7 1 0 8))");
    assert_eq!(three.nodes[1].coordinates, [1.0, 0.0, 8.0]);
    let message = error_text("(10 (1 1 3 1)(0 0 1 0 0 1 5))");
    assert!(message.contains("neither 2 nor 3"), "{message}");
    // The header wins over the dimension section.
    let header = parsed("(2 2)\n(10 (1 1 1 1 3)(4 5 6))");
    assert_eq!(header.nodes[0].coordinates, [4.0, 5.0, 6.0]);
    let message = error_text("(10 (1 1 1 1 4)(4 5 6 7))");
    assert!(message.contains("4 coordinates per node"), "{message}");
}

#[test]
fn declaration_sections_carry_no_data() {
    let mesh = parsed("(10 (0 1 ff 0 3))\n(12 (0 1 ff 0))\n(13 (0 1 1ff 0))\n(12 (0 1 4 0 0))\n");
    assert!(mesh.nodes.is_empty() && mesh.faces.is_empty() && mesh.cells.is_empty());
}

#[test]
fn malformed_coordinates_name_the_line_and_the_token() {
    let message = error_text("(2 3)\n(10 (1 1 2 1 3)(\n0 0 0\n1 abc 0\n))\n");
    assert!(
        message.contains("line 4") && message.contains("\"abc\""),
        "{message}"
    );
    let non_finite = error_text("(10 (1 1 1 1 3)(0 NaN 0))");
    assert!(non_finite.contains("finite"), "{non_finite}");
    let infinite = error_text("(10 (1 1 1 1 3)(0 inf 0))");
    assert!(infinite.contains("finite"), "{infinite}");
}

#[test]
fn node_count_must_match_the_header() {
    let few = error_text("(10 (1 1 3 1 3)(\n0 0 0\n1 0 0\n))\n");
    assert!(
        few.contains("declares 3 nodes but its data ends after 2"),
        "{few}"
    );
    let short_row = error_text("(10 (1 1 2 1 3)(\n0 0 0\n1 0\n))\n");
    assert!(short_row.contains("ends after 1"), "{short_row}");
    let many = error_text("(10 (1 1 1 1 3)(\n0 0 0\n1 0 0\n))\n");
    assert!(many.contains("more coordinates than the 1 nodes"), "{many}");
    let no_data = error_text("(10 (1 1 2 1 3))\n");
    assert!(no_data.contains("no data block"), "{no_data}");
}

#[test]
fn face_data_must_match_the_header() {
    let truncated = error_text("(13 (1 1 2 2 3)(\n1 2 3 1 0\n))\n");
    assert!(truncated.contains("ends before face 2"), "{truncated}");
    let surplus = error_text("(13 (1 1 1 2 3)(\n1 2 3 1 0\n4 5 6 1 0\n))\n");
    assert!(surplus.contains("more data than the 1 faces"), "{surplus}");
    let cut_row = error_text("(13 (1 1 1 2 3)(\n1 2 3 1\n))\n");
    assert!(cut_row.contains("ends before face 1"), "{cut_row}");
    let bad_number = error_text("(13 (1 1 1 2 3)(\n1 2 x 1 0\n))\n");
    assert!(
        bad_number.contains("hexadecimal") && bad_number.contains("line 2"),
        "{bad_number}"
    );
}

#[test]
fn faces_that_cannot_be_represented_are_errors_not_truncated() {
    let pentagon = error_text("(13 (1 1 1 3 5)(\n5 1 2 3 4 5 1 0\n))\n");
    assert!(
        pentagon.contains("5 nodes") && pentagon.contains("can be represented"),
        "{pentagon}"
    );
    let point = error_text("(13 (1 1 1 3 0)(\n1 7 1 0\n))\n");
    assert!(point.contains("1 nodes"), "{point}");
    let unsupported = error_text("(13 (1 1 1 3 6)(\n1 2 3 4 5 6 1 0\n))\n");
    assert!(unsupported.contains("face type 6"), "{unsupported}");
    let missing = error_text("(13 (1 1 1 3))\n");
    assert!(missing.contains("no face type"), "{missing}");
}

#[test]
fn mixed_cell_zone_lists_the_type_of_every_cell() {
    let mesh = parsed("(12 (7 1 4 1 0)(\n2 5 6 1\n))\n");
    let types: Vec<FluentCellType> = mesh.cells.iter().map(|c| c.cell_type.clone()).collect();
    assert_eq!(
        types,
        vec![
            FluentCellType::Tet,
            FluentCellType::Pyramid,
            FluentCellType::Wedge,
            FluentCellType::Tri
        ]
    );
    let ids: Vec<usize> = mesh.cells.iter().map(|c| c.id).collect();
    assert_eq!(ids, vec![1, 2, 3, 4]);
    assert!(mesh.cells.iter().all(|c| c.zone_id == 7));
    let short = error_text("(12 (7 1 4 1 0)(2 5 6))");
    assert!(short.contains("lists the types of only 3 cells"), "{short}");
    let long = error_text("(12 (7 1 3 1 0)(2 5 6 1))");
    assert!(long.contains("more cell types"), "{long}");
    let unlisted = error_text("(12 (7 1 3 1 0))");
    assert!(unlisted.contains("no data block"), "{unlisted}");
}

#[test]
fn uniform_cell_zones_and_unsupported_cell_types() {
    let mesh = parsed("(12 (2 1 a 1 4))\n(12 (3 b c 1 5))\n(12 (4 d d 1 6))\n");
    assert_eq!(mesh.cell_count(), 13);
    assert_eq!(mesh.cells[0].cell_type, FluentCellType::Hex);
    assert_eq!(mesh.cells[9].id, 10);
    assert_eq!(mesh.cells[10].cell_type, FluentCellType::Pyramid);
    assert_eq!(mesh.cells[11].zone_id, 3);
    assert_eq!(mesh.cells[12].id, 13);
    assert_eq!(mesh.cells[12].cell_type, FluentCellType::Wedge);
    assert_eq!(mesh.cells[12].zone_id, 4);
    let polyhedral = error_text("(12 (2 1 4 1 7))");
    assert!(polyhedral.contains("polyhedral"), "{polyhedral}");
    let unknown = error_text("(12 (2 1 4 1 9))");
    assert!(
        unknown.contains("not a known cell element type"),
        "{unknown}"
    );
    let no_type = error_text("(12 (2 1 4 1))");
    assert!(no_type.contains("no element type"), "{no_type}");
}

#[test]
fn absurd_index_ranges_are_rejected_without_allocating() {
    let cells = error_text("(12 (1 1 ffffffff 1 2))");
    assert!(cells.contains("limit"), "{cells}");
    let faces = error_text("(13 (1 1 ffffffff 2 2)(\n1 2 1 0\n))");
    assert!(faces.contains("ends before face 2"), "{faces}");
    let nodes = error_text("(10 (1 1 ffffffff 1 3)(\n0 0 0\n))");
    assert!(nodes.contains("data ends after 1"), "{nodes}");
}

#[test]
fn index_ranges_must_start_at_one_and_not_run_backwards() {
    for header in ["(10 (1 0 2 1 3)", "(10 (1 3 2 1 3)"] {
        let message = error_text(&format!("{header}(0 0 0 0 0 0))"));
        assert!(message.contains("indices"), "{header}: {message}");
    }
}

#[test]
fn header_problems_are_reported() {
    let bad_hex = error_text("(10 (1 1 zz 1 3)(0 0 0))");
    assert!(bad_hex.contains("hexadecimal"), "{bad_hex}");
    let too_short = error_text("(10 (1 1))");
    assert!(too_short.contains("2 fields"), "{too_short}");
    let too_long = error_text("(13 (1 1 1 2 2 9)(\n1 2 1 0\n))");
    assert!(too_long.contains("6 fields"), "{too_long}");
    let no_header = error_text("(10 3)");
    assert!(no_header.contains("no (zone first last"), "{no_header}");
    let bad_dimension = error_text("(2 4)");
    assert!(bad_dimension.contains("2 or 3"), "{bad_dimension}");
}

#[test]
fn binary_sections_are_rejected_with_a_clean_parse_error() {
    let mut bytes = b"(2 3)\n(3010 (1 1 4 1 3)(".to_vec();
    bytes.extend_from_slice(&[0xff, 0xfe, b'(', b')', 0x00, 0x80, b')']);
    bytes.extend_from_slice(b"End of Binary Section  3010)\n");
    let error = FluentReader::parse(bytes.as_slice()).unwrap_err();
    assert!(matches!(error, IoError::Parse(_)), "{error:?}");
    let message = error.to_string();
    assert!(
        message.contains("3010") && message.contains("binary") && message.contains("line 2"),
        "{message}"
    );
    for number in [2010, 2012, 2013, 3012, 3013, 2018, 3058] {
        let message = error_text(&format!("({number} (1 1 4 1 3)(\nbinary"));
        assert!(
            message.contains(&number.to_string()) && message.contains("binary"),
            "{number}: {message}"
        );
    }
}

#[test]
fn unknown_and_comment_sections_are_skipped() {
    let text = "\
(0 \"Comment with a ) paren and a ( too\")
(1 \"Title\")
(37 (1 2)(foo (bar) \"baz)\"))
(18 (b 7 9 3 4)(1 2 3))
(59 (1 5 2 3)(4 5 6 7))
(10 (1 1 1 1 3)(1 2 3))
(2 3)
";
    let mesh = parsed(text);
    assert_eq!(mesh.nodes.len(), 1);
    assert_eq!(mesh.nodes[0].coordinates, [1.0, 2.0, 3.0]);
}

#[test]
fn zone_sections_in_both_numberings() {
    let text = "\
(39 (4 wall wall-4 1)(
(condition . 1)
(Fx . ((a . 1) (b . 2)))
))
(45 (10 velocity-inlet vel-in)())
(45 (c pressure-outlet outlet-c)())
(45 (4 symmetry duplicate-of-4)())
(45 (11 solid solid-11)(
))
";
    let mesh = parsed(text);
    assert_eq!(
        mesh.zones,
        vec![
            (4, FluentZoneType::Wall),
            (10, FluentZoneType::Inlet),
            (12, FluentZoneType::Outlet),
            (11, FluentZoneType::Solid)
        ]
    );
    let bad = error_text("(45 (x fluid f)())");
    assert!(bad.contains("not a zone id"), "{bad}");
    let short = error_text("(45 (3))");
    assert!(short.contains("zone id and a zone type"), "{short}");
    let headerless = error_text("(45 3 fluid)");
    assert!(headerless.contains("no (zone-id"), "{headerless}");
}

#[test]
fn structural_errors_name_the_line() {
    let stray = error_text("(2 3)\nhello\n");
    assert!(
        stray.contains("line 2") && stray.contains("expected '('"),
        "{stray}"
    );
    let unterminated = error_text("(2 3)\n(10 (1 1 2 1 3)(\n0 0 0\n");
    assert!(unterminated.contains("unterminated"), "{unterminated}");
    let no_number = error_text("(\"comment\")");
    assert!(no_number.contains("its number"), "{no_number}");
    let trailing = error_text("(2 3 4)");
    assert!(trailing.contains("should end"), "{trailing}");
    let junk_after_data = error_text("(10 (1 1 1 1 3)(0 0 0) junk)");
    assert!(junk_after_data.contains("should end"), "{junk_after_data}");
}

#[test]
fn empty_input_and_byte_order_mark_are_fine() {
    assert_eq!(parsed("").node_count(), 0);
    assert_eq!(parsed("  \n\t\n").face_count(), 0);
    let mut bytes = vec![0xEF, 0xBB, 0xBF];
    bytes.extend_from_slice(b"(10 (1 1 1 1 3)(1 2 3))");
    let mesh = FluentReader::parse(bytes.as_slice()).unwrap();
    assert_eq!(mesh.nodes[0].coordinates, [1.0, 2.0, 3.0]);
}

#[test]
fn files_are_read_from_disk() {
    let path = std::env::temp_dir().join("oxiphysics_fluent_real_layout.msh");
    std::fs::write(&path, SQUARE_2D).unwrap();
    let loaded = FluentMesh::read(path.to_str().unwrap()).unwrap();
    assert_eq!(loaded.node_count(), 4);
    assert_eq!(loaded.face_count(), 5);
    assert_eq!(loaded.cell_count(), 2);
    let _ = std::fs::remove_file(&path);
}

#[test]
fn writer_output_reads_back_with_every_field() {
    let mut mesh = FluentMesh::new();
    mesh.nodes.push(FluentNode::new(1, [0.0, 0.0, 0.0]));
    mesh.nodes.push(FluentNode::new(2, [1.0, 0.0, 0.5]));
    mesh.nodes.push(FluentNode::new(3, [1.0, 1.0, -0.25]));
    mesh.nodes.push(FluentNode::new(4, [0.0, 1.0, 2.0e-3]));
    mesh.cells.push(FluentCell::new(1, FluentCellType::Hex, 1));
    mesh.cells.push(FluentCell::new(2, FluentCellType::Hex, 1));
    mesh.cells
        .push(FluentCell::new(3, FluentCellType::Wedge, 2));
    mesh.faces.push(FluentFace::new(
        1,
        FluentFaceType::Quad,
        vec![1, 2, 3, 4],
        1,
        2,
    ));
    mesh.faces.push(FluentFace::new(
        2,
        FluentFaceType::Triangle,
        vec![1, 2, 3],
        1,
        0,
    ));
    mesh.faces
        .push(FluentFace::new(3, FluentFaceType::Line, vec![3, 4], 0, 3));
    mesh.add_zone(1, FluentZoneType::Fluid);
    mesh.add_zone(2, FluentZoneType::Solid);
    let text = FluentWriter::new(&mesh).to_string().unwrap();
    let back = parsed(&text);
    assert_eq!(back.nodes.len(), 4);
    for (written, read) in mesh.nodes.iter().zip(&back.nodes) {
        assert_eq!(written.id, read.id);
        for axis in 0..3 {
            assert!((written.coordinates[axis] - read.coordinates[axis]).abs() < 1e-9);
        }
    }
    assert_eq!(back.cells, mesh.cells);
    assert_eq!(back.faces, mesh.faces);
    assert_eq!(back.zones, mesh.zones);
}

#[test]
fn writer_opens_data_groups_on_the_header_line_like_fluent() {
    let mut mesh = FluentMesh::new();
    for (id, coordinates) in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    .into_iter()
    .enumerate()
    {
        mesh.nodes.push(FluentNode::new(id + 1, coordinates));
    }
    mesh.faces.push(FluentFace::new(
        1,
        FluentFaceType::Triangle,
        vec![1, 2, 3],
        1,
        0,
    ));
    mesh.faces.push(FluentFace::new(
        2,
        FluentFaceType::Triangle,
        vec![1, 2, 4],
        1,
        0,
    ));
    let text = FluentWriter::new(&mesh).to_string().unwrap();
    assert!(text.contains("(2 3)\n"), "{text}");
    assert!(text.contains("(10 (1 1 4 1 3)(\n"), "{text}");
    assert!(text.contains("(13 (1 1 2 2 0)(\n"), "{text}");
    assert!(
        !text.lines().any(|line| line.trim() == "("),
        "no data group may open on a line of its own:\n{text}"
    );
    let back = parsed(&text);
    assert_eq!(back.node_count(), 4);
    assert_eq!(back.faces, mesh.faces);
}

#[test]
fn writer_gives_a_planar_mesh_two_coordinates_per_node() {
    let mut mesh = FluentMesh::new();
    mesh.nodes.push(FluentNode::new(1, [0.0, 0.0, 0.0]));
    mesh.nodes.push(FluentNode::new(2, [1.5, -2.0, 0.0]));
    mesh.nodes.push(FluentNode::new(3, [0.25, 4.0, 0.0]));
    let text = FluentWriter::new(&mesh).to_string().unwrap();
    assert!(text.contains("(2 2)\n"), "{text}");
    assert!(text.contains("(10 (1 1 3 1 2)(\n"), "{text}");
    let rows: Vec<&str> = text
        .lines()
        .skip_while(|line| !line.starts_with("(10 (1 "))
        .skip(1)
        .take(3)
        .collect();
    assert_eq!(rows.len(), 3);
    for row in &rows {
        assert_eq!(row.split_whitespace().count(), 2, "{row}");
    }
    let back = parsed(&text);
    assert_eq!(back.nodes.len(), 3);
    assert_eq!(back.nodes[1].coordinates, [1.5, -2.0, 0.0]);
    assert_eq!(back.nodes[2].coordinates, [0.25, 4.0, 0.0]);
}

#[test]
fn cell_and_face_type_helpers() {
    for code in 1..=6u32 {
        let cell_type = FluentCellType::try_from_code(code).unwrap();
        assert_eq!(cell_type.code(), code);
        assert_eq!(FluentCellType::from_code(code), cell_type);
    }
    assert_eq!(FluentCellType::try_from_code(0), None);
    assert_eq!(FluentCellType::try_from_code(7), None);
    assert_eq!(FluentCellType::from_code(7), FluentCellType::Tet);
    assert_eq!(
        FluentFaceType::from_node_count(2),
        Some(FluentFaceType::Line)
    );
    assert_eq!(
        FluentFaceType::from_node_count(3),
        Some(FluentFaceType::Triangle)
    );
    assert_eq!(
        FluentFaceType::from_node_count(4),
        Some(FluentFaceType::Quad)
    );
    assert_eq!(FluentFaceType::from_node_count(1), None);
    assert_eq!(FluentFaceType::from_node_count(5), None);
}
