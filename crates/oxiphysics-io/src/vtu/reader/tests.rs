// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the VTU reader.

use super::*;
use crate::mesh_io::base64_encode_buffer;
use crate::vtu::{
    VTK_HEX, VTK_LINE, VTK_QUAD, VTK_TET, VTK_TRIANGLE, VtuAnnotation, VtuMultiPiece, VtuWriter,
    VtuWriterLegacy, write_vtu_with_annotations,
};
use oxiphysics_core::math::Vec3;

// ---- fixture helpers ---------------------------------------------------

fn array(name: &str, ty: &str, components: usize, format: &str, content: &str) -> String {
    format!(
        "<DataArray type=\"{ty}\" Name=\"{name}\" NumberOfComponents=\"{components}\" \
         format=\"{format}\">{content}</DataArray>"
    )
}

fn points_array(ty: &str, format: &str, content: &str) -> String {
    format!(
        "<DataArray type=\"{ty}\" NumberOfComponents=\"3\" format=\"{format}\">{content}</DataArray>"
    )
}

/// A one-triangle document; `vtk_attrs` are extra `VTKFile` attributes.
fn triangle_doc(
    points: &str,
    connectivity: &str,
    offsets: &str,
    types: &str,
    vtk_attrs: &str,
) -> String {
    format!(
        "<?xml version=\"1.0\"?>\n\
         <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\"{vtk_attrs}>\n\
         <UnstructuredGrid>\n<Piece NumberOfPoints=\"3\" NumberOfCells=\"1\">\n\
         <Points>{points}</Points>\n\
         <Cells>{connectivity}{offsets}{types}</Cells>\n\
         </Piece>\n</UnstructuredGrid>\n</VTKFile>\n"
    )
}

fn ascii_triangle_doc() -> String {
    triangle_doc(
        &points_array("Float64", "ascii", "0 0 0 1 0 0 0 1 0"),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    )
}

fn f32_bytes(values: &[f32], big: bool) -> Vec<u8> {
    values
        .iter()
        .flat_map(|v| {
            if big {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            }
        })
        .collect()
}
fn f64_bytes(values: &[f64], big: bool) -> Vec<u8> {
    values
        .iter()
        .flat_map(|v| {
            if big {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            }
        })
        .collect()
}
fn i32_bytes(values: &[i32], big: bool) -> Vec<u8> {
    values
        .iter()
        .flat_map(|v| {
            if big {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            }
        })
        .collect()
}
fn i64_bytes(values: &[i64], big: bool) -> Vec<u8> {
    values
        .iter()
        .flat_map(|v| {
            if big {
                v.to_be_bytes()
            } else {
                v.to_le_bytes()
            }
        })
        .collect()
}

/// Length header for `payload`.
fn header_bytes(payload: &[u8], header64: bool, big: bool) -> Vec<u8> {
    match (header64, big) {
        (false, false) => (payload.len() as u32).to_le_bytes().to_vec(),
        (false, true) => (payload.len() as u32).to_be_bytes().to_vec(),
        (true, false) => (payload.len() as u64).to_le_bytes().to_vec(),
        (true, true) => (payload.len() as u64).to_be_bytes().to_vec(),
    }
}

/// VTK inline binary layout: header and payload are base64-encoded
/// separately and concatenated.
fn vtk_binary(payload: &[u8], header64: bool, big: bool) -> String {
    format!(
        "{}{}",
        base64_encode_buffer(&header_bytes(payload, header64, big)),
        base64_encode_buffer(payload)
    )
}

const TET_POINTS: [f32; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// One tetrahedron with Float32 points, a Float64 point scalar and a
/// Float32 cell vector, everything inline-binary.
fn tetra_binary_doc(header64: bool, big: bool) -> String {
    let enc = |payload: Vec<u8>| vtk_binary(&payload, header64, big);
    let (conn, offs) = if header64 {
        (
            enc(i64_bytes(&[0, 1, 2, 3], big)),
            enc(i64_bytes(&[4], big)),
        )
    } else {
        (
            enc(i32_bytes(&[0, 1, 2, 3], big)),
            enc(i32_bytes(&[4], big)),
        )
    };
    let int_type = if header64 { "Int64" } else { "Int32" };
    let header_attr = if header64 { "UInt64" } else { "UInt32" };
    let order = if big { "BigEndian" } else { "LittleEndian" };
    format!(
        "<?xml version=\"1.0\"?>\n\
         <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"{order}\" header_type=\"{header_attr}\">\n\
         <UnstructuredGrid><Piece NumberOfPoints=\"4\" NumberOfCells=\"1\">\n\
         <PointData>{pd}</PointData><CellData>{cd}</CellData>\n\
         <Points>{pts}</Points>\n\
         <Cells>{c}{o}{t}</Cells></Piece></UnstructuredGrid></VTKFile>\n",
        pd = array(
            "temperature",
            "Float64",
            1,
            "binary",
            &enc(f64_bytes(&[100.0, 200.0, 300.0, 400.0], big))
        ),
        cd = array(
            "flux",
            "Float32",
            3,
            "binary",
            &enc(f32_bytes(&[0.5, 1.5, -2.5], big))
        ),
        pts = points_array("Float32", "binary", &enc(f32_bytes(&TET_POINTS, big))),
        c = array("connectivity", int_type, 1, "binary", &conn),
        o = array("offsets", int_type, 1, "binary", &offs),
        t = array("types", "UInt8", 1, "binary", &enc(vec![10u8])),
    )
}

fn assert_tetra(reader: &VtuReader) {
    assert_eq!(reader.num_points(), 4);
    assert_eq!(reader.points()[1], [1.0, 0.0, 0.0]);
    assert_eq!(reader.points()[3], [0.0, 0.0, 1.0]);
    assert_eq!(reader.num_cells(), 1);
    assert_eq!(reader.cell_types(), &[10u8]);
    assert_eq!(reader.connectivity(), &[0usize, 1, 2, 3]);
    assert_eq!(reader.offsets(), &[4usize]);
    assert_eq!(reader.cell_nodes(0), Some(&[0usize, 1, 2, 3][..]));
    let temperature = reader.point_data_array("temperature").unwrap();
    assert_eq!(temperature.num_components, 1);
    assert_eq!(temperature.values, vec![100.0, 200.0, 300.0, 400.0]);
    let flux = reader.cell_data_array("flux").unwrap();
    assert_eq!(flux.num_components, 3);
    assert_eq!(flux.num_tuples(), 1);
    assert_eq!(flux.values, vec![0.5, 1.5, -2.5]);
    assert_eq!(reader.point_data_names(), &["temperature".to_string()]);
    assert_eq!(reader.cell_data_names(), &["flux".to_string()]);
}

fn parse_error(xml: &str) -> String {
    match VtuReader::from_xml(xml) {
        Err(Error::Parse(message)) => message,
        Err(other) => panic!("expected Error::Parse, got {other:?}"),
        Ok(_) => panic!("expected an error"),
    }
}

// ---- writer -> reader round trips ---------------------------------------

fn mixed_writer() -> VtuWriter {
    let mut w = VtuWriter::new();
    w.add_points(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
        [2.5, 0.125, -3.75],
    ]);
    w.add_hex(0, 1, 2, 3, 4, 5, 6, 7);
    w.add_tet(0, 1, 3, 4);
    w.add_quad(0, 1, 2, 3);
    w.add_triangle(5, 6, 8);
    w.add_cell(vec![7, 8], VTK_LINE);
    w.add_point_data_scalar("pressure", &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.5]);
    w.add_point_data_vector(
        "velocity",
        &[
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.5, 0.5, 0.5],
            [-1.0, -2.0, -3.0],
        ],
    );
    w.add_cell_data_scalar("quality", &[0.9, 0.8, 0.7, 0.6, 0.5]);
    w.add_cell_data_vector(
        "normal",
        &[
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.6, 0.8],
        ],
    );
    w
}

#[test]
fn roundtrip_preserves_connectivity_offsets_and_types() {
    let w = mixed_writer();
    let reader = VtuReader::from_xml(&w.to_xml()).unwrap();
    assert_eq!(reader.points(), w.points());
    assert_eq!(reader.num_cells(), 5);
    assert_eq!(reader.cell_types(), w.cell_types());
    assert_eq!(reader.connectivity(), w.cells_flat().as_slice());
    assert_eq!(reader.offsets(), &[8usize, 12, 16, 19, 21]);
    for (cell, nodes) in w.cells_per_cell().iter().enumerate() {
        assert_eq!(
            reader.cell_nodes(cell),
            Some(nodes.as_slice()),
            "cell {cell}"
        );
    }
    assert_eq!(reader.cell_nodes(5), None);
    assert_eq!(reader.cell_types()[0], VTK_HEX);
    assert_eq!(reader.cell_types()[1], VTK_TET);
    assert_eq!(reader.cell_types()[2], VTK_QUAD);
    assert_eq!(reader.cell_types()[3], VTK_TRIANGLE);
    assert_eq!(reader.cell_types()[4], VTK_LINE);
}

#[test]
fn roundtrip_preserves_point_and_cell_data_values() {
    let reader = VtuReader::from_xml(&mixed_writer().to_xml()).unwrap();
    let pressure = reader.point_data_array("pressure").unwrap();
    assert_eq!(pressure.num_components, 1);
    assert_eq!(pressure.values.len(), 9);
    assert_eq!(pressure.values[8], 9.5);
    let velocity = reader.point_data_array("velocity").unwrap();
    assert_eq!(velocity.num_components, 3);
    assert_eq!(velocity.num_tuples(), 9);
    assert_eq!(&velocity.values[24..], &[-1.0, -2.0, -3.0]);
    let quality = reader.cell_data_array("quality").unwrap();
    assert_eq!(quality.values, vec![0.9, 0.8, 0.7, 0.6, 0.5]);
    let normal = reader.cell_data_array("normal").unwrap();
    assert_eq!(normal.num_tuples(), 5);
    assert_eq!(&normal.values[12..], &[0.0, 0.6, 0.8]);
    assert_eq!(
        reader.point_data_names(),
        &["pressure".to_string(), "velocity".to_string()]
    );
    assert_eq!(
        reader.cell_data_names(),
        &["quality".to_string(), "normal".to_string()]
    );
    assert!(reader.point_data_array("missing").is_none());
    assert!(reader.cell_data_array("pressure").is_none());
}

#[test]
fn roundtrip_empty_grid() {
    let reader = VtuReader::from_xml(&VtuWriter::new().to_xml()).unwrap();
    assert_eq!(reader.num_points(), 0);
    assert_eq!(reader.num_cells(), 0);
    assert!(reader.connectivity().is_empty());
    assert!(reader.offsets().is_empty());
    assert!(reader.cell_types().is_empty());
    assert!(reader.bounding_box().is_none());
}

#[test]
fn roundtrip_points_only_grid_keeps_bounding_box() {
    let mut w = VtuWriter::new();
    w.add_points(&[[-1.0, -2.0, -3.0], [4.0, 5.0, 6.0]]);
    let reader = VtuReader::from_xml(&w.to_xml()).unwrap();
    let (min, max) = reader.bounding_box().unwrap();
    assert_eq!(min, [-1.0, -2.0, -3.0]);
    assert_eq!(max, [4.0, 5.0, 6.0]);
}

#[test]
fn legacy_writer_file_is_read_back_from_disk() {
    let path = std::env::temp_dir().join("oxiphysics_vtu_reader_legacy.vtu");
    let pts = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 1.0, 1.0),
    ];
    let cells = vec![[0, 1, 2, 3], [1, 2, 3, 4]];
    let temperature = [10.0, 20.0, 30.0, 40.0, 50.0];
    let volume = [0.25, 0.75];
    VtuWriterLegacy::write(
        path.to_str().unwrap_or(""),
        &pts,
        &cells,
        &[("temperature", temperature.as_slice())],
        &[("volume", volume.as_slice())],
    )
    .unwrap();
    let reader = VtuReader::from_file(&path).unwrap();
    std::fs::remove_file(&path).ok();
    assert_eq!(reader.num_points(), 5);
    assert_eq!(reader.points()[4], [1.0, 1.0, 1.0]);
    assert_eq!(reader.cell_types(), &[10u8, 10]);
    assert_eq!(reader.connectivity(), &[0usize, 1, 2, 3, 1, 2, 3, 4]);
    assert_eq!(reader.offsets(), &[4usize, 8]);
    assert_eq!(
        reader.point_data_array("temperature").unwrap().values,
        temperature.to_vec()
    );
    assert_eq!(
        reader.cell_data_array("volume").unwrap().values,
        volume.to_vec()
    );
}

#[test]
fn xml_vtu_writer_output_with_int64_arrays_is_read_back() {
    use crate::vtk_writer::{VtkCellTypeW, VtkUnstructuredGrid, XmlVtuWriter};
    let mut grid = VtkUnstructuredGrid::new();
    for p in [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
    ] {
        grid.add_point(p);
    }
    grid.add_cell(vec![0, 1, 2, 3], VtkCellTypeW::Tetrahedron);
    grid.add_cell(vec![1, 2, 3, 4], VtkCellTypeW::Tetrahedron);
    grid.add_cell(vec![0, 4], VtkCellTypeW::Line);
    grid.add_point_scalars("temperature", vec![1.0, 2.0, 3.0, 4.0, 5.0]);
    grid.add_point_vectors("velocity", vec![[1.0, 0.0, 0.0]; 5]);
    grid.add_cell_scalars("volume", vec![0.25, 0.5, 0.0]);
    let reader = VtuReader::from_xml(&XmlVtuWriter::write(&grid)).unwrap();
    assert_eq!(reader.num_points(), 5);
    assert_eq!(reader.cell_types(), &[10u8, 10, 3]);
    assert_eq!(reader.connectivity(), &[0usize, 1, 2, 3, 1, 2, 3, 4, 0, 4]);
    assert_eq!(reader.offsets(), &[4usize, 8, 10]);
    assert_eq!(
        reader.point_data_array("temperature").unwrap().values,
        vec![1.0, 2.0, 3.0, 4.0, 5.0]
    );
    assert_eq!(
        reader.point_data_array("velocity").unwrap().num_components,
        3
    );
    assert_eq!(
        reader.cell_data_array("volume").unwrap().values,
        vec![0.25, 0.5, 0.0]
    );
}

#[test]
fn annotated_output_still_parses() {
    let mut w = VtuWriter::new();
    w.add_points(&[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    w.add_triangle(0, 1, 2);
    let annotations = vec![VtuAnnotation {
        key: "author".into(),
        value: "OxiPhysics".into(),
    }];
    let xml = write_vtu_with_annotations(&w, &annotations);
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert_eq!(reader.num_points(), 3);
    assert_eq!(reader.connectivity(), &[0usize, 1, 2]);
}

// ---- hand written fixtures ---------------------------------------------

#[test]
fn float32_points_fixture_is_read() {
    let xml = triangle_doc(
        &points_array("Float32", "ascii", "0 0 0  1.5 0 0  0 -2.25 0"),
        &array("connectivity", "Int64", 1, "ascii", "0 1 2"),
        &array("offsets", "Int64", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert_eq!(
        reader.points(),
        &[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, -2.25, 0.0]]
    );
    assert_eq!(reader.connectivity(), &[0usize, 1, 2]);
    assert_eq!(reader.cell_types(), &[VTK_TRIANGLE]);
}

#[test]
fn single_line_document_is_parsed() {
    let one_line = ascii_triangle_doc().replace('\n', "");
    assert!(!one_line.contains('\n'));
    let reader = VtuReader::from_xml(&one_line).unwrap();
    assert_eq!(reader.num_points(), 3);
    assert_eq!(reader.points()[1], [1.0, 0.0, 0.0]);
    assert_eq!(reader.connectivity(), &[0usize, 1, 2]);
}

#[test]
fn attribute_order_and_extra_attributes_do_not_matter() {
    let xml = "<VTKFile byte_order=\"LittleEndian\" version=\"1.0\" type=\"UnstructuredGrid\">\
        <UnstructuredGrid><Piece NumberOfCells=\"1\" NumberOfPoints=\"3\" Extra=\"x\">\
        <Points><DataArray format=\"ascii\" NumberOfComponents=\"3\" RangeMin=\"0\" type=\"Float64\">\
        0 0 0\n1 0 0\n0\n1\n0</DataArray></Points>\
        <Cells>\
        <DataArray format=\"ascii\" type=\"Int32\" Name=\"offsets\">3</DataArray>\
        <DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\">5</DataArray>\
        <DataArray Name=\"connectivity\" format=\"ascii\" type=\"Int32\">0 1 2</DataArray>\
        </Cells></Piece></UnstructuredGrid></VTKFile>";
    let reader = VtuReader::from_xml(xml).unwrap();
    // The coordinate list may wrap at arbitrary places.
    assert_eq!(
        reader.points(),
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    );
    assert_eq!(reader.connectivity(), &[0usize, 1, 2]);
    assert_eq!(reader.offsets(), &[3usize]);
}

#[test]
fn cdata_and_character_references_in_arrays_are_decoded() {
    let xml = triangle_doc(
        &points_array("Float64", "ascii", "<![CDATA[0 0 0 1 0 0 0 1 0]]>"),
        &array("connectivity", "Int32", 1, "ascii", "0&#32;1&#32;2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert_eq!(reader.points()[2], [0.0, 1.0, 0.0]);
    assert_eq!(reader.connectivity(), &[0usize, 1, 2]);
}

#[test]
fn named_entity_inside_numeric_data_is_an_error() {
    let xml = triangle_doc(
        &points_array("Float64", "ascii", "0 0 0 1 0 0 0 1 &amp; 0"),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    assert!(parse_error(&xml).contains("entity"));
}

#[test]
fn data_array_names_are_unescaped_and_unnamed_arrays_skipped() {
    let xml = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\
        <UnstructuredGrid><Piece NumberOfPoints=\"2\" NumberOfCells=\"0\">\
        <PointData>\
        <DataArray type=\"Float64\" Name=\"a&amp;b &lt;1&gt;\" format=\"ascii\">1 2</DataArray>\
        <DataArray type=\"Float64\" format=\"ascii\">8 9</DataArray>\
        </PointData>\
        <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">0 0 0 1 1 1</DataArray></Points>\
        <Cells></Cells></Piece></UnstructuredGrid></VTKFile>";
    let reader = VtuReader::from_xml(xml).unwrap();
    assert_eq!(reader.point_data_names(), &["a&b <1>".to_string()]);
    assert_eq!(
        reader.point_data_array("a&b <1>").unwrap().values,
        vec![1.0, 2.0]
    );
}

#[test]
fn field_data_and_unknown_cell_arrays_are_ignored() {
    let xml = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\
        <UnstructuredGrid><FieldData>\
        <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" offset=\"0\"/>\
        </FieldData>\
        <Piece NumberOfPoints=\"1\" NumberOfCells=\"1\">\
        <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">1 2 3</DataArray></Points>\
        <Cells>\
        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">0</DataArray>\
        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">1</DataArray>\
        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">1</DataArray>\
        <DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">not numbers here</DataArray>\
        </Cells></Piece></UnstructuredGrid></VTKFile>";
    let reader = VtuReader::from_xml(xml).unwrap();
    assert_eq!(reader.points(), &[[1.0, 2.0, 3.0]]);
    assert_eq!(reader.cell_types(), &[1u8]);
}

#[test]
fn offsets_with_leading_zero_are_accepted() {
    let xml = triangle_doc(
        &points_array("Float64", "ascii", "0 0 0 1 0 0 0 1 0"),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "0 3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert_eq!(reader.offsets(), &[3usize]);
    assert_eq!(reader.cell_nodes(0), Some(&[0usize, 1, 2][..]));
}

#[test]
fn nan_and_infinity_tokens_are_decoded() {
    let xml = triangle_doc(
        &points_array("Float64", "ascii", "nan inf -inf 0 0 0 1 1 1"),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert!(reader.points()[0][0].is_nan());
    assert_eq!(reader.points()[0][1], f64::INFINITY);
    assert_eq!(reader.points()[0][2], f64::NEG_INFINITY);
}

// ---- inline binary -----------------------------------------------------

#[test]
fn binary_arrays_with_uint32_header_little_endian() {
    let reader = VtuReader::from_xml(&tetra_binary_doc(false, false)).unwrap();
    assert_tetra(&reader);
}

#[test]
fn binary_arrays_with_uint64_header_and_int64_connectivity() {
    let reader = VtuReader::from_xml(&tetra_binary_doc(true, false)).unwrap();
    assert_tetra(&reader);
}

#[test]
fn binary_arrays_big_endian() {
    assert_tetra(&VtuReader::from_xml(&tetra_binary_doc(false, true)).unwrap());
    assert_tetra(&VtuReader::from_xml(&tetra_binary_doc(true, true)).unwrap());
}

#[test]
fn binary_fixture_from_an_independent_encoder_is_read() {
    // The base64 literals were produced by Python's `base64` and `struct`
    // modules (header and payload encoded separately, as VTK does), so the
    // layout is checked against something other than this crate's encoder.
    let xml = "<?xml version=\"1.0\"?>\n\
        <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n\
        <UnstructuredGrid><Piece NumberOfPoints=\"4\" NumberOfCells=\"1\">\n\
        <PointData><DataArray type=\"Float64\" Name=\"temperature\" format=\"binary\">\
        IAAAAA==AAAAAAAAWUAAAAAAAABpQAAAAAAAwHJAAAAAAAAAeUA=</DataArray></PointData>\n\
        <Points><DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\
        MAAAAA==AAAAAAAAAAAAAAAAAACAPwAAAAAAAAAAAAAAAAAAgD8AAAAAAAAAAAAAAAAAAIA/\
        </DataArray></Points>\n\
        <Cells>\
        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">EAAAAA==AAAAAAEAAAACAAAAAwAAAA==</DataArray>\
        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">BAAAAA==BAAAAA==</DataArray>\
        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">AQAAAA==Cg==</DataArray>\
        </Cells></Piece></UnstructuredGrid></VTKFile>\n";
    let reader = VtuReader::from_xml(xml).unwrap();
    assert_eq!(reader.num_points(), 4);
    assert_eq!(
        reader.points(),
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]
        ]
    );
    assert_eq!(reader.connectivity(), &[0usize, 1, 2, 3]);
    assert_eq!(reader.offsets(), &[4usize]);
    assert_eq!(reader.cell_types(), &[10u8]);
    let temperature = reader.point_data_array("temperature").unwrap();
    assert_eq!(temperature.values, vec![100.0, 200.0, 300.0, 400.0]);
}

#[test]
fn binary_header_and_payload_in_one_base64_stream_is_accepted() {
    let payload = f64_bytes(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0], false);
    let mut all = header_bytes(&payload, false, false);
    all.extend_from_slice(&payload);
    let xml = triangle_doc(
        &points_array("Float64", "binary", &base64_encode_buffer(&all)),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert_eq!(reader.points()[1], [1.0, 0.0, 0.0]);
    assert_eq!(reader.points()[2], [0.0, 1.0, 0.0]);
}

#[test]
fn binary_base64_text_may_be_wrapped_over_lines() {
    let payload = f64_bytes(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0], false);
    let text = vtk_binary(&payload, false, false);
    let wrapped: String = text
        .as_bytes()
        .chunks(16)
        .map(|c| format!("   {}\n", String::from_utf8_lossy(c)))
        .collect();
    let xml = triangle_doc(
        &points_array("Float64", "binary", &wrapped),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert_eq!(reader.points()[1], [1.0, 0.0, 0.0]);
}

#[test]
fn binary_with_wrong_length_header_is_an_error() {
    let payload = f64_bytes(&[0.0; 9], false);
    let mut header = header_bytes(&payload, false, false);
    header[0] = header[0].wrapping_add(8);
    let text = format!(
        "{}{}",
        base64_encode_buffer(&header),
        base64_encode_buffer(&payload)
    );
    let xml = triangle_doc(
        &points_array("Float64", "binary", &text),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    assert!(parse_error(&xml).contains("length header"));
}

#[test]
fn truncated_binary_payload_is_an_error() {
    let payload = f64_bytes(&[0.0; 9], false);
    let text = vtk_binary(&payload, false, false);
    let truncated = &text[..text.len() - 8];
    let xml = triangle_doc(
        &points_array("Float64", "binary", truncated),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    assert!(!parse_error(&xml).is_empty());
}

#[test]
fn invalid_base64_is_an_error() {
    let xml = triangle_doc(
        &points_array("Float64", "binary", "!!!!!!!!@@@@####"),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    assert!(!parse_error(&xml).is_empty());
}

#[test]
fn binary_payload_that_is_not_a_whole_number_of_elements_is_an_error() {
    // 5 payload bytes cannot hold Float32 values.
    let payload = [1u8, 2, 3, 4, 5];
    let xml = triangle_doc(
        &points_array("Float32", "binary", &vtk_binary(&payload, false, false)),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    assert!(parse_error(&xml).contains("multiple"));
}

#[test]
fn empty_binary_array_is_an_empty_array() {
    let doc = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\
        <UnstructuredGrid><Piece NumberOfPoints=\"0\" NumberOfCells=\"0\">\
        <Points><DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\"></DataArray></Points>\
        <Cells/></Piece></UnstructuredGrid></VTKFile>";
    let reader = VtuReader::from_xml(doc).unwrap();
    assert_eq!(reader.num_points(), 0);
}

// ---- unsupported encodings and malformed input --------------------------

#[test]
fn appended_data_arrays_are_rejected_not_silently_dropped() {
    let xml = triangle_doc(
        "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>",
        "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"36\"/>",
        "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"52\"/>",
        "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"60\"/>",
        " header_type=\"UInt64\"",
    );
    let xml = xml.replace(
        "</VTKFile>",
        "<AppendedData encoding=\"base64\">_AAAA</AppendedData></VTKFile>",
    );
    let message = parse_error(&xml);
    assert!(message.contains("appended"), "{message}");
}

#[test]
fn compressed_binary_data_is_rejected() {
    let payload = f64_bytes(&[0.0; 9], false);
    let xml = triangle_doc(
        &points_array("Float64", "binary", &vtk_binary(&payload, false, false)),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        " compressor=\"vtkZLibDataCompressor\"",
    );
    let message = parse_error(&xml);
    assert!(message.contains("compress"), "{message}");
}

#[test]
fn compressor_attribute_does_not_affect_ascii_arrays() {
    let doc = ascii_triangle_doc().replace(
        "byte_order=\"LittleEndian\"",
        "byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\"",
    );
    let reader = VtuReader::from_xml(&doc).unwrap();
    assert_eq!(reader.num_points(), 3);
}

#[test]
fn other_file_types_and_non_vtk_documents_are_rejected() {
    let polydata = ascii_triangle_doc().replace("UnstructuredGrid\"", "PolyData\"");
    assert!(parse_error(&polydata).contains("PolyData"));
    assert!(parse_error("<root><child/></root>").contains("VTKFile"));
    assert!(parse_error("").contains("VTKFile"));
    assert!(parse_error("plain text, not xml").contains("VTKFile"));
    assert!(parse_error("<VTKFile version=\"1.0\"></VTKFile>").contains("type"));
}

#[test]
fn malformed_xml_and_truncated_documents_are_errors() {
    let doc = ascii_triangle_doc();
    // Mismatched closing tag.
    assert!(!parse_error(&doc.replace("</Points>", "</Pointz>")).is_empty());
    // Cut off in the middle of the document.
    let cut = &doc[..doc.find("<Cells>").unwrap()];
    assert!(parse_error(cut).contains("unclosed") || parse_error(cut).contains("end"));
    // Garbage inside an attribute.
    assert!(!parse_error(&doc.replace("Name=\"offsets\"", "Name=offsets")).is_empty());
}

#[test]
fn invalid_vtkfile_attributes_are_errors() {
    let doc = ascii_triangle_doc();
    let big = doc.replace("LittleEndian", "MiddleEndian");
    assert!(parse_error(&big).contains("byte_order"));
    let header = doc.replace("<VTKFile ", "<VTKFile header_type=\"UInt16\" ");
    assert!(parse_error(&header).contains("header_type"));
}

#[test]
fn non_numeric_ascii_token_is_an_error() {
    let xml = triangle_doc(
        &points_array("Float64", "ascii", "0 0 0 1 0 0 0 one 0"),
        &array("connectivity", "Int32", 1, "ascii", "0 1 2"),
        &array("offsets", "Int32", 1, "ascii", "3"),
        &array("types", "UInt8", 1, "ascii", "5"),
        "",
    );
    let message = parse_error(&xml);
    assert!(message.contains("one"), "{message}");
}

#[test]
fn data_array_attribute_problems_are_errors() {
    let base = ascii_triangle_doc();
    let no_type = base.replacen(" type=\"Float64\"", "", 1);
    assert!(parse_error(&no_type).contains("type"));
    let bad_type = base.replacen("Float64", "Float128", 1);
    assert!(parse_error(&bad_type).contains("Float128"));
    let two_components = base.replacen("NumberOfComponents=\"3\"", "NumberOfComponents=\"2\"", 1);
    assert!(parse_error(&two_components).contains("3 components"));
    let zero_components = base.replacen("NumberOfComponents=\"3\"", "NumberOfComponents=\"0\"", 1);
    assert!(parse_error(&zero_components).contains("zero"));
    let bad_format = base.replacen("format=\"ascii\"", "format=\"hex\"", 1);
    assert!(parse_error(&bad_format).contains("hex"));
    let tuples = base.replacen(
        "NumberOfComponents=\"3\"",
        "NumberOfComponents=\"3\" NumberOfTuples=\"4\"",
        1,
    );
    assert!(parse_error(&tuples).contains("tuples"));
}

#[test]
fn piece_counts_are_validated_against_the_arrays() {
    let doc = ascii_triangle_doc();
    // More points declared than present.
    let points = doc.replace("NumberOfPoints=\"3\"", "NumberOfPoints=\"4\"");
    assert!(parse_error(&points).contains("points"));
    // More cells declared than present.
    let cells = doc.replace("NumberOfCells=\"1\"", "NumberOfCells=\"2\"");
    assert!(parse_error(&cells).contains("cells"));
    // Missing counts.
    let no_count = doc.replace(" NumberOfCells=\"1\"", "");
    assert!(parse_error(&no_count).contains("NumberOfCells"));
    let bad_count = doc.replace("NumberOfPoints=\"3\"", "NumberOfPoints=\"three\"");
    assert!(parse_error(&bad_count).contains("NumberOfPoints"));
}

#[test]
fn connectivity_problems_are_errors() {
    let doc = ascii_triangle_doc();
    // Point index outside the piece.
    let outside = doc.replace(">0 1 2<", ">0 1 3<");
    assert!(parse_error(&outside).contains("point 3"));
    // Last offset disagrees with the connectivity length.
    let offset = doc.replace(">3<", ">2<");
    assert!(parse_error(&offset).contains("offset"));
    // Negative and fractional indices.
    assert!(parse_error(&doc.replace(">0 1 2<", ">0 1 -2<")).contains("not a valid index"));
    assert!(parse_error(&doc.replace(">0 1 2<", ">0 1 1.5<")).contains("not a valid index"));
    // A cell type outside 0..=255.
    assert!(parse_error(&doc.replace(">5<", ">300<")).contains("cell type"));
}

#[test]
fn decreasing_offsets_are_an_error() {
    let xml = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\
        <UnstructuredGrid><Piece NumberOfPoints=\"3\" NumberOfCells=\"2\">\
        <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">0 0 0 1 0 0 0 1 0</DataArray></Points>\
        <Cells>\
        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">0 1 2 0 1</DataArray>\
        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">5 3</DataArray>\
        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">5 3</DataArray>\
        </Cells></Piece></UnstructuredGrid></VTKFile>";
    assert!(parse_error(xml).contains("decreases"));
}

#[test]
fn data_array_lengths_must_match_the_piece() {
    let mut w = VtuWriter::new();
    w.add_points(&[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    // Writer does not validate: two values for three points.
    w.add_point_data_scalar("short", &[1.0, 2.0]);
    let message = parse_error(&w.to_xml());
    assert!(message.contains("short"), "{message}");
}

// ---- several pieces ----------------------------------------------------

fn piece_writer(offset: f64, with_data: bool) -> VtuWriter {
    let mut w = VtuWriter::new();
    w.add_points(&[
        [offset, 0.0, 0.0],
        [offset + 1.0, 0.0, 0.0],
        [offset, 1.0, 0.0],
        [offset, 0.0, 1.0],
    ]);
    w.add_tet(0, 1, 2, 3);
    w.add_triangle(0, 1, 2);
    if with_data {
        w.add_point_data_scalar("id", &[offset, offset + 1.0, offset + 2.0, offset + 3.0]);
        w.add_cell_data_scalar("tag", &[offset * 10.0, offset * 10.0 + 1.0]);
    }
    w
}

#[test]
fn multiple_pieces_are_concatenated_with_shifted_indices() {
    let a = piece_writer(0.0, true);
    let b = piece_writer(10.0, true);
    let xml = VtuMultiPiece::write(&[&a, &b]);
    let reader = VtuReader::from_xml(&xml).unwrap();
    assert_eq!(reader.num_points(), 8);
    assert_eq!(reader.num_cells(), 4);
    assert_eq!(reader.points()[4], [10.0, 0.0, 0.0]);
    assert_eq!(
        reader.connectivity(),
        &[0usize, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7, 4, 5, 6]
    );
    assert_eq!(reader.offsets(), &[4usize, 7, 11, 14]);
    assert_eq!(reader.cell_nodes(2), Some(&[4usize, 5, 6, 7][..]));
    assert_eq!(
        reader.cell_types(),
        &[VTK_TET, VTK_TRIANGLE, VTK_TET, VTK_TRIANGLE]
    );
    assert_eq!(
        reader.point_data_array("id").unwrap().values,
        vec![0.0, 1.0, 2.0, 3.0, 10.0, 11.0, 12.0, 13.0]
    );
    assert_eq!(
        reader.cell_data_array("tag").unwrap().values,
        vec![0.0, 1.0, 100.0, 101.0]
    );
    assert_eq!(reader.point_data_names(), &["id".to_string()]);
}

#[test]
fn data_array_missing_from_a_later_piece_is_an_error() {
    let a = piece_writer(0.0, true);
    let b = piece_writer(10.0, false);
    let message = parse_error(&VtuMultiPiece::write(&[&a, &b]));
    assert!(
        message.contains("id") || message.contains("tag"),
        "{message}"
    );
}

#[test]
fn data_array_missing_from_the_first_piece_is_an_error() {
    let a = piece_writer(0.0, false);
    let b = piece_writer(10.0, true);
    let message = parse_error(&VtuMultiPiece::write(&[&a, &b]));
    assert!(
        message.contains("id") || message.contains("tag"),
        "{message}"
    );
}

#[test]
fn empty_pieces_do_not_disturb_the_merge() {
    let empty = VtuWriter::new();
    let b = piece_writer(10.0, true);
    let reader = VtuReader::from_xml(&VtuMultiPiece::write(&[&empty, &b, &empty])).unwrap();
    assert_eq!(reader.num_points(), 4);
    assert_eq!(reader.connectivity(), &[0usize, 1, 2, 3, 0, 1, 2]);
    assert_eq!(reader.point_data_array("id").unwrap().values.len(), 4);
}

// ---- file access --------------------------------------------------------

#[test]
fn from_file_reports_missing_files_as_io_errors() {
    let path = std::env::temp_dir().join("oxiphysics_vtu_reader_does_not_exist.vtu");
    std::fs::remove_file(&path).ok();
    assert!(matches!(VtuReader::from_file(&path), Err(Error::Io(_))));
}

#[test]
fn from_file_rejects_raw_appended_binary_files_with_a_clear_message() {
    let path = std::env::temp_dir().join("oxiphysics_vtu_reader_raw_appended.vtu");
    let mut bytes = b"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" \
        version=\"1.0\"><UnstructuredGrid></UnstructuredGrid>\n\
        <AppendedData encoding=\"raw\">\n   _"
        .to_vec();
    bytes.extend_from_slice(&[0xff, 0xfe, 0x00, 0x80, 0xc3, 0x28]);
    bytes.extend_from_slice(b"\n</AppendedData>\n</VTKFile>\n");
    std::fs::write(&path, &bytes).unwrap();
    let result = VtuReader::from_file(&path);
    std::fs::remove_file(&path).ok();
    match result {
        Err(Error::Parse(message)) => assert!(message.contains("appended"), "{message}"),
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn from_file_rejects_other_non_utf8_content() {
    let path = std::env::temp_dir().join("oxiphysics_vtu_reader_not_utf8.vtu");
    std::fs::write(&path, [0xffu8, 0xfe, 0xfd]).unwrap();
    let result = VtuReader::from_file(&path);
    std::fs::remove_file(&path).ok();
    assert!(matches!(result, Err(Error::Parse(_))));
}

// ---- read_cell_data -----------------------------------------------------

#[test]
fn read_cell_data_works_on_single_line_and_vector_arrays() {
    let xml = mixed_writer().to_xml();
    let flat = xml.replace('\n', " ");
    let normal = VtuReader::read_cell_data(&flat, "normal").unwrap();
    assert_eq!(normal.len(), 15);
    assert_eq!(&normal[12..], &[0.0, 0.6, 0.8]);
    let quality = VtuReader::read_cell_data(&xml, "quality").unwrap();
    assert_eq!(quality, vec![0.9, 0.8, 0.7, 0.6, 0.5]);
}

#[test]
fn read_cell_data_reports_missing_arrays_and_bad_documents() {
    let xml = mixed_writer().to_xml();
    let err = VtuReader::read_cell_data(&xml, "nope").unwrap_err();
    assert!(err.contains("nope"));
    // Point data is not cell data.
    assert!(VtuReader::read_cell_data(&xml, "pressure").is_err());
    assert!(VtuReader::read_cell_data("not xml", "x").is_err());
}

#[test]
fn read_cell_data_decodes_binary_arrays() {
    let reader_doc = tetra_binary_doc(false, false);
    let values = VtuReader::read_cell_data(&reader_doc, "flux").unwrap();
    assert_eq!(values, vec![0.5, 1.5, -2.5]);
}

// ---- low-level helpers --------------------------------------------------

#[test]
fn index_and_cell_type_conversion_reject_invalid_values() {
    assert_eq!(to_index(0.0, "x").unwrap(), 0);
    assert_eq!(to_index(41.0, "x").unwrap(), 41);
    assert!(to_index(-1.0, "x").is_err());
    assert!(to_index(0.5, "x").is_err());
    assert!(to_index(f64::NAN, "x").is_err());
    assert!(to_index(f64::INFINITY, "x").is_err());
    assert!(to_index(1.0e300, "x").is_err());
    assert_eq!(to_cell_type(12.0).unwrap(), 12);
    assert_eq!(to_cell_type(255.0).unwrap(), 255);
    assert!(to_cell_type(256.0).is_err());
    assert!(to_cell_type(-1.0).is_err());
    assert!(to_cell_type(5.5).is_err());
}

#[test]
fn scalar_types_decode_in_both_byte_orders() {
    let cases: Vec<(ScalarType, Vec<u8>, Vec<u8>, f64)> = vec![
        (ScalarType::Int8, vec![0xfe], vec![0xfe], -2.0),
        (ScalarType::UInt8, vec![0xfe], vec![0xfe], 254.0),
        (
            ScalarType::Int16,
            (-300i16).to_le_bytes().to_vec(),
            (-300i16).to_be_bytes().to_vec(),
            -300.0,
        ),
        (
            ScalarType::UInt16,
            60000u16.to_le_bytes().to_vec(),
            60000u16.to_be_bytes().to_vec(),
            60000.0,
        ),
        (
            ScalarType::Int32,
            (-70000i32).to_le_bytes().to_vec(),
            (-70000i32).to_be_bytes().to_vec(),
            -70000.0,
        ),
        (
            ScalarType::UInt32,
            4_000_000_000u32.to_le_bytes().to_vec(),
            4_000_000_000u32.to_be_bytes().to_vec(),
            4_000_000_000.0,
        ),
        (
            ScalarType::Int64,
            (-5_000_000_000i64).to_le_bytes().to_vec(),
            (-5_000_000_000i64).to_be_bytes().to_vec(),
            -5.0e9,
        ),
        (
            ScalarType::UInt64,
            7_000_000_000u64.to_le_bytes().to_vec(),
            7_000_000_000u64.to_be_bytes().to_vec(),
            7.0e9,
        ),
        (
            ScalarType::Float32,
            2.5f32.to_le_bytes().to_vec(),
            2.5f32.to_be_bytes().to_vec(),
            2.5,
        ),
        (
            ScalarType::Float64,
            (-0.125f64).to_le_bytes().to_vec(),
            (-0.125f64).to_be_bytes().to_vec(),
            -0.125,
        ),
    ];
    for (ty, le, be, expected) in cases {
        assert_eq!(le.len(), ty.byte_size(), "{ty:?}");
        assert_eq!(ty.decode(&le, false), Some(expected), "{ty:?} LE");
        assert_eq!(ty.decode(&be, true), Some(expected), "{ty:?} BE");
    }
    // A chunk that is too short is rejected instead of panicking.
    assert_eq!(ScalarType::Float64.decode(&[1, 2, 3], false), None);
}
