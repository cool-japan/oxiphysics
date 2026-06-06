//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{PvtuPiece, VtkCellTypeW, VtkUnstructuredGrid, VtuValidationResult};

/// Build a [`VtkUnstructuredGrid`] from a tetrahedral mesh.
pub fn tet_mesh_to_vtk(points: &[[f64; 3]], tets: &[[usize; 4]]) -> VtkUnstructuredGrid {
    let mut grid = VtkUnstructuredGrid::new();
    for &p in points {
        grid.add_point(p);
    }
    for &[a, b, c, d] in tets {
        grid.add_cell(vec![a, b, c, d], VtkCellTypeW::Tetrahedron);
    }
    grid
}
/// Build a [`VtkUnstructuredGrid`] from a triangle mesh.
pub fn tri_mesh_to_vtk(points: &[[f64; 3]], tris: &[[usize; 3]]) -> VtkUnstructuredGrid {
    let mut grid = VtkUnstructuredGrid::new();
    for &p in points {
        grid.add_point(p);
    }
    for &[a, b, c] in tris {
        grid.add_cell(vec![a, b, c], VtkCellTypeW::Triangle);
    }
    grid
}
/// Build a [`VtkUnstructuredGrid`] from a hexahedral mesh.
pub fn hex_mesh_to_vtk(points: &[[f64; 3]], hexes: &[[usize; 8]]) -> VtkUnstructuredGrid {
    let mut grid = VtkUnstructuredGrid::new();
    for &p in points {
        grid.add_point(p);
    }
    for &[a, b, c, d, e, f, g, h] in hexes {
        grid.add_cell(vec![a, b, c, d, e, f, g, h], VtkCellTypeW::Hexahedron);
    }
    grid
}
/// Validate data arrays attached to a `VtkUnstructuredGrid`.
///
/// Returns a list of error strings (empty if valid).
pub fn validate_vtk_grid(grid: &VtkUnstructuredGrid) -> Vec<String> {
    let mut errors = Vec::new();
    let n_pts = grid.n_points();
    let n_cells = grid.n_cells();
    for arr in &grid.point_data {
        if arr.len() != n_pts {
            errors.push(format!(
                "PointData '{}': len {} != n_points {}",
                arr.name(),
                arr.len(),
                n_pts
            ));
        }
    }
    for arr in &grid.cell_data {
        if arr.len() != n_cells {
            errors.push(format!(
                "CellData '{}': len {} != n_cells {}",
                arr.name(),
                arr.len(),
                n_cells
            ));
        }
    }
    for (ci, conn) in grid.cells.iter().enumerate() {
        if conn.is_empty() {
            errors.push(format!("Cell {ci} has empty connectivity"));
        }
        for &idx in conn {
            if idx >= n_pts {
                errors.push(format!("Cell {ci}: index {idx} >= n_points {n_pts}"));
            }
        }
    }
    errors
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::FormattedAsciiWriter;
    use crate::LegacyVtkWriter;
    use crate::ParticleVtkExporter;

    use crate::VtkBinaryWriter;
    use crate::VtkCompression;
    use crate::VtkParallelWriter;
    use crate::VtkPartition;

    use crate::XmlVtuWriter;
    use crate::vtk_writer::types::*;
    #[test]
    fn test_add_point_and_n_points() {
        let mut grid = VtkUnstructuredGrid::new();
        assert_eq!(grid.n_points(), 0);
        let idx0 = grid.add_point([1.0, 2.0, 3.0]);
        let idx1 = grid.add_point([4.0, 5.0, 6.0]);
        assert_eq!(idx0, 0);
        assert_eq!(idx1, 1);
        assert_eq!(grid.n_points(), 2);
    }
    #[test]
    fn test_legacy_writer_contains_dataset_unstructured_grid() {
        let mut grid = VtkUnstructuredGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_cell(vec![0], VtkCellTypeW::Vertex);
        let output = LegacyVtkWriter::write(&grid);
        assert!(
            output.contains("DATASET UNSTRUCTURED_GRID"),
            "expected 'DATASET UNSTRUCTURED_GRID' in output"
        );
    }
    #[test]
    fn test_xml_vtu_writer_contains_vtkfile() {
        let mut grid = VtkUnstructuredGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_cell(vec![0], VtkCellTypeW::Vertex);
        let output = XmlVtuWriter::write(&grid);
        assert!(output.contains("<VTKFile"), "expected '<VTKFile' in output");
    }
    #[test]
    fn test_encode_base64_known_string() {
        let encoded = XmlVtuWriter::encode_base64(b"Man");
        assert_eq!(encoded, "TWFu");
        let encoded2 = XmlVtuWriter::encode_base64(b"Ma");
        assert_eq!(encoded2, "TWE=");
        let encoded3 = XmlVtuWriter::encode_base64(b"M");
        assert_eq!(encoded3, "TQ==");
    }
    #[test]
    fn test_export_particles_contains_points() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let output = ParticleVtkExporter::export_particles(&positions, None, None);
        assert!(
            output.contains("POINTS"),
            "expected 'POINTS' in particle output"
        );
    }
    #[test]
    fn test_tet_mesh_to_vtk_cell_type_is_tetrahedron() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let tets = [[0, 1, 2, 3]];
        let grid = tet_mesh_to_vtk(&pts, &tets);
        assert_eq!(grid.n_cells(), 1);
        assert_eq!(grid.cell_types[0], VtkCellTypeW::Tetrahedron);
    }
    #[test]
    fn test_write_point_data_contains_field_name() {
        let mut grid = VtkUnstructuredGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_cell(vec![0, 1], VtkCellTypeW::Line);
        grid.add_point_scalars("temperature", vec![300.0, 350.0]);
        let output = LegacyVtkWriter::write(&grid);
        assert!(
            output.contains("temperature"),
            "expected field name 'temperature' in POINT_DATA section"
        );
    }
    #[test]
    fn test_encode_f64_le_single() {
        let v = 1.0_f64;
        let bytes = VtkBinaryWriter::encode_f64_le(&[v]);
        assert_eq!(bytes.len(), 8);
        let decoded = f64::from_le_bytes(bytes.try_into().unwrap());
        assert!((decoded - v).abs() < 1e-15);
    }
    #[test]
    fn test_encode_f64_le_multiple() {
        let vals = [1.0_f64, 2.0, 3.0];
        let bytes = VtkBinaryWriter::encode_f64_le(&vals);
        assert_eq!(bytes.len(), 24);
        let decoded: Vec<f64> = bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        for (a, b) in vals.iter().zip(decoded.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }
    #[test]
    fn test_encode_i64_le() {
        let vals = [42_i64, -1, 0];
        let bytes = VtkBinaryWriter::encode_i64_le(&vals);
        assert_eq!(bytes.len(), 24);
        let decoded: Vec<i64> = bytes
            .chunks_exact(8)
            .map(|c| i64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        assert_eq!(decoded, vals);
    }
    #[test]
    fn test_binary_data_array_xml_contains_name() {
        let xml = VtkBinaryWriter::binary_data_array_xml("pressure", 1, &[1.0, 2.0, 3.0]);
        assert!(xml.contains("Name=\"pressure\""));
        assert!(xml.contains("format=\"binary\""));
    }
    #[test]
    fn test_rle_encode_decode_roundtrip() {
        let data = vec![1.0, 1.0, 1.0, 2.0, 3.0, 3.0];
        let runs = VtkCompression::rle_encode_f64(&data);
        let decoded = VtkCompression::rle_decode_f64(&runs);
        assert_eq!(decoded.len(), data.len());
        for (a, b) in data.iter().zip(decoded.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }
    #[test]
    fn test_rle_encode_all_same() {
        let data = vec![5.0_f64; 100];
        let runs = VtkCompression::rle_encode_f64(&data);
        assert_eq!(runs.len(), 1, "All-same data should produce 1 run");
        assert_eq!(runs[0].1, 100);
    }
    #[test]
    fn test_rle_encode_empty() {
        let runs = VtkCompression::rle_encode_f64(&[]);
        assert!(runs.is_empty());
    }
    #[test]
    fn test_rle_compression_ratio() {
        let r = VtkCompression::compression_ratio(100, 10);
        assert!((r - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_delta_encode_decode_roundtrip() {
        let data = vec![10u8, 12, 11, 15, 20, 18];
        let (deltas, first) = VtkCompression::delta_encode_u8(&data);
        let decoded = VtkCompression::delta_decode_u8(&deltas, first);
        assert_eq!(decoded, data);
    }
    #[test]
    fn test_delta_encode_empty() {
        let (deltas, first) = VtkCompression::delta_encode_u8(&[]);
        assert!(deltas.is_empty());
        assert_eq!(first, 0);
    }
    #[test]
    fn test_pvtu_contains_pieces() {
        let partitions = vec![
            VtkPartition {
                rank: 0,
                filename: "p0.vtu".to_owned(),
                point_offset: 0,
                n_points: 50,
                n_cells: 10,
            },
            VtkPartition {
                rank: 1,
                filename: "p1.vtu".to_owned(),
                point_offset: 50,
                n_points: 50,
                n_cells: 10,
            },
        ];
        let pvtu = VtkParallelWriter::write_pvtu(&partitions, &[("velocity", 3)], &[]);
        assert!(pvtu.contains("<VTKFile type=\"PUnstructuredGrid\""));
        assert!(pvtu.contains("Source=\"p0.vtu\""));
        assert!(pvtu.contains("Source=\"p1.vtu\""));
        assert!(pvtu.contains("Name=\"velocity\""));
    }
    #[test]
    fn test_partition_points_balanced() {
        let parts = VtkParallelWriter::partition_points(10, 3);
        assert_eq!(parts.len(), 3);
        let total: usize = parts.iter().map(|(_, c)| c).sum();
        assert_eq!(total, 10);
        assert_eq!(parts[0].1, 4);
        assert_eq!(parts[1].1, 3);
        assert_eq!(parts[2].1, 3);
    }
    #[test]
    fn test_partition_points_zero_ranks() {
        let parts = VtkParallelWriter::partition_points(10, 0);
        assert!(parts.is_empty());
    }
    #[test]
    fn test_partition_points_zero_points() {
        let parts = VtkParallelWriter::partition_points(0, 4);
        assert!(parts.is_empty());
    }
    #[test]
    fn test_validate_vtk_grid_valid() {
        let mut grid = VtkUnstructuredGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_cell(vec![0, 1], VtkCellTypeW::Line);
        grid.add_point_scalars("s", vec![1.0, 2.0]);
        let errors = validate_vtk_grid(&grid);
        assert!(
            errors.is_empty(),
            "Valid grid should have no errors: {:?}",
            errors
        );
    }
    #[test]
    fn test_validate_vtk_grid_bad_point_data() {
        let mut grid = VtkUnstructuredGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.add_point([1.0, 0.0, 0.0]);
        grid.add_cell(vec![0, 1], VtkCellTypeW::Line);
        grid.add_point_scalars("bad", vec![1.0]);
        let errors = validate_vtk_grid(&grid);
        assert!(!errors.is_empty());
        assert!(errors[0].contains("bad"));
    }
    #[test]
    fn test_validate_vtk_grid_bad_connectivity() {
        let mut grid = VtkUnstructuredGrid::new();
        grid.add_point([0.0, 0.0, 0.0]);
        grid.cells.push(vec![0, 99]);
        grid.cell_types.push(VtkCellTypeW::Line);
        let errors = validate_vtk_grid(&grid);
        assert!(!errors.is_empty());
    }
    #[test]
    fn test_formatted_writer_precision() {
        let fw = FormattedAsciiWriter::with_precision(3);
        let s = fw.format_f64(std::f64::consts::PI);
        assert_eq!(s, "3.142");
    }
    #[test]
    fn test_formatted_writer_point() {
        let fw = FormattedAsciiWriter::with_precision(2);
        let s = fw.format_point([1.0, 2.5, 2.54321]);
        assert_eq!(s, "1.00 2.50 2.54");
    }
    #[test]
    fn test_formatted_writer_write_grid() {
        let fw = FormattedAsciiWriter::with_precision(4);
        let mut grid = VtkUnstructuredGrid::new();
        grid.add_point([1.23456789, 0.0, 0.0]);
        grid.add_cell(vec![0], VtkCellTypeW::Vertex);
        let output = fw.write_grid(&grid);
        assert!(output.contains("DATASET UNSTRUCTURED_GRID"));
        assert!(
            output.contains("1.2346") || output.contains("1.234"),
            "Precision not applied: {}",
            output
        );
    }
    #[test]
    fn test_n_points_cell_type() {
        assert_eq!(VtkCellTypeW::Vertex.n_points(), 1);
        assert_eq!(VtkCellTypeW::Line.n_points(), 2);
        assert_eq!(VtkCellTypeW::Triangle.n_points(), 3);
        assert_eq!(VtkCellTypeW::Tetrahedron.n_points(), 4);
        assert_eq!(VtkCellTypeW::Hexahedron.n_points(), 8);
    }
}
/// Write a parallel VTK XML (PVTU) master file referencing multiple piece VTU files.
pub fn write_pvtu(
    path: &str,
    pieces: &[PvtuPiece],
    point_array_names: &[&str],
    cell_array_names: &[&str],
) -> std::io::Result<()> {
    use std::io::Write;
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, r#"<?xml version="1.0"?>"#)?;
    writeln!(
        w,
        r#"<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">"#
    )?;
    writeln!(w, r#"  <PUnstructuredGrid GhostLevel="0">"#)?;
    if !point_array_names.is_empty() {
        writeln!(w, r#"    <PPointData>"#)?;
        for name in point_array_names {
            writeln!(w, r#"      <PDataArray type="Float64" Name="{name}"/>"#)?;
        }
        writeln!(w, r#"    </PPointData>"#)?;
    }
    if !cell_array_names.is_empty() {
        writeln!(w, r#"    <PCellData>"#)?;
        for name in cell_array_names {
            writeln!(w, r#"      <PDataArray type="Float64" Name="{name}"/>"#)?;
        }
        writeln!(w, r#"    </PCellData>"#)?;
    }
    writeln!(w, r#"    <PPoints>"#)?;
    writeln!(
        w,
        r#"      <PDataArray type="Float64" NumberOfComponents="3"/>"#
    )?;
    writeln!(w, r#"    </PPoints>"#)?;
    for piece in pieces {
        writeln!(w, r#"    <Piece Source="{}"/>"#, piece.filename)?;
    }
    writeln!(w, r#"  </PUnstructuredGrid>"#)?;
    writeln!(w, r#"</VTKFile>"#)?;
    w.flush()?;
    Ok(())
}
/// Encode bytes to Base64 (standard alphabet, no line breaks).
///
/// This is a pure-Rust implementation of the Base64 encoding used by
/// VTK's AppendedData section.
pub fn base64_encode(data: &[u8]) -> String {
    pub(super) const ALPHABET: &[u8] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = Vec::with_capacity(data.len().div_ceil(3) * 4);
    let mut chunks = data.chunks_exact(3);
    for chunk in chunks.by_ref() {
        let b0 = chunk[0] as u32;
        let b1 = chunk[1] as u32;
        let b2 = chunk[2] as u32;
        let combined = (b0 << 16) | (b1 << 8) | b2;
        out.push(ALPHABET[((combined >> 18) & 0x3F) as usize]);
        out.push(ALPHABET[((combined >> 12) & 0x3F) as usize]);
        out.push(ALPHABET[((combined >> 6) & 0x3F) as usize]);
        out.push(ALPHABET[(combined & 0x3F) as usize]);
    }
    let rem = chunks.remainder();
    match rem.len() {
        1 => {
            let b0 = rem[0] as u32;
            let combined = b0 << 16;
            out.push(ALPHABET[((combined >> 18) & 0x3F) as usize]);
            out.push(ALPHABET[((combined >> 12) & 0x3F) as usize]);
            out.push(b'=');
            out.push(b'=');
        }
        2 => {
            let b0 = rem[0] as u32;
            let b1 = rem[1] as u32;
            let combined = (b0 << 16) | (b1 << 8);
            out.push(ALPHABET[((combined >> 18) & 0x3F) as usize]);
            out.push(ALPHABET[((combined >> 12) & 0x3F) as usize]);
            out.push(ALPHABET[((combined >> 6) & 0x3F) as usize]);
            out.push(b'=');
        }
        _ => {}
    }
    String::from_utf8(out).unwrap_or_default()
}
/// Encode an array of f64 values to a Base64 string (little-endian bytes).
pub fn base64_encode_f64(data: &[f64]) -> String {
    let mut bytes = Vec::with_capacity(data.len() * 8);
    for &v in data {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    base64_encode(&bytes)
}
/// Map a VTK cell type integer to a human-readable name.
pub fn vtk_cell_type_name(type_id: u8) -> &'static str {
    match type_id {
        1 => "VTK_VERTEX",
        2 => "VTK_POLY_VERTEX",
        3 => "VTK_LINE",
        4 => "VTK_POLY_LINE",
        5 => "VTK_TRIANGLE",
        6 => "VTK_TRIANGLE_STRIP",
        7 => "VTK_POLYGON",
        8 => "VTK_PIXEL",
        9 => "VTK_QUAD",
        10 => "VTK_TETRA",
        11 => "VTK_VOXEL",
        12 => "VTK_HEXAHEDRON",
        13 => "VTK_WEDGE",
        14 => "VTK_PYRAMID",
        21 => "VTK_QUADRATIC_EDGE",
        22 => "VTK_QUADRATIC_TRIANGLE",
        23 => "VTK_QUADRATIC_QUAD",
        24 => "VTK_QUADRATIC_TETRA",
        25 => "VTK_QUADRATIC_HEXAHEDRON",
        _ => "VTK_UNKNOWN",
    }
}
/// Return the number of points for a given VTK cell type id.
///
/// Returns `None` for variable-length cell types or unknown types.
pub fn vtk_cell_n_points(type_id: u8) -> Option<usize> {
    match type_id {
        1 => Some(1),
        3 => Some(2),
        5 => Some(3),
        9 => Some(4),
        10 => Some(4),
        12 => Some(8),
        13 => Some(6),
        14 => Some(5),
        21 => Some(3),
        22 => Some(6),
        23 => Some(8),
        24 => Some(10),
        25 => Some(20),
        _ => None,
    }
}
/// Validate a VTU XML string content (in-memory).
///
/// Checks:
/// - XML has VTKFile root element
/// - type="UnstructuredGrid"
/// - Contains UnstructuredGrid and Piece elements
/// - NumberOfPoints and NumberOfCells are present
pub fn validate_vtu_xml(content: &str) -> VtuValidationResult {
    let mut result = VtuValidationResult::ok();
    if !content.contains("VTKFile") {
        result.add_issue("Missing VTKFile root element");
    }
    if !content.contains(r#"type="UnstructuredGrid""#) {
        result.add_issue("VTKFile is not of type UnstructuredGrid");
    }
    if !content.contains("UnstructuredGrid") {
        result.add_issue("Missing UnstructuredGrid element");
    }
    if !content.contains("Piece") {
        result.add_issue("Missing Piece element");
    }
    if !content.contains("NumberOfPoints") {
        result.add_issue("Missing NumberOfPoints attribute");
    }
    if !content.contains("NumberOfCells") {
        result.add_issue("Missing NumberOfCells attribute");
    }
    if !content.contains("Points") {
        result.add_issue("Missing Points element");
    }
    if !content.contains("Cells") {
        result.add_issue("Missing Cells element");
    }
    result
}
#[cfg(test)]
mod tests_vtk_writer_extended {
    use super::*;
    use crate::TimeStepWriter;

    use crate::VtuXmlWriter;
    use crate::vtk_writer::types::*;
    #[test]
    fn test_vtu_xml_writer_basic() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cells = vec![vec![0_usize, 1, 2]];
        let cell_types = vec![VtkCellTypeW::Triangle];
        let path = std::env::temp_dir().join("test_vtu_xml_basic.vtu");
        VtuXmlWriter::write(
            path.to_str().unwrap_or(""),
            &points,
            &cells,
            &cell_types,
            &[],
            &[],
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("UnstructuredGrid"));
        assert!(content.contains("NumberOfPoints=\"3\""));
        assert!(content.contains("NumberOfCells=\"1\""));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_vtu_xml_writer_with_point_data() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cells = vec![vec![0_usize, 1, 2]];
        let cell_types = vec![VtkCellTypeW::Triangle];
        let pressure = vec![1.0, 2.0, 3.0];
        let path = std::env::temp_dir().join("test_vtu_xml_pdata.vtu");
        VtuXmlWriter::write(
            path.to_str().unwrap_or(""),
            &points,
            &cells,
            &cell_types,
            &[("pressure", &pressure)],
            &[],
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PointData"));
        assert!(content.contains("pressure"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_vtu_xml_writer_with_cell_data() {
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cells = vec![vec![0_usize, 1, 2]];
        let cell_types = vec![VtkCellTypeW::Triangle];
        let stress = vec![5.0];
        let path = std::env::temp_dir().join("test_vtu_xml_cdata.vtu");
        VtuXmlWriter::write(
            path.to_str().unwrap_or(""),
            &points,
            &cells,
            &cell_types,
            &[],
            &[("stress", &stress)],
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("CellData"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_pvtu_basic() {
        let pieces = vec![
            PvtuPiece::new("sim_p0_000000.vtu"),
            PvtuPiece::new("sim_p1_000000.vtu"),
        ];
        let path = std::env::temp_dir().join("test_pvtu.pvtu");
        write_pvtu(
            path.to_str().unwrap_or(""),
            &pieces,
            &["velocity", "pressure"],
            &[],
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PUnstructuredGrid"));
        assert!(content.contains("sim_p0_000000.vtu"));
        assert!(content.contains("velocity"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_pvtu_empty_pieces() {
        let path = std::env::temp_dir().join("test_pvtu_empty.pvtu");
        write_pvtu(path.to_str().unwrap_or(""), &[], &[], &[]).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PUnstructuredGrid"));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_base64_encode_empty() {
        assert_eq!(base64_encode(b""), "");
    }
    #[test]
    fn test_base64_encode_man() {
        assert_eq!(base64_encode(b"Man"), "TWFu");
    }
    #[test]
    fn test_base64_encode_hello() {
        assert_eq!(base64_encode(b"hello"), "aGVsbG8=");
    }
    #[test]
    fn test_base64_encode_two_bytes() {
        assert_eq!(base64_encode(b"Ma"), "TWE=");
    }
    #[test]
    fn test_base64_encode_f64_length() {
        let data = vec![1.0_f64, 2.0, 3.0];
        let encoded = base64_encode_f64(&data);
        assert_eq!(encoded.len(), 32);
    }
    #[test]
    fn test_base64_encode_f64_deterministic() {
        let data = vec![std::f64::consts::PI];
        let e1 = base64_encode_f64(&data);
        let e2 = base64_encode_f64(&data);
        assert_eq!(e1, e2);
    }
    #[test]
    fn test_vtk_cell_type_name_known() {
        assert_eq!(vtk_cell_type_name(5), "VTK_TRIANGLE");
        assert_eq!(vtk_cell_type_name(10), "VTK_TETRA");
        assert_eq!(vtk_cell_type_name(12), "VTK_HEXAHEDRON");
    }
    #[test]
    fn test_vtk_cell_type_name_unknown() {
        assert_eq!(vtk_cell_type_name(200), "VTK_UNKNOWN");
    }
    #[test]
    fn test_vtk_cell_n_points_known() {
        assert_eq!(vtk_cell_n_points(5), Some(3));
        assert_eq!(vtk_cell_n_points(10), Some(4));
        assert_eq!(vtk_cell_n_points(12), Some(8));
        assert_eq!(vtk_cell_n_points(1), Some(1));
    }
    #[test]
    fn test_vtk_cell_n_points_unknown() {
        assert_eq!(vtk_cell_n_points(200), None);
    }
    #[test]
    fn test_timestep_writer_filenames() {
        let tmpdir = std::env::temp_dir();
        let w = TimeStepWriter::new(tmpdir.join("sim").to_str().unwrap_or(""), "flow");
        assert_eq!(w.vtu_filename(0), "flow_000000.vtu");
        assert_eq!(w.vtu_filename(42), "flow_000042.vtu");
    }
    #[test]
    fn test_timestep_writer_register() {
        let tmpdir = std::env::temp_dir();
        let mut w = TimeStepWriter::new(tmpdir.to_str().unwrap_or(""), "sim");
        w.register_step(0.0, 0);
        w.register_step(0.01, 1);
        assert_eq!(w.n_steps(), 2);
    }
    #[test]
    fn test_timestep_writer_pvd() {
        let tmpdir = std::env::temp_dir();
        let mut w = TimeStepWriter::new(tmpdir.to_str().unwrap_or(""), "flow");
        w.register_step(0.0, 0);
        w.register_step(0.01, 1);
        w.write_pvd("test_ts.pvd").unwrap();
        let pvd_path = tmpdir.join("test_ts.pvd");
        let content = std::fs::read_to_string(&pvd_path).unwrap();
        assert!(content.contains("flow_000000.vtu"));
        assert!(content.contains("flow_000001.vtu"));
        std::fs::remove_file(&pvd_path).ok();
    }
    #[test]
    fn test_validate_vtu_xml_valid() {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cells = vec![vec![0_usize, 1, 2]];
        let cell_types = vec![VtkCellTypeW::Triangle];
        let path = std::env::temp_dir().join("test_validate_vtu.vtu");
        VtuXmlWriter::write(
            path.to_str().unwrap_or(""),
            &pts,
            &cells,
            &cell_types,
            &[],
            &[],
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let result = validate_vtu_xml(&content);
        assert!(result.is_valid, "Issues: {:?}", result.issues);
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_validate_vtu_xml_missing_root() {
        let content = "<foo></foo>";
        let result = validate_vtu_xml(content);
        assert!(!result.is_valid);
        assert!(result.issues.iter().any(|i| i.contains("VTKFile")));
    }
    #[test]
    fn test_validate_vtu_xml_wrong_type() {
        let content = r#"<VTKFile type="PolyData"></VTKFile>"#;
        let result = validate_vtu_xml(content);
        assert!(!result.is_valid);
        assert!(result.issues.iter().any(|i| i.contains("UnstructuredGrid")));
    }
    #[test]
    fn test_pvtu_piece_new() {
        let p = PvtuPiece::new("sim_p0.vtu");
        assert_eq!(p.filename, "sim_p0.vtu");
    }
}
#[cfg(test)]
mod tests_vtk_writer_new {

    use crate::vtk_writer::VtkWriter;
    use crate::vtk_writer::types::*;
    fn triangle_mesh() -> (Vec<[f64; 3]>, Vec<Vec<usize>>, Vec<VtkCellTypeW>) {
        let pts = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let cells = vec![vec![0_usize, 1, 2]];
        let ct = vec![VtkCellTypeW::Triangle];
        (pts, cells, ct)
    }
    #[test]
    fn test_write_vector_field_creates_file() {
        let (pts, cells, ct) = triangle_mesh();
        let vecs = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let path = std::env::temp_dir().join("test_vtk_vec_field.vtu");
        VtkWriter::write_unstructured_vector_field(
            path.to_str().unwrap_or(""),
            &pts,
            &cells,
            &ct,
            "velocity",
            &vecs,
        )
        .unwrap();
        assert!(path.exists());
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_vector_field_contains_field_name() {
        let (pts, cells, ct) = triangle_mesh();
        let vecs = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let path = std::env::temp_dir().join("test_vtk_vec_name.vtu");
        VtkWriter::write_unstructured_vector_field(
            path.to_str().unwrap_or(""),
            &pts,
            &cells,
            &ct,
            "myVelocity",
            &vecs,
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("myVelocity"), "Field name not in file");
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_vector_field_correct_point_count() {
        let (pts, cells, ct) = triangle_mesh();
        let vecs = vec![[0.0, 0.0, 0.0]; 3];
        let path = std::env::temp_dir().join("test_vtk_vec_pts.vtu");
        VtkWriter::write_unstructured_vector_field(
            path.to_str().unwrap_or(""),
            &pts,
            &cells,
            &ct,
            "v",
            &vecs,
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("NumberOfPoints=\"3\""));
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_streamlines_creates_file() {
        let sl = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 2.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 2.0, 0.0]],
        ];
        let path = std::env::temp_dir().join("test_streamlines.vtk");
        VtkWriter::write_streamlines(path.to_str().unwrap_or(""), &sl).unwrap();
        assert!(path.exists());
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_streamlines_contains_polydata() {
        let sl = vec![vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]];
        let path = std::env::temp_dir().join("test_streamlines_pd.vtk");
        VtkWriter::write_streamlines(path.to_str().unwrap_or(""), &sl).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("POLYDATA"), "Expected POLYDATA keyword");
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_streamlines_correct_point_count() {
        let sl = vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let path = std::env::temp_dir().join("test_streamlines_npts.vtk");
        VtkWriter::write_streamlines(path.to_str().unwrap_or(""), &sl).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("POINTS 5 double"), "Expected 5 points");
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_rectilinear_grid_creates_file() {
        let xs: Vec<f64> = (0..4).map(|i| i as f64).collect();
        let ys: Vec<f64> = (0..3).map(|i| i as f64).collect();
        let zs = vec![0.0_f64];
        let path = std::env::temp_dir().join("test_rect_grid.vtk");
        VtkWriter::write_rectilinear_grid(path.to_str().unwrap_or(""), &xs, &ys, &zs, None, &[])
            .unwrap();
        assert!(path.exists());
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_rectilinear_grid_dimensions_in_file() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0];
        let zs = vec![0.0];
        let path = std::env::temp_dir().join("test_rect_grid_dim.vtk");
        VtkWriter::write_rectilinear_grid(path.to_str().unwrap_or(""), &xs, &ys, &zs, None, &[])
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(
            content.contains("DIMENSIONS 3 2 1"),
            "Dimensions not correct"
        );
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_rectilinear_grid_with_scalars() {
        let xs = vec![0.0, 1.0];
        let ys = vec![0.0, 1.0];
        let zs = vec![0.0];
        let scalars = vec![1.0_f64, 2.0, 3.0, 4.0];
        let path = std::env::temp_dir().join("test_rect_grid_scalar.vtk");
        VtkWriter::write_rectilinear_grid(
            path.to_str().unwrap_or(""),
            &xs,
            &ys,
            &zs,
            Some("pressure"),
            &scalars,
        )
        .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("pressure"), "Scalar field name missing");
        assert!(content.contains("POINT_DATA 4"), "POINT_DATA count wrong");
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_write_rectilinear_grid_rectilinear_keyword() {
        let xs = vec![0.0_f64; 2];
        let ys = vec![0.0_f64; 2];
        let zs = vec![0.0_f64; 1];
        let path = std::env::temp_dir().join("test_rect_grid_kw.vtk");
        VtkWriter::write_rectilinear_grid(path.to_str().unwrap_or(""), &xs, &ys, &zs, None, &[])
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(
            content.contains("RECTILINEAR_GRID"),
            "Missing RECTILINEAR_GRID keyword"
        );
        std::fs::remove_file(&path).ok();
    }
}
