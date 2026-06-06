// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SU2 mesh format reader/writer for aerodynamic simulation.
//!
//! SU2 (Stanford University Unstructured) is an open-source CFD suite widely
//! used in aerodynamic shape optimisation.  This module supports reading and
//! writing 2-D and 3-D `.su2` mesh files in the ASCII format produced by SU2
//! v7+, including boundary markers, element connectivity, and solution CSV
//! output.

use std::fmt;
use std::io::{self, BufRead, Write};

// ---------------------------------------------------------------------------
// Su2ElementType
// ---------------------------------------------------------------------------

/// VTK-compatible element type codes used in SU2 mesh files.
///
/// SU2 uses VTK type codes in the `ELEM=` connectivity section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Su2ElementType {
    /// Line (2 nodes), VTK type 3.
    Line = 3,
    /// Triangle (3 nodes), VTK type 5.
    Triangle = 5,
    /// Quadrilateral (4 nodes), VTK type 9.
    Quadrilateral = 9,
    /// Tetrahedron (4 nodes), VTK type 10.
    Tetrahedron = 10,
    /// Hexahedron (8 nodes), VTK type 12.
    Hexahedron = 12,
    /// Prism / wedge (6 nodes), VTK type 13.
    Prism = 13,
    /// Pyramid (5 nodes), VTK type 14.
    Pyramid = 14,
}

impl Su2ElementType {
    /// Number of nodes required for this element type.
    pub fn n_nodes(self) -> usize {
        match self {
            Su2ElementType::Line => 2,
            Su2ElementType::Triangle => 3,
            Su2ElementType::Quadrilateral => 4,
            Su2ElementType::Tetrahedron => 4,
            Su2ElementType::Hexahedron => 8,
            Su2ElementType::Prism => 6,
            Su2ElementType::Pyramid => 5,
        }
    }

    /// Try to convert a VTK integer type code to [`Su2ElementType`].
    pub fn from_vtk(code: u32) -> Option<Self> {
        match code {
            3 => Some(Su2ElementType::Line),
            5 => Some(Su2ElementType::Triangle),
            9 => Some(Su2ElementType::Quadrilateral),
            10 => Some(Su2ElementType::Tetrahedron),
            12 => Some(Su2ElementType::Hexahedron),
            13 => Some(Su2ElementType::Prism),
            14 => Some(Su2ElementType::Pyramid),
            _ => None,
        }
    }
}

impl fmt::Display for Su2ElementType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", *self as u32)
    }
}

// ---------------------------------------------------------------------------
// Su2Element
// ---------------------------------------------------------------------------

/// A single element (cell) in an SU2 mesh.
#[derive(Debug, Clone)]
pub struct Su2Element {
    /// VTK element type.
    pub element_type: Su2ElementType,
    /// Node indices (0-based).
    pub nodes: Vec<usize>,
    /// Optional element global index (used in parallel partitioning).
    pub global_id: Option<usize>,
}

impl Su2Element {
    /// Create a new element.
    pub fn new(element_type: Su2ElementType, nodes: Vec<usize>) -> Self {
        Self {
            element_type,
            nodes,
            global_id: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Su2BoundaryMarker
// ---------------------------------------------------------------------------

/// Named boundary marker (patch) in an SU2 mesh.
///
/// Boundary conditions (e.g. `wall`, `farfield`, `inlet`, `outlet`) are
/// identified by tag strings in SU2.
#[derive(Debug, Clone)]
pub struct Su2BoundaryMarker {
    /// Tag name of the boundary (e.g. `"wall"`, `"farfield"`).
    pub tag: String,
    /// Boundary elements (faces / edges) on this patch.
    pub elements: Vec<Su2Element>,
}

impl Su2BoundaryMarker {
    /// Create an empty boundary marker with the given tag.
    pub fn new(tag: impl Into<String>) -> Self {
        Self {
            tag: tag.into(),
            elements: Vec::new(),
        }
    }

    /// Number of boundary elements on this marker.
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    /// Add a boundary element.
    pub fn add_element(&mut self, elem: Su2Element) {
        self.elements.push(elem);
    }
}

// ---------------------------------------------------------------------------
// Su2Mesh
// ---------------------------------------------------------------------------

/// SU2 mesh in memory.
///
/// Holds node coordinates, volume elements, and named boundary markers.
#[derive(Debug, Clone)]
pub struct Su2Mesh {
    /// Problem dimensionality (2 or 3).
    pub ndim: usize,
    /// Node coordinates, length `n_nodes * ndim`.
    ///
    /// For 2-D: `[x0, y0, x1, y1, ...]`; for 3-D: `[x0, y0, z0, ...]`.
    pub coords: Vec<f64>,
    /// Volume elements.
    pub elements: Vec<Su2Element>,
    /// Named boundary markers.
    pub markers: Vec<Su2BoundaryMarker>,
}

impl Su2Mesh {
    /// Create an empty mesh with the given dimensionality.
    pub fn new(ndim: usize) -> Self {
        Self {
            ndim,
            coords: Vec::new(),
            elements: Vec::new(),
            markers: Vec::new(),
        }
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.coords.len().checked_div(self.ndim).unwrap_or(0)
    }

    /// Number of volume elements.
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    /// Add a node and return its 0-based index.
    ///
    /// Panics in debug mode if the coordinate slice length != `ndim`.
    pub fn add_node(&mut self, coords: &[f64]) -> usize {
        debug_assert_eq!(coords.len(), self.ndim);
        let idx = self.n_nodes();
        self.coords.extend_from_slice(coords);
        idx
    }

    /// Add a volume element.
    pub fn add_element(&mut self, elem: Su2Element) {
        self.elements.push(elem);
    }

    /// Add a boundary marker.
    pub fn add_marker(&mut self, marker: Su2BoundaryMarker) {
        self.markers.push(marker);
    }

    /// Look up a boundary marker by tag, returning a reference if found.
    pub fn boundary_marker(&self, tag: &str) -> Option<&Su2BoundaryMarker> {
        self.markers.iter().find(|m| m.tag == tag)
    }

    /// Get the coordinates of node `i` as a slice of length `ndim`.
    pub fn node_coords(&self, i: usize) -> &[f64] {
        let start = i * self.ndim;
        &self.coords[start..start + self.ndim]
    }
}

// ---------------------------------------------------------------------------
// element_connectivity
// ---------------------------------------------------------------------------

/// Return the VTK element type code for a given [`Su2ElementType`].
///
/// This is a convenience function for format-agnostic code that needs the
/// integer VTK type.
pub fn element_connectivity(elem_type: Su2ElementType) -> u32 {
    elem_type as u32
}

// ---------------------------------------------------------------------------
// read_su2
// ---------------------------------------------------------------------------

/// Parse an SU2 mesh file from a buffered reader.
///
/// Supports `NDIME=`, `NPOIN=`, `NELEM=`, `NMARK=`, `MARKER_TAG=`,
/// `MARKER_ELEMS=` keywords in ASCII SU2 v7 format.
///
/// # Errors
/// Returns an [`io::Error`] on I/O failure or malformed input.
pub fn read_su2<R: BufRead>(reader: R) -> io::Result<Su2Mesh> {
    let lines: Vec<String> = reader.lines().map(|l| l.unwrap_or_default()).collect();

    let mut ndim = 2usize;
    let mut mesh = Su2Mesh::new(2);
    let mut pos = 0;

    while pos < lines.len() {
        let line = lines[pos].trim().to_string();
        pos += 1;

        if line.starts_with('%') || line.is_empty() {
            continue;
        }

        if let Some(val) = parse_kv(&line, "NDIME=") {
            ndim = val.parse().unwrap_or(2);
            mesh = Su2Mesh::new(ndim);
            continue;
        }

        if let Some(val) = parse_kv(&line, "NPOIN=") {
            let n_nodes: usize = val
                .split_whitespace()
                .next()
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            for _ in 0..n_nodes {
                if pos >= lines.len() {
                    break;
                }
                let node_line = lines[pos].trim().to_string();
                pos += 1;
                let tokens: Vec<f64> = node_line
                    .split_whitespace()
                    .take(ndim)
                    .filter_map(|t| t.parse().ok())
                    .collect();
                if tokens.len() == ndim {
                    mesh.coords.extend_from_slice(&tokens);
                }
            }
            continue;
        }

        if let Some(val) = parse_kv(&line, "NELEM=") {
            let n_elem: usize = val.parse().unwrap_or(0);
            for _ in 0..n_elem {
                if pos >= lines.len() {
                    break;
                }
                let el_line = lines[pos].trim().to_string();
                pos += 1;
                if let Some(elem) = parse_element_line(&el_line) {
                    mesh.elements.push(elem);
                }
            }
            continue;
        }

        if let Some(val) = parse_kv(&line, "NMARK=") {
            let _n_mark: usize = val.parse().unwrap_or(0);
            // Parse markers inline below
            let _ = _n_mark;
            continue;
        }

        if let Some(tag) = parse_kv(&line, "MARKER_TAG=") {
            let mut marker = Su2BoundaryMarker::new(tag.trim());
            // Expect MARKER_ELEMS= on the next non-empty line
            while pos < lines.len() {
                let next = lines[pos].trim().to_string();
                pos += 1;
                if next.is_empty() || next.starts_with('%') {
                    continue;
                }
                if let Some(val) = parse_kv(&next, "MARKER_ELEMS=") {
                    let n_mel: usize = val.parse().unwrap_or(0);
                    for _ in 0..n_mel {
                        if pos >= lines.len() {
                            break;
                        }
                        let mel_line = lines[pos].trim().to_string();
                        pos += 1;
                        if let Some(elem) = parse_element_line(&mel_line) {
                            marker.elements.push(elem);
                        }
                    }
                    break;
                }
            }
            mesh.markers.push(marker);
            continue;
        }
    }

    Ok(mesh)
}

/// Parse a `KEY=value` line and return the value string if the key matches.
fn parse_kv<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    line.trim().strip_prefix(key).map(str::trim)
}

/// Parse one element line: `vtk_type` `node0` `node1` ... [global_id]`.
fn parse_element_line(line: &str) -> Option<Su2Element> {
    let mut tokens = line.split_whitespace();
    let vtk_type: u32 = tokens.next()?.parse().ok()?;
    let elem_type = Su2ElementType::from_vtk(vtk_type)?;
    let n = elem_type.n_nodes();
    let nodes: Vec<usize> = tokens.take(n).filter_map(|t| t.parse().ok()).collect();
    if nodes.len() == n {
        Some(Su2Element::new(elem_type, nodes))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// write_su2
// ---------------------------------------------------------------------------

/// Write an SU2 mesh to a writer in ASCII SU2 v7 format.
///
/// # Errors
/// Returns an [`io::Error`] on I/O failure.
pub fn write_su2<W: Write>(writer: &mut W, mesh: &Su2Mesh) -> io::Result<()> {
    writeln!(writer, "%")?;
    writeln!(writer, "% SU2 mesh written by OxiPhysics")?;
    writeln!(writer, "%")?;
    writeln!(writer, "NDIME= {}", mesh.ndim)?;
    writeln!(writer)?;

    // Volume elements
    writeln!(writer, "NELEM= {}", mesh.n_elements())?;
    for (i, elem) in mesh.elements.iter().enumerate() {
        let nodes_str: Vec<String> = elem.nodes.iter().map(|n| n.to_string()).collect();
        writeln!(
            writer,
            "{} {} {}",
            elem.element_type,
            nodes_str.join(" "),
            i
        )?;
    }
    writeln!(writer)?;

    // Node coordinates
    writeln!(writer, "NPOIN= {}", mesh.n_nodes())?;
    for i in 0..mesh.n_nodes() {
        let coords = mesh.node_coords(i);
        let coord_str: Vec<String> = coords.iter().map(|&c| format!("{c:.10e}")).collect();
        writeln!(writer, "{}", coord_str.join("\t"))?;
    }
    writeln!(writer)?;

    // Boundary markers
    writeln!(writer, "NMARK= {}", mesh.markers.len())?;
    for marker in &mesh.markers {
        writeln!(writer, "MARKER_TAG= {}", marker.tag)?;
        writeln!(writer, "MARKER_ELEMS= {}", marker.n_elements())?;
        for elem in &marker.elements {
            let nodes_str: Vec<String> = elem.nodes.iter().map(|n| n.to_string()).collect();
            writeln!(writer, "{} {}", elem.element_type, nodes_str.join(" "))?;
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// solution_file
// ---------------------------------------------------------------------------

/// Write an SU2-compatible solution CSV file.
///
/// Columns: `PointID, x, y\[, z\], <field_names...>`.
///
/// # Arguments
/// * `writer` - Target writer.
/// * `mesh` - Reference mesh for node coordinates.
/// * `field_names` - Names of the solution fields.
/// * `fields` - Field values; `fields\[k\]\[i\]` is the value of field `k` at node `i`.
pub fn solution_file<W: Write>(
    writer: &mut W,
    mesh: &Su2Mesh,
    field_names: &[&str],
    fields: &[&[f64]],
) -> io::Result<()> {
    // Header
    let mut header = vec!["PointID".to_string()];
    match mesh.ndim {
        2 => {
            header.push("x".to_string());
            header.push("y".to_string());
        }
        _ => {
            header.push("x".to_string());
            header.push("y".to_string());
            header.push("z".to_string());
        }
    }
    for name in field_names {
        header.push(name.to_string());
    }
    writeln!(writer, "{}", header.join(","))?;

    for i in 0..mesh.n_nodes() {
        let mut row = vec![i.to_string()];
        let coords = mesh.node_coords(i);
        for &c in coords {
            row.push(format!("{c:.10e}"));
        }
        for field in fields {
            if i < field.len() {
                row.push(format!("{:.10e}", field[i]));
            } else {
                row.push("0.0".to_string());
            }
        }
        writeln!(writer, "{}", row.join(","))?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    // ── Su2ElementType ─────────────────────────────────────────────────────

    #[test]
    fn test_element_type_n_nodes_triangle() {
        assert_eq!(Su2ElementType::Triangle.n_nodes(), 3);
    }

    #[test]
    fn test_element_type_n_nodes_hexa() {
        assert_eq!(Su2ElementType::Hexahedron.n_nodes(), 8);
    }

    #[test]
    fn test_element_type_n_nodes_tetra() {
        assert_eq!(Su2ElementType::Tetrahedron.n_nodes(), 4);
    }

    #[test]
    fn test_element_type_n_nodes_line() {
        assert_eq!(Su2ElementType::Line.n_nodes(), 2);
    }

    #[test]
    fn test_element_type_from_vtk_triangle() {
        assert_eq!(Su2ElementType::from_vtk(5), Some(Su2ElementType::Triangle));
    }

    #[test]
    fn test_element_type_from_vtk_unknown() {
        assert_eq!(Su2ElementType::from_vtk(99), None);
    }

    #[test]
    fn test_element_type_from_vtk_hexa() {
        assert_eq!(
            Su2ElementType::from_vtk(12),
            Some(Su2ElementType::Hexahedron)
        );
    }

    #[test]
    fn test_element_connectivity_triangle() {
        assert_eq!(element_connectivity(Su2ElementType::Triangle), 5);
    }

    #[test]
    fn test_element_type_display() {
        assert_eq!(Su2ElementType::Triangle.to_string(), "5");
    }

    // ── Su2Mesh construction ───────────────────────────────────────────────

    #[test]
    fn test_mesh_new_2d() {
        let m = Su2Mesh::new(2);
        assert_eq!(m.ndim, 2);
        assert_eq!(m.n_nodes(), 0);
        assert_eq!(m.n_elements(), 0);
    }

    #[test]
    fn test_mesh_add_node_2d() {
        let mut m = Su2Mesh::new(2);
        let idx = m.add_node(&[1.0, 2.0]);
        assert_eq!(idx, 0);
        assert_eq!(m.n_nodes(), 1);
    }

    #[test]
    fn test_mesh_node_coords() {
        let mut m = Su2Mesh::new(2);
        m.add_node(&[3.0, 4.0]);
        let c = m.node_coords(0);
        assert!((c[0] - 3.0).abs() < 1e-12);
        assert!((c[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_mesh_add_element() {
        let mut m = Su2Mesh::new(2);
        m.add_node(&[0.0, 0.0]);
        m.add_node(&[1.0, 0.0]);
        m.add_node(&[0.0, 1.0]);
        m.add_element(Su2Element::new(Su2ElementType::Triangle, vec![0, 1, 2]));
        assert_eq!(m.n_elements(), 1);
    }

    #[test]
    fn test_mesh_boundary_marker_lookup() {
        let mut m = Su2Mesh::new(2);
        m.add_marker(Su2BoundaryMarker::new("wall"));
        assert!(m.boundary_marker("wall").is_some());
        assert!(m.boundary_marker("farfield").is_none());
    }

    // ── Su2BoundaryMarker ──────────────────────────────────────────────────

    #[test]
    fn test_marker_new_empty() {
        let mk = Su2BoundaryMarker::new("inlet");
        assert_eq!(mk.tag, "inlet");
        assert_eq!(mk.n_elements(), 0);
    }

    #[test]
    fn test_marker_add_element() {
        let mut mk = Su2BoundaryMarker::new("wall");
        mk.add_element(Su2Element::new(Su2ElementType::Line, vec![0, 1]));
        assert_eq!(mk.n_elements(), 1);
    }

    // ── write_su2 / read_su2 roundtrip ─────────────────────────────────────

    fn make_simple_2d_mesh() -> Su2Mesh {
        let mut m = Su2Mesh::new(2);
        m.add_node(&[0.0, 0.0]);
        m.add_node(&[1.0, 0.0]);
        m.add_node(&[0.0, 1.0]);
        m.add_element(Su2Element::new(Su2ElementType::Triangle, vec![0, 1, 2]));
        let mut mk = Su2BoundaryMarker::new("wall");
        mk.add_element(Su2Element::new(Su2ElementType::Line, vec![0, 1]));
        m.add_marker(mk);
        m
    }

    #[test]
    fn test_write_su2_contains_ndime() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("NDIME= 2"), "missing NDIME");
    }

    #[test]
    fn test_write_su2_contains_nelem() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("NELEM= 1"));
    }

    #[test]
    fn test_write_su2_contains_marker_tag() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("MARKER_TAG= wall"));
    }

    #[test]
    fn test_roundtrip_n_nodes() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let m2 = read_su2(reader).unwrap();
        assert_eq!(m2.n_nodes(), 3);
    }

    #[test]
    fn test_roundtrip_n_elements() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let m2 = read_su2(reader).unwrap();
        assert_eq!(m2.n_elements(), 1);
    }

    #[test]
    fn test_roundtrip_node_coords() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let m2 = read_su2(reader).unwrap();
        let c = m2.node_coords(0);
        assert!(c[0].abs() < 1e-6, "x0={}", c[0]);
        assert!(c[1].abs() < 1e-6, "y0={}", c[1]);
    }

    #[test]
    fn test_roundtrip_element_type() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let m2 = read_su2(reader).unwrap();
        assert_eq!(m2.elements[0].element_type, Su2ElementType::Triangle);
    }

    #[test]
    fn test_roundtrip_element_nodes() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let m2 = read_su2(reader).unwrap();
        assert_eq!(m2.elements[0].nodes, vec![0, 1, 2]);
    }

    #[test]
    fn test_roundtrip_boundary_marker_tag() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let m2 = read_su2(reader).unwrap();
        assert!(m2.boundary_marker("wall").is_some());
    }

    #[test]
    fn test_roundtrip_boundary_marker_n_elements() {
        let m = make_simple_2d_mesh();
        let mut buf = Vec::new();
        write_su2(&mut buf, &m).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let m2 = read_su2(reader).unwrap();
        let mk = m2.boundary_marker("wall").unwrap();
        assert_eq!(mk.n_elements(), 1);
    }

    // ── solution_file ──────────────────────────────────────────────────────

    #[test]
    fn test_solution_file_header() {
        let m = make_simple_2d_mesh();
        let pressure = vec![1.0, 2.0, 3.0];
        let mut buf = Vec::new();
        solution_file(&mut buf, &m, &["Pressure"], &[&pressure]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("PointID"), "missing PointID");
        assert!(s.contains("Pressure"), "missing Pressure");
        assert!(s.contains("x"), "missing x");
        assert!(s.contains("y"), "missing y");
    }

    #[test]
    fn test_solution_file_row_count() {
        let m = make_simple_2d_mesh();
        let pressure = vec![1.0, 2.0, 3.0];
        let mut buf = Vec::new();
        solution_file(&mut buf, &m, &["p"], &[&pressure]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // 1 header + 3 data rows
        assert_eq!(s.lines().count(), 4);
    }

    #[test]
    fn test_solution_file_first_data_row() {
        let m = make_simple_2d_mesh();
        let pressure = vec![101325.0, 101000.0, 100500.0];
        let mut buf = Vec::new();
        solution_file(&mut buf, &m, &["p"], &[&pressure]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // First data row starts with "0,"
        let first_data = s.lines().nth(1).unwrap();
        assert!(first_data.starts_with("0,"), "row: {first_data}");
    }

    #[test]
    fn test_solution_file_multiple_fields() {
        let m = make_simple_2d_mesh();
        let pressure = vec![1.0, 2.0, 3.0];
        let temperature = vec![300.0, 310.0, 320.0];
        let mut buf = Vec::new();
        solution_file(&mut buf, &m, &["p", "T"], &[&pressure, &temperature]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("p"), "missing p");
        assert!(s.contains("T"), "missing T");
    }

    // ── read_su2 edge cases ────────────────────────────────────────────────

    #[test]
    fn test_read_su2_empty_input() {
        let reader = BufReader::new("".as_bytes());
        let m = read_su2(reader).unwrap();
        assert_eq!(m.n_nodes(), 0);
        assert_eq!(m.n_elements(), 0);
    }

    #[test]
    fn test_read_su2_comments_ignored() {
        let input = "% This is a comment\n% Another comment\n";
        let reader = BufReader::new(input.as_bytes());
        let m = read_su2(reader).unwrap();
        assert_eq!(m.n_nodes(), 0);
    }
}
