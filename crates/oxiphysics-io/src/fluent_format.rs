// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! ANSYS Fluent mesh format I/O (ASCII `.msh`).
//!
//! Supports reading and writing Fluent mesh files containing nodes, faces,
//! cells, and zone assignments.  The implementation covers the major section
//! types used in Fluent ASCII mesh files:
//!
//! * `(0 …)` — comment
//! * `(2 …)` — dimension
//! * `(10 …)` — node section
//! * `(12 …)` — cell section
//! * `(13 …)` — face section
//! * `(39 …)`, `(45 …)` — zone section
//!
//! # Reading
//!
//! [`FluentReader`] reads the layout Fluent itself writes, for example
//!
//! ```text
//! (2 3)
//! (10 (0 1 2d4 0 3))
//! (10 (1 1 2d4 1 3)(
//! 0.0 0.0 0.0
//! ...
//! ))
//! (13 (4 1 b 3 0)(
//! 3 1 2 3 1 0
//! ...
//! ))
//! ```
//!
//! The reader follows the parenthesis structure, not the lines, so the data
//! of a section may start on the header line or on a line of its own and may
//! be broken at any white space.  Header numbers are hexadecimal, each node
//! zone is read with the coordinate count of its header (two-dimensional
//! meshes have two), nodes and faces take their ids from the index range of
//! their zone, and zones of face type `0` (mixed) or `5` (polygonal) are
//! read with the node count in front of every face.
//!
//! Binary meshes (sections 2010, 3010, 2012, ...) and the content
//! [`FluentMesh`] has no place for (faces with more than four nodes,
//! polyhedral cells) are rejected with an error rather than being misread.

use std::fmt::Write as FmtWrite;
use std::fs;
use std::io::{self, BufRead};

use crate::Error as IoError;

#[cfg(test)]
mod reader_tests;
mod scan;
mod sections;

// ── FluentZoneType ────────────────────────────────────────────────────────────

/// Fluent zone type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FluentZoneType {
    /// Fluid volume zone.
    Fluid,
    /// Solid volume zone.
    Solid,
    /// Wall boundary.
    Wall,
    /// Inlet boundary.
    Inlet,
    /// Outlet boundary.
    Outlet,
    /// Symmetry boundary.
    Symmetry,
    /// Interior (internal) face zone.
    Interior,
}

impl FluentZoneType {
    /// Fluent keyword string for this zone type.
    pub fn as_str(&self) -> &str {
        match self {
            FluentZoneType::Fluid => "fluid",
            FluentZoneType::Solid => "solid",
            FluentZoneType::Wall => "wall",
            FluentZoneType::Inlet => "velocity-inlet",
            FluentZoneType::Outlet => "pressure-outlet",
            FluentZoneType::Symmetry => "symmetry",
            FluentZoneType::Interior => "interior",
        }
    }

    /// Parse a Fluent keyword string into a zone type.
    pub fn from_keyword(s: &str) -> Self {
        Self::from(s)
    }
}

impl From<&str> for FluentZoneType {
    fn from(s: &str) -> Self {
        match s.trim() {
            "fluid" => FluentZoneType::Fluid,
            "solid" => FluentZoneType::Solid,
            "wall" => FluentZoneType::Wall,
            "velocity-inlet" | "inlet" => FluentZoneType::Inlet,
            "pressure-outlet" | "outlet" => FluentZoneType::Outlet,
            "symmetry" => FluentZoneType::Symmetry,
            "interior" => FluentZoneType::Interior,
            _ => FluentZoneType::Interior,
        }
    }
}

// ── FluentNode ────────────────────────────────────────────────────────────────

/// A mesh node (vertex) in a Fluent mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct FluentNode {
    /// 1-based node identifier.
    pub id: usize,
    /// Node coordinates `[x, y, z]`.
    pub coordinates: [f64; 3],
}

impl FluentNode {
    /// Create a new node.
    pub fn new(id: usize, coordinates: [f64; 3]) -> Self {
        Self { id, coordinates }
    }
}

// ── FluentFaceType ────────────────────────────────────────────────────────────

/// Face element topology type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FluentFaceType {
    /// 2-node line (2-D edge).
    Line,
    /// 3-node triangle.
    Triangle,
    /// 4-node quadrilateral.
    Quad,
}

impl FluentFaceType {
    /// Fluent face type integer code.
    pub fn code(&self) -> u32 {
        match self {
            FluentFaceType::Line => 2,
            FluentFaceType::Triangle => 3,
            FluentFaceType::Quad => 4,
        }
    }

    /// Parse a Fluent face type code.
    pub fn from_code(code: u32) -> Self {
        match code {
            2 => FluentFaceType::Line,
            3 => FluentFaceType::Triangle,
            4 => FluentFaceType::Quad,
            _ => FluentFaceType::Triangle,
        }
    }

    /// Number of nodes for this face type.
    pub fn node_count(&self) -> usize {
        match self {
            FluentFaceType::Line => 2,
            FluentFaceType::Triangle => 3,
            FluentFaceType::Quad => 4,
        }
    }

    /// The face type of a face with `count` nodes, or `None` when there is
    /// none (anything but 2, 3 or 4 nodes).
    pub fn from_node_count(count: usize) -> Option<Self> {
        match count {
            2 => Some(FluentFaceType::Line),
            3 => Some(FluentFaceType::Triangle),
            4 => Some(FluentFaceType::Quad),
            _ => None,
        }
    }
}

// ── FluentFace ────────────────────────────────────────────────────────────────

/// A mesh face (edge in 2-D, polygon in 3-D) in a Fluent mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct FluentFace {
    /// 1-based face identifier.
    pub id: usize,
    /// Face topology type.
    pub face_type: FluentFaceType,
    /// Indices of the nodes forming this face (1-based).
    pub node_ids: Vec<usize>,
    /// Left (owner) cell id, or 0 for boundary faces.
    pub left_cell: usize,
    /// Right (neighbour) cell id, or 0 for boundary faces.
    pub right_cell: usize,
}

impl FluentFace {
    /// Create a new face.
    pub fn new(
        id: usize,
        face_type: FluentFaceType,
        node_ids: Vec<usize>,
        left_cell: usize,
        right_cell: usize,
    ) -> Self {
        Self {
            id,
            face_type,
            node_ids,
            left_cell,
            right_cell,
        }
    }

    /// True if this is a boundary face (one adjacent cell is 0).
    pub fn is_boundary(&self) -> bool {
        self.left_cell == 0 || self.right_cell == 0
    }
}

// ── FluentCellType ────────────────────────────────────────────────────────────

/// Cell element topology type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FluentCellType {
    /// 2-D triangle.
    Tri,
    /// 3-D tetrahedron.
    Tet,
    /// 2-D quadrilateral.
    Quad,
    /// 3-D hexahedron.
    Hex,
    /// 3-D pyramid (quadrilateral base, one apex).
    Pyramid,
    /// 3-D wedge (triangular prism).
    Wedge,
}

impl FluentCellType {
    /// Fluent cell type integer code.
    pub fn code(&self) -> u32 {
        match self {
            FluentCellType::Tri => 1,
            FluentCellType::Tet => 2,
            FluentCellType::Quad => 3,
            FluentCellType::Hex => 4,
            FluentCellType::Pyramid => 5,
            FluentCellType::Wedge => 6,
        }
    }

    /// The cell type of a Fluent cell type code, or `None` for a code this
    /// type has no variant for (`0` mixed, `7` polyhedral, anything unknown).
    pub fn try_from_code(code: u32) -> Option<Self> {
        match code {
            1 => Some(FluentCellType::Tri),
            2 => Some(FluentCellType::Tet),
            3 => Some(FluentCellType::Quad),
            4 => Some(FluentCellType::Hex),
            5 => Some(FluentCellType::Pyramid),
            6 => Some(FluentCellType::Wedge),
            _ => None,
        }
    }

    /// Parse a Fluent cell type code, falling back to [`FluentCellType::Tet`]
    /// for codes that have no variant.  Prefer
    /// [`FluentCellType::try_from_code`], which reports them.
    pub fn from_code(code: u32) -> Self {
        Self::try_from_code(code).unwrap_or(FluentCellType::Tet)
    }
}

// ── FluentCell ────────────────────────────────────────────────────────────────

/// A mesh cell (volume element) in a Fluent mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct FluentCell {
    /// 1-based cell identifier.
    pub id: usize,
    /// Element topology type.
    pub cell_type: FluentCellType,
    /// Zone this cell belongs to (zone id).
    pub zone_id: usize,
}

impl FluentCell {
    /// Create a new cell.
    pub fn new(id: usize, cell_type: FluentCellType, zone_id: usize) -> Self {
        Self {
            id,
            cell_type,
            zone_id,
        }
    }
}

// ── FluentMesh ────────────────────────────────────────────────────────────────

/// A complete Fluent mesh.
#[derive(Debug, Clone, Default)]
pub struct FluentMesh {
    /// All nodes.
    pub nodes: Vec<FluentNode>,
    /// All faces.
    pub faces: Vec<FluentFace>,
    /// All cells.
    pub cells: Vec<FluentCell>,
    /// Zone list: `(zone_id, zone_type)`.
    pub zones: Vec<(usize, FluentZoneType)>,
}

impl FluentMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a zone to the mesh.
    pub fn add_zone(&mut self, zone_id: usize, zone_type: FluentZoneType) {
        self.zones.push((zone_id, zone_type));
    }

    /// Write the mesh to a Fluent ASCII `.msh` file at `path`.
    pub fn write(&self, path: &str) -> crate::Result<()> {
        let writer = FluentWriter::new(self);
        writer.write_to_path(path)
    }

    /// Read a Fluent ASCII `.msh` file from `path`.
    pub fn read(path: &str) -> crate::Result<Self> {
        FluentReader::read_from_path(path)
    }

    /// Number of nodes in the mesh.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of faces in the mesh.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Number of cells in the mesh.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Look up a zone type by id.
    pub fn zone_type(&self, zone_id: usize) -> Option<&FluentZoneType> {
        self.zones
            .iter()
            .find(|(id, _)| *id == zone_id)
            .map(|(_, t)| t)
    }
}

// ── FluentWriter ─────────────────────────────────────────────────────────────

/// Writes a `FluentMesh` to an ASCII Fluent `.msh` file.
pub struct FluentWriter<'a> {
    mesh: &'a FluentMesh,
}

impl<'a> FluentWriter<'a> {
    /// Create a writer for the given mesh.
    pub fn new(mesh: &'a FluentMesh) -> Self {
        Self { mesh }
    }

    /// Serialize the mesh to a `String`.
    pub fn to_string(&self) -> crate::Result<String> {
        let mut buf = String::new();
        self.write_comment(&mut buf)?;
        self.write_dimension(&mut buf)?;
        self.write_nodes(&mut buf)?;
        self.write_cells(&mut buf)?;
        self.write_faces(&mut buf)?;
        self.write_zones(&mut buf)?;
        Ok(buf)
    }

    /// Write the mesh to a file.
    pub fn write_to_path(&self, path: &str) -> crate::Result<()> {
        let content = self.to_string()?;
        fs::write(path, content)?;
        Ok(())
    }

    fn write_comment(&self, buf: &mut String) -> crate::Result<()> {
        writeln!(buf, "(0 \"OxiPhysics Fluent mesh export\")")?;
        Ok(())
    }

    /// Coordinates per node: 3 when any node has a z coordinate, else 2.
    fn dimension(&self) -> usize {
        if self
            .mesh
            .nodes
            .iter()
            .any(|n| n.coordinates[2].abs() > 1e-15)
        {
            3
        } else {
            2
        }
    }

    fn write_dimension(&self, buf: &mut String) -> crate::Result<()> {
        let dim = self.dimension();
        writeln!(buf, "(2 {dim})")?;
        Ok(())
    }

    fn write_nodes(&self, buf: &mut String) -> crate::Result<()> {
        let n = self.mesh.nodes.len();
        if n == 0 {
            return Ok(());
        }
        // Header: (10 (zone-id first last type nd)); the data group opens on
        // the header line, as in the files Fluent writes itself, and a
        // two-dimensional mesh has two coordinates per node.
        let dim = self.dimension();
        writeln!(buf, "(10 (0 1 {n:x} 0))")?;
        writeln!(buf, "(10 (1 1 {n:x} 1 {dim})(")?;
        for node in &self.mesh.nodes {
            let [x, y, z] = node.coordinates;
            if dim == 3 {
                writeln!(buf, "{x:.10e} {y:.10e} {z:.10e}")?;
            } else {
                writeln!(buf, "{x:.10e} {y:.10e}")?;
            }
        }
        writeln!(buf, "))")?;
        Ok(())
    }

    fn write_cells(&self, buf: &mut String) -> crate::Result<()> {
        let n = self.mesh.cells.len();
        if n == 0 {
            return Ok(());
        }
        writeln!(buf, "(12 (0 1 {n:x} 0))")?;
        // Group cells by zone
        let mut zones: Vec<usize> = self.mesh.cells.iter().map(|c| c.zone_id).collect();
        zones.sort_unstable();
        zones.dedup();
        for zone_id in zones {
            let zone_cells: Vec<&FluentCell> = self
                .mesh
                .cells
                .iter()
                .filter(|c| c.zone_id == zone_id)
                .collect();
            let first = zone_cells.first().map(|c| c.id).unwrap_or(1);
            let last = zone_cells.last().map(|c| c.id).unwrap_or(1);
            let cell_type = zone_cells.first().map(|c| c.cell_type.code()).unwrap_or(1);
            writeln!(buf, "(12 ({zone_id:x} {first:x} {last:x} 1 {cell_type}))")?;
        }
        Ok(())
    }

    fn write_faces(&self, buf: &mut String) -> crate::Result<()> {
        let n = self.mesh.faces.len();
        if n == 0 {
            return Ok(());
        }
        writeln!(buf, "(13 (0 1 {n:x} 0))")?;
        writeln!(buf, "(13 (1 1 {n:x} 2 0)(")?;
        for face in &self.mesh.faces {
            // format: <face_type> <n1> <n2> ... <left_cell> <right_cell>
            write!(buf, "{:x}", face.face_type.code())?;
            for &nid in &face.node_ids {
                write!(buf, " {nid:x}")?;
            }
            writeln!(buf, " {:x} {:x}", face.left_cell, face.right_cell)?;
        }
        writeln!(buf, "))")?;
        Ok(())
    }

    fn write_zones(&self, buf: &mut String) -> crate::Result<()> {
        for (zone_id, zone_type) in &self.mesh.zones {
            writeln!(
                buf,
                "(45 ({zone_id} {} zone-{zone_id} ()))",
                zone_type.as_str()
            )?;
        }
        Ok(())
    }
}

// ── FluentReader ──────────────────────────────────────────────────────────────

/// Parses an ANSYS Fluent ASCII `.msh` file into a `FluentMesh`.
///
/// See the [module documentation](self) for the layout that is understood.
pub struct FluentReader;

impl FluentReader {
    /// Read a mesh from an ASCII Fluent `.msh` file.
    ///
    /// # Errors
    ///
    /// Fails when the file cannot be read and with [`IoError::Parse`] when its
    /// content is not an ASCII Fluent mesh [`FluentMesh`] can hold; see
    /// [`FluentReader::parse`].
    pub fn read_from_path(path: &str) -> crate::Result<FluentMesh> {
        let file = fs::File::open(path)?;
        let reader = io::BufReader::new(file);
        Self::parse(reader)
    }

    /// Parse from any `BufRead` source.
    ///
    /// # Errors
    ///
    /// Returns [`IoError::Parse`], naming the line, for binary sections, for a
    /// section that is not closed, for a number that cannot be read, for node,
    /// face or cell data that does not match the index range of its header, for
    /// faces with more than four nodes and for polyhedral cells.
    pub fn parse<R: BufRead>(mut reader: R) -> crate::Result<FluentMesh> {
        let mut bytes = Vec::new();
        io::Read::read_to_end(&mut reader, &mut bytes).map_err(IoError::Io)?;
        sections::parse_mesh(&bytes)
    }
}

// ── impl std::fmt::Write passthrough ─────────────────────────────────────────

impl From<std::fmt::Error> for IoError {
    fn from(e: std::fmt::Error) -> Self {
        IoError::General(e.to_string())
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── FluentZoneType ───────────────────────────────────────────────────────

    #[test]
    fn test_zone_type_as_str_fluid() {
        assert_eq!(FluentZoneType::Fluid.as_str(), "fluid");
    }

    #[test]
    fn test_zone_type_as_str_wall() {
        assert_eq!(FluentZoneType::Wall.as_str(), "wall");
    }

    #[test]
    fn test_zone_type_roundtrip_fluid() {
        let zt = FluentZoneType::Fluid;
        assert_eq!(FluentZoneType::from_keyword(zt.as_str()), zt);
    }

    #[test]
    fn test_zone_type_roundtrip_all() {
        let types = [
            FluentZoneType::Fluid,
            FluentZoneType::Solid,
            FluentZoneType::Wall,
            FluentZoneType::Symmetry,
            FluentZoneType::Interior,
        ];
        for zt in &types {
            assert_eq!(&FluentZoneType::from_keyword(zt.as_str()), zt);
        }
    }

    #[test]
    fn test_zone_type_unknown_defaults_to_interior() {
        assert_eq!(
            FluentZoneType::from_keyword("unknown-zone"),
            FluentZoneType::Interior
        );
    }

    // ── FluentNode ───────────────────────────────────────────────────────────

    #[test]
    fn test_fluent_node_creation() {
        let node = FluentNode::new(1, [1.0, 2.0, 3.0]);
        assert_eq!(node.id, 1);
        assert_eq!(node.coordinates, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_fluent_node_zero_coords() {
        let node = FluentNode::new(5, [0.0; 3]);
        assert_eq!(node.coordinates, [0.0; 3]);
    }

    // ── FluentFaceType ───────────────────────────────────────────────────────

    #[test]
    fn test_face_type_code_roundtrip() {
        let types = [
            FluentFaceType::Line,
            FluentFaceType::Triangle,
            FluentFaceType::Quad,
        ];
        for ft in &types {
            assert_eq!(&FluentFaceType::from_code(ft.code()), ft);
        }
    }

    #[test]
    fn test_face_type_node_count_line() {
        assert_eq!(FluentFaceType::Line.node_count(), 2);
    }

    #[test]
    fn test_face_type_node_count_tri() {
        assert_eq!(FluentFaceType::Triangle.node_count(), 3);
    }

    #[test]
    fn test_face_type_node_count_quad() {
        assert_eq!(FluentFaceType::Quad.node_count(), 4);
    }

    // ── FluentFace ───────────────────────────────────────────────────────────

    #[test]
    fn test_fluent_face_boundary_detection() {
        let face = FluentFace::new(1, FluentFaceType::Triangle, vec![1, 2, 3], 0, 5);
        assert!(face.is_boundary());
    }

    #[test]
    fn test_fluent_face_interior() {
        let face = FluentFace::new(2, FluentFaceType::Triangle, vec![1, 2, 3], 4, 5);
        assert!(!face.is_boundary());
    }

    #[test]
    fn test_fluent_face_node_ids_preserved() {
        let ids = vec![10, 20, 30];
        let face = FluentFace::new(1, FluentFaceType::Triangle, ids.clone(), 1, 2);
        assert_eq!(face.node_ids, ids);
    }

    // ── FluentCellType ───────────────────────────────────────────────────────

    #[test]
    fn test_cell_type_code_roundtrip() {
        let types = [
            FluentCellType::Tri,
            FluentCellType::Tet,
            FluentCellType::Quad,
            FluentCellType::Hex,
        ];
        for ct in &types {
            assert_eq!(&FluentCellType::from_code(ct.code()), ct);
        }
    }

    // ── FluentCell ───────────────────────────────────────────────────────────

    #[test]
    fn test_fluent_cell_creation() {
        let cell = FluentCell::new(1, FluentCellType::Tet, 2);
        assert_eq!(cell.id, 1);
        assert_eq!(cell.cell_type, FluentCellType::Tet);
        assert_eq!(cell.zone_id, 2);
    }

    // ── FluentMesh ───────────────────────────────────────────────────────────

    #[test]
    fn test_fluent_mesh_empty() {
        let mesh = FluentMesh::new();
        assert_eq!(mesh.node_count(), 0);
        assert_eq!(mesh.face_count(), 0);
        assert_eq!(mesh.cell_count(), 0);
    }

    #[test]
    fn test_fluent_mesh_add_zone() {
        let mut mesh = FluentMesh::new();
        mesh.add_zone(1, FluentZoneType::Fluid);
        assert_eq!(mesh.zones.len(), 1);
        assert_eq!(mesh.zone_type(1), Some(&FluentZoneType::Fluid));
    }

    #[test]
    fn test_fluent_mesh_zone_type_not_found() {
        let mesh = FluentMesh::new();
        assert_eq!(mesh.zone_type(99), None);
    }

    #[test]
    fn test_fluent_mesh_counts() {
        let mut mesh = FluentMesh::new();
        mesh.nodes.push(FluentNode::new(1, [0.0; 3]));
        mesh.cells.push(FluentCell::new(1, FluentCellType::Tet, 1));
        mesh.faces.push(FluentFace::new(
            1,
            FluentFaceType::Triangle,
            vec![1, 2, 3],
            0,
            1,
        ));
        assert_eq!(mesh.node_count(), 1);
        assert_eq!(mesh.cell_count(), 1);
        assert_eq!(mesh.face_count(), 1);
    }

    // ── FluentWriter ─────────────────────────────────────────────────────────

    fn make_simple_mesh() -> FluentMesh {
        let mut mesh = FluentMesh::new();
        mesh.nodes.push(FluentNode::new(1, [0.0, 0.0, 0.0]));
        mesh.nodes.push(FluentNode::new(2, [1.0, 0.0, 0.0]));
        mesh.nodes.push(FluentNode::new(3, [0.5, 1.0, 0.0]));
        mesh.cells.push(FluentCell::new(1, FluentCellType::Tri, 1));
        mesh.faces
            .push(FluentFace::new(1, FluentFaceType::Line, vec![1, 2], 0, 1));
        mesh.faces
            .push(FluentFace::new(2, FluentFaceType::Line, vec![2, 3], 0, 1));
        mesh.faces
            .push(FluentFace::new(3, FluentFaceType::Line, vec![3, 1], 0, 1));
        mesh.add_zone(1, FluentZoneType::Fluid);
        mesh
    }

    #[test]
    fn test_writer_produces_comment_section() {
        let mesh = make_simple_mesh();
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        assert!(out.contains("(0 "));
    }

    #[test]
    fn test_writer_produces_node_section() {
        let mesh = make_simple_mesh();
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        assert!(out.contains("(10 "));
    }

    #[test]
    fn test_writer_produces_cell_section() {
        let mesh = make_simple_mesh();
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        assert!(out.contains("(12 "));
    }

    #[test]
    fn test_writer_produces_face_section() {
        let mesh = make_simple_mesh();
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        assert!(out.contains("(13 "));
    }

    #[test]
    fn test_writer_produces_zone_section() {
        let mesh = make_simple_mesh();
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        assert!(out.contains("(45 "));
    }

    #[test]
    fn test_writer_node_count_hex() {
        // 3 nodes → hex "3" in node header
        let mesh = make_simple_mesh();
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        // "3" hex = 3 decimal
        assert!(out.contains("(10 (0 1 3 0))"));
    }

    // ── Write/Read roundtrip ─────────────────────────────────────────────────

    #[test]
    fn test_write_read_roundtrip_node_count() {
        let mesh = make_simple_mesh();
        let path = std::env::temp_dir().join("oxiphysics_fluent_test_roundtrip.msh");
        mesh.write(path.to_str().unwrap_or(""))
            .expect("write failed");
        let loaded = FluentMesh::read(path.to_str().unwrap_or("")).expect("read failed");
        assert_eq!(loaded.node_count(), mesh.node_count());
    }

    #[test]
    fn test_write_read_roundtrip_node_coords() {
        let mesh = make_simple_mesh();
        let path = std::env::temp_dir().join("oxiphysics_fluent_test_coords.msh");
        mesh.write(path.to_str().unwrap_or("")).unwrap();
        let loaded = FluentMesh::read(path.to_str().unwrap_or("")).unwrap();
        for (orig, loaded_node) in mesh.nodes.iter().zip(loaded.nodes.iter()) {
            for k in 0..3 {
                let diff = (orig.coordinates[k] - loaded_node.coordinates[k]).abs();
                assert!(diff < 1e-6, "coord diff={diff}");
            }
        }
    }

    #[test]
    fn test_write_read_roundtrip_cell_count() {
        let mesh = make_simple_mesh();
        let path = std::env::temp_dir().join("oxiphysics_fluent_test_cells.msh");
        mesh.write(path.to_str().unwrap_or("")).unwrap();
        let loaded = FluentMesh::read(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(loaded.cell_count(), mesh.cell_count());
    }

    #[test]
    fn test_write_read_roundtrip_zone() {
        let mesh = make_simple_mesh();
        let path = std::env::temp_dir().join("oxiphysics_fluent_test_zones.msh");
        mesh.write(path.to_str().unwrap_or("")).unwrap();
        let loaded = FluentMesh::read(path.to_str().unwrap_or("")).unwrap();
        assert!(!loaded.zones.is_empty());
        assert_eq!(loaded.zones[0].1, FluentZoneType::Fluid);
    }

    #[test]
    fn test_write_creates_file() {
        let mesh = make_simple_mesh();
        let path = std::env::temp_dir().join("oxiphysics_fluent_write_check.msh");
        mesh.write(path.to_str().unwrap_or("")).unwrap();
        assert!(path.exists());
    }

    #[test]
    fn test_read_nonexistent_file_errors() {
        let path = std::env::temp_dir().join("oxiphysics_does_not_exist.msh");
        let result = FluentMesh::read(path.to_str().unwrap_or(""));
        assert!(result.is_err());
    }

    // ── 3-D mesh ─────────────────────────────────────────────────────────────

    #[test]
    fn test_3d_mesh_dimension() {
        let mut mesh = FluentMesh::new();
        mesh.nodes.push(FluentNode::new(1, [0.0, 0.0, 1.0]));
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        assert!(out.contains("(2 3)"));
    }

    #[test]
    fn test_2d_mesh_dimension() {
        let mut mesh = FluentMesh::new();
        mesh.nodes.push(FluentNode::new(1, [1.0, 2.0, 0.0]));
        let writer = FluentWriter::new(&mesh);
        let out = writer.to_string().unwrap();
        assert!(out.contains("(2 2)"));
    }

    #[test]
    fn test_multiple_zones() {
        let mut mesh = FluentMesh::new();
        mesh.add_zone(1, FluentZoneType::Fluid);
        mesh.add_zone(2, FluentZoneType::Wall);
        assert_eq!(mesh.zones.len(), 2);
        assert_eq!(mesh.zone_type(2), Some(&FluentZoneType::Wall));
    }

    #[test]
    fn test_hex_cell_type_preserved() {
        let cell = FluentCell::new(1, FluentCellType::Hex, 1);
        assert_eq!(cell.cell_type, FluentCellType::Hex);
        assert_eq!(cell.cell_type.code(), 4);
    }

    #[test]
    fn test_quad_face_node_ids() {
        let face = FluentFace::new(1, FluentFaceType::Quad, vec![1, 2, 3, 4], 1, 2);
        assert_eq!(face.node_ids.len(), 4);
        assert!(!face.is_boundary());
    }
}
