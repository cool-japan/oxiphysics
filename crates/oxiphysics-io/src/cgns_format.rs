#![allow(clippy::manual_strip, clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CGNS (CFD General Notation System) format I/O — text-based subset.
//!
//! Provides readers and writers for a simplified text representation of CGNS
//! files, covering structured and unstructured zones, flow solutions, and
//! multi-base hierarchies.

use std::collections::HashMap;
use std::fs;
use std::io::{self, Write as IoWrite};

use crate::Error as IoError;

// ---------------------------------------------------------------------------
// CgnsNode
// ---------------------------------------------------------------------------

/// A generic CGNS tree node carrying metadata and optional data payload.
#[derive(Debug, Clone)]
pub struct CgnsNode {
    /// Unique integer identifier for this node.
    pub id: u32,
    /// Node name string (e.g. `"GridCoordinates"`).
    pub name: String,
    /// CGNS label string (e.g. `"GridCoordinates_t"`).
    pub label: String,
    /// Data-type descriptor (e.g. `"R8"`, `"I4"`).
    pub data_type: String,
    /// Floating-point data payload.
    pub data: Vec<f64>,
}

impl CgnsNode {
    /// Create a new node with the given metadata and no data.
    pub fn new(
        id: u32,
        name: impl Into<String>,
        label: impl Into<String>,
        data_type: impl Into<String>,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            label: label.into(),
            data_type: data_type.into(),
            data: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// ZoneType
// ---------------------------------------------------------------------------

/// CGNS zone topology type.
#[derive(Debug, Clone, PartialEq)]
pub enum ZoneType {
    /// Structured (curvilinear) grid topology.
    Structured,
    /// Unstructured (arbitrary connectivity) grid topology.
    Unstructured,
}

impl ZoneType {
    /// Canonical string representation used in the text format.
    pub fn as_str(&self) -> &'static str {
        match self {
            ZoneType::Structured => "Structured",
            ZoneType::Unstructured => "Unstructured",
        }
    }

    /// Parse from the string representation used in the text format.
    pub fn from_str(s: &str) -> Self {
        match s.trim() {
            "Unstructured" => ZoneType::Unstructured,
            _ => ZoneType::Structured,
        }
    }
}

// ---------------------------------------------------------------------------
// CgnsZone
// ---------------------------------------------------------------------------

/// A CGNS zone: a block of coordinates and optional element connectivity.
#[derive(Debug, Clone)]
pub struct CgnsZone {
    /// Zone name string.
    pub zone_name: String,
    /// Topology type of this zone.
    pub zone_type: ZoneType,
    /// X-coordinates of grid points.
    pub x: Vec<f64>,
    /// Y-coordinates of grid points.
    pub y: Vec<f64>,
    /// Z-coordinates of grid points.
    pub z: Vec<f64>,
    /// Flattened element connectivity array (node IDs, 1-based).
    pub elements: Vec<usize>,
}

impl CgnsZone {
    /// Create an empty zone with the given name and type.
    pub fn new(zone_name: impl Into<String>, zone_type: ZoneType) -> Self {
        Self {
            zone_name: zone_name.into(),
            zone_type,
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            elements: Vec::new(),
        }
    }

    /// Set the coordinate arrays for this zone.
    pub fn set_coords(&mut self, x: Vec<f64>, y: Vec<f64>, z: Vec<f64>) {
        self.x = x;
        self.y = y;
        self.z = z;
    }

    /// Number of grid points in this zone.
    pub fn n_points(&self) -> usize {
        self.x.len()
    }
}

// ---------------------------------------------------------------------------
// CgnsBase
// ---------------------------------------------------------------------------

/// A CGNS base node: collects zones under a named base with dimension info.
#[derive(Debug, Clone)]
pub struct CgnsBase {
    /// Base name string.
    pub base_name: String,
    /// Cell dimension (1, 2, or 3).
    pub cell_dim: u8,
    /// Physical dimension (1, 2, or 3).
    pub phys_dim: u8,
    /// Zones belonging to this base.
    pub zones: Vec<CgnsZone>,
}

impl CgnsBase {
    /// Create a new base with the given name and dimensions.
    pub fn new(base_name: impl Into<String>, cell_dim: u8, phys_dim: u8) -> Self {
        Self {
            base_name: base_name.into(),
            cell_dim,
            phys_dim,
            zones: Vec::new(),
        }
    }

    /// Append a zone to this base.
    pub fn add_zone(&mut self, zone: CgnsZone) {
        self.zones.push(zone);
    }
}

// ---------------------------------------------------------------------------
// SolutionLocation
// ---------------------------------------------------------------------------

/// Location at which a flow solution is defined.
#[derive(Debug, Clone, PartialEq)]
pub enum SolutionLocation {
    /// Data is defined at grid vertices.
    Vertex,
    /// Data is defined at element/cell centres.
    CellCenter,
}

impl SolutionLocation {
    /// Canonical string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            SolutionLocation::Vertex => "Vertex",
            SolutionLocation::CellCenter => "CellCenter",
        }
    }

    /// Parse from canonical string representation.
    pub fn from_str(s: &str) -> Self {
        match s.trim() {
            "CellCenter" => SolutionLocation::CellCenter,
            _ => SolutionLocation::Vertex,
        }
    }
}

// ---------------------------------------------------------------------------
// FlowSolution
// ---------------------------------------------------------------------------

/// A named flow solution attached to a zone.
#[derive(Debug, Clone)]
pub struct FlowSolution {
    /// Name of this solution (e.g. `"FlowSolution"`).
    pub solution_name: String,
    /// Location at which solution fields are defined.
    pub location: SolutionLocation,
    /// Named scalar fields (e.g. `"Pressure"`, `"VelocityX"`).
    pub fields: HashMap<String, Vec<f64>>,
}

impl FlowSolution {
    /// Create a new empty flow solution.
    pub fn new(solution_name: impl Into<String>, location: SolutionLocation) -> Self {
        Self {
            solution_name: solution_name.into(),
            location,
            fields: HashMap::new(),
        }
    }

    /// Insert or replace a named field.
    pub fn add_field(&mut self, name: impl Into<String>, values: Vec<f64>) {
        self.fields.insert(name.into(), values);
    }
}

// ---------------------------------------------------------------------------
// CgnsFile
// ---------------------------------------------------------------------------

/// Top-level CGNS file container with multi-base support.
#[derive(Debug, Clone)]
pub struct CgnsFile {
    /// All bases stored in this file.
    pub bases: Vec<CgnsBase>,
}

impl CgnsFile {
    /// Create an empty CGNS file.
    pub fn new() -> Self {
        Self { bases: Vec::new() }
    }

    /// Append a base to the file.
    pub fn add_base(&mut self, base: CgnsBase) {
        self.bases.push(base);
    }

    /// Write the file in the simplified text CGNS format to `path`.
    pub fn write_text(&self, path: &str) -> Result<(), IoError> {
        let mut f = fs::File::create(path).map_err(IoError::Io)?;
        writeln!(f, "CGNS_TEXT_V1").map_err(IoError::Io)?;
        writeln!(f, "N_BASES {}", self.bases.len()).map_err(IoError::Io)?;
        for base in &self.bases {
            writeln!(
                f,
                "BASE {} {} {}",
                base.base_name, base.cell_dim, base.phys_dim
            )
            .map_err(IoError::Io)?;
            writeln!(f, "N_ZONES {}", base.zones.len()).map_err(IoError::Io)?;
            for zone in &base.zones {
                writeln!(f, "ZONE {} {}", zone.zone_name, zone.zone_type.as_str())
                    .map_err(IoError::Io)?;
                writeln!(f, "N_POINTS {}", zone.x.len()).map_err(IoError::Io)?;
                // Coords
                let x_strs: Vec<String> = zone.x.iter().map(|v| v.to_string()).collect();
                writeln!(f, "X {}", x_strs.join(" ")).map_err(IoError::Io)?;
                let y_strs: Vec<String> = zone.y.iter().map(|v| v.to_string()).collect();
                writeln!(f, "Y {}", y_strs.join(" ")).map_err(IoError::Io)?;
                let z_strs: Vec<String> = zone.z.iter().map(|v| v.to_string()).collect();
                writeln!(f, "Z {}", z_strs.join(" ")).map_err(IoError::Io)?;
                // Elements
                let elem_strs: Vec<String> = zone.elements.iter().map(|v| v.to_string()).collect();
                writeln!(
                    f,
                    "ELEMENTS {} {}",
                    zone.elements.len(),
                    elem_strs.join(" ")
                )
                .map_err(IoError::Io)?;
                writeln!(f, "END_ZONE").map_err(IoError::Io)?;
            }
            writeln!(f, "END_BASE").map_err(IoError::Io)?;
        }
        writeln!(f, "END_FILE").map_err(IoError::Io)?;
        Ok(())
    }

    /// Read a CGNS file from the simplified text format at `path`.
    pub fn read_text(path: &str) -> Result<Self, IoError> {
        let text = fs::read_to_string(path).map_err(IoError::Io)?;
        CgnsReader::parse(&text)
    }
}

impl Default for CgnsFile {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CgnsWriter
// ---------------------------------------------------------------------------

/// High-level writer that assembles a structured grid with a flow solution
/// and delegates to `CgnsFile::write_text`.
#[derive(Debug)]
pub struct CgnsWriter {
    /// The file being assembled.
    pub file: CgnsFile,
}

impl CgnsWriter {
    /// Create a new writer backed by an empty `CgnsFile`.
    pub fn new() -> Self {
        Self {
            file: CgnsFile::new(),
        }
    }

    /// Write a structured grid zone (with optional flow solution) to `path`.
    ///
    /// * `base_name`  — name of the CGNS base node.
    /// * `zone_name`  — name of the zone.
    /// * `x`, `y`, `z` — coordinate arrays of equal length.
    /// * `solution`   — optional flow solution to attach.
    /// * `path`       — output file path.
    #[allow(clippy::too_many_arguments)]
    pub fn write_structured_grid(
        &mut self,
        base_name: &str,
        zone_name: &str,
        x: Vec<f64>,
        y: Vec<f64>,
        z: Vec<f64>,
        solution: Option<FlowSolution>,
        path: &str,
    ) -> Result<(), IoError> {
        let mut base = CgnsBase::new(base_name, 3, 3);
        let mut zone = CgnsZone::new(zone_name, ZoneType::Structured);
        zone.set_coords(x, y, z);
        base.add_zone(zone);
        self.file.add_base(base);

        // Write main file
        self.file.write_text(path)?;

        // If a solution is provided, write a companion solution file
        if let Some(sol) = solution {
            let sol_path = format!("{path}.sol");
            let mut sf = fs::File::create(&sol_path).map_err(IoError::Io)?;
            writeln!(
                sf,
                "SOLUTION {} {}",
                sol.solution_name,
                sol.location.as_str()
            )
            .map_err(IoError::Io)?;
            for (name, vals) in &sol.fields {
                let strs: Vec<String> = vals.iter().map(|v| v.to_string()).collect();
                writeln!(sf, "FIELD {} {}", name, strs.join(" ")).map_err(IoError::Io)?;
            }
            writeln!(sf, "END_SOLUTION").map_err(IoError::Io)?;
        }
        Ok(())
    }
}

impl Default for CgnsWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CgnsReader
// ---------------------------------------------------------------------------

/// Parser for the simplified CGNS text format produced by `CgnsFile::write_text`.
#[derive(Debug)]
pub struct CgnsReader;

impl CgnsReader {
    /// Parse CGNS text format from a string slice.
    pub fn parse(text: &str) -> Result<CgnsFile, IoError> {
        let mut file = CgnsFile::new();
        let mut lines = text.lines().peekable();

        // Expect header
        match lines.next() {
            Some(l) if l.trim() == "CGNS_TEXT_V1" => {}
            _ => {
                return Err(IoError::Io(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "missing CGNS_TEXT_V1 header",
                )));
            }
        }

        while let Some(line) = lines.next() {
            let line = line.trim();
            if line.starts_with("N_BASES") {
                // Just consume — we infer count from BASE tokens
                continue;
            }
            if line.starts_with("BASE ") {
                let parts: Vec<&str> = line.splitn(4, ' ').collect();
                let base_name = parts.get(1).copied().unwrap_or("Base").to_string();
                let cell_dim: u8 = parts.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
                let phys_dim: u8 = parts.get(3).and_then(|s| s.parse().ok()).unwrap_or(3);
                let mut base = CgnsBase::new(base_name, cell_dim, phys_dim);

                // Read zones until END_BASE
                'base_loop: loop {
                    match lines.next() {
                        None => break 'base_loop,
                        Some(l) => {
                            let l = l.trim();
                            if l == "END_BASE" {
                                break 'base_loop;
                            }
                            if l.starts_with("N_ZONES") {
                                continue;
                            }
                            if l.starts_with("ZONE ") {
                                let zparts: Vec<&str> = l.splitn(3, ' ').collect();
                                let zone_name =
                                    zparts.get(1).copied().unwrap_or("Zone").to_string();
                                let zone_type = ZoneType::from_str(
                                    zparts.get(2).copied().unwrap_or("Structured"),
                                );
                                let mut zone = CgnsZone::new(zone_name, zone_type);

                                // Read zone contents until END_ZONE
                                'zone_loop: loop {
                                    match lines.next() {
                                        None => break 'zone_loop,
                                        Some(zl) => {
                                            let zl = zl.trim();
                                            if zl == "END_ZONE" {
                                                break 'zone_loop;
                                            }
                                            if zl.starts_with("N_POINTS") {
                                                continue;
                                            }
                                            if zl.starts_with("X ") {
                                                zone.x = Self::parse_floats(&zl[2..]);
                                            } else if zl.starts_with("Y ") {
                                                zone.y = Self::parse_floats(&zl[2..]);
                                            } else if zl.starts_with("Z ") {
                                                zone.z = Self::parse_floats(&zl[2..]);
                                            } else if zl.starts_with("ELEMENTS ") {
                                                // Format: ELEMENTS <count> <v1> <v2> ...
                                                let ep: Vec<&str> = zl.splitn(3, ' ').collect();
                                                if let Some(data_str) = ep.get(2) {
                                                    zone.elements = data_str
                                                        .split_whitespace()
                                                        .filter_map(|s| s.parse().ok())
                                                        .collect();
                                                }
                                            }
                                        }
                                    }
                                }
                                base.add_zone(zone);
                            }
                        }
                    }
                }
                file.add_base(base);
            }
        }
        Ok(file)
    }

    /// Parse a space-separated list of f64 values from a string slice.
    fn parse_floats(s: &str) -> Vec<f64> {
        s.split_whitespace()
            .filter_map(|tok| tok.parse().ok())
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_path(name: &str) -> String {
        format!("/tmp/cgns_test_{name}")
    }

    // ── CgnsNode tests ───────────────────────────────────────────────────

    #[test]
    fn test_cgns_node_new() {
        let node = CgnsNode::new(1, "GridCoords", "GridCoordinates_t", "R8");
        assert_eq!(node.id, 1);
        assert_eq!(node.name, "GridCoords");
        assert_eq!(node.label, "GridCoordinates_t");
        assert_eq!(node.data_type, "R8");
        assert!(node.data.is_empty());
    }

    #[test]
    fn test_cgns_node_data_storage() {
        let mut node = CgnsNode::new(2, "Pressure", "DataArray_t", "R8");
        node.data = vec![1.0, 2.0, 3.0];
        assert_eq!(node.data.len(), 3);
        assert!((node.data[1] - 2.0).abs() < 1e-12);
    }

    // ── ZoneType tests ───────────────────────────────────────────────────

    #[test]
    fn test_zone_type_as_str() {
        assert_eq!(ZoneType::Structured.as_str(), "Structured");
        assert_eq!(ZoneType::Unstructured.as_str(), "Unstructured");
    }

    #[test]
    fn test_zone_type_from_str_structured() {
        assert_eq!(ZoneType::from_str("Structured"), ZoneType::Structured);
        assert_eq!(ZoneType::from_str("anything"), ZoneType::Structured);
    }

    #[test]
    fn test_zone_type_from_str_unstructured() {
        assert_eq!(ZoneType::from_str("Unstructured"), ZoneType::Unstructured);
    }

    // ── CgnsZone tests ───────────────────────────────────────────────────

    #[test]
    fn test_zone_creation() {
        let zone = CgnsZone::new("Zone1", ZoneType::Structured);
        assert_eq!(zone.zone_name, "Zone1");
        assert_eq!(zone.zone_type, ZoneType::Structured);
        assert!(zone.x.is_empty());
    }

    #[test]
    fn test_zone_set_coords() {
        let mut zone = CgnsZone::new("Z", ZoneType::Structured);
        zone.set_coords(vec![0.0, 1.0], vec![0.0, 0.0], vec![0.0, 0.0]);
        assert_eq!(zone.n_points(), 2);
        assert!((zone.x[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_zone_n_points() {
        let mut zone = CgnsZone::new("Z", ZoneType::Unstructured);
        zone.x = vec![1.0, 2.0, 3.0];
        assert_eq!(zone.n_points(), 3);
    }

    #[test]
    fn test_zone_elements() {
        let mut zone = CgnsZone::new("Z", ZoneType::Unstructured);
        zone.elements = vec![1, 2, 3, 4];
        assert_eq!(zone.elements.len(), 4);
    }

    // ── CgnsBase tests ───────────────────────────────────────────────────

    #[test]
    fn test_base_creation() {
        let base = CgnsBase::new("Base1", 3, 3);
        assert_eq!(base.base_name, "Base1");
        assert_eq!(base.cell_dim, 3);
        assert_eq!(base.phys_dim, 3);
        assert!(base.zones.is_empty());
    }

    #[test]
    fn test_base_add_zone() {
        let mut base = CgnsBase::new("Base1", 3, 3);
        base.add_zone(CgnsZone::new("Z1", ZoneType::Structured));
        base.add_zone(CgnsZone::new("Z2", ZoneType::Unstructured));
        assert_eq!(base.zones.len(), 2);
    }

    // ── SolutionLocation tests ───────────────────────────────────────────

    #[test]
    fn test_solution_location_as_str() {
        assert_eq!(SolutionLocation::Vertex.as_str(), "Vertex");
        assert_eq!(SolutionLocation::CellCenter.as_str(), "CellCenter");
    }

    #[test]
    fn test_solution_location_from_str() {
        assert_eq!(
            SolutionLocation::from_str("Vertex"),
            SolutionLocation::Vertex
        );
        assert_eq!(
            SolutionLocation::from_str("CellCenter"),
            SolutionLocation::CellCenter
        );
        assert_eq!(
            SolutionLocation::from_str("other"),
            SolutionLocation::Vertex
        );
    }

    // ── FlowSolution tests ───────────────────────────────────────────────

    #[test]
    fn test_flow_solution_new() {
        let sol = FlowSolution::new("FlowSolution", SolutionLocation::Vertex);
        assert_eq!(sol.solution_name, "FlowSolution");
        assert_eq!(sol.location, SolutionLocation::Vertex);
        assert!(sol.fields.is_empty());
    }

    #[test]
    fn test_flow_solution_add_field() {
        let mut sol = FlowSolution::new("Sol", SolutionLocation::CellCenter);
        sol.add_field("Pressure", vec![1.0, 2.0, 3.0]);
        assert!(sol.fields.contains_key("Pressure"));
        assert_eq!(sol.fields["Pressure"].len(), 3);
    }

    #[test]
    fn test_flow_solution_multiple_fields() {
        let mut sol = FlowSolution::new("Sol", SolutionLocation::Vertex);
        sol.add_field("P", vec![1.0]);
        sol.add_field("T", vec![300.0]);
        sol.add_field("U", vec![0.5]);
        assert_eq!(sol.fields.len(), 3);
    }

    // ── CgnsFile write/read roundtrip tests ─────────────────────────────

    #[test]
    fn test_file_write_read_empty() {
        let file = CgnsFile::new();
        let path = tmp_path("empty");
        file.write_text(&path).unwrap();
        let restored = CgnsFile::read_text(&path).unwrap();
        assert_eq!(restored.bases.len(), 0);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_file_write_read_single_base() {
        let mut file = CgnsFile::new();
        let base = CgnsBase::new("MainBase", 3, 3);
        file.add_base(base);
        let path = tmp_path("single_base");
        file.write_text(&path).unwrap();
        let restored = CgnsFile::read_text(&path).unwrap();
        assert_eq!(restored.bases.len(), 1);
        assert_eq!(restored.bases[0].base_name, "MainBase");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_file_write_read_coords_roundtrip() {
        let mut file = CgnsFile::new();
        let mut base = CgnsBase::new("B", 3, 3);
        let mut zone = CgnsZone::new("Z1", ZoneType::Structured);
        zone.set_coords(vec![0.0, 1.0, 2.0], vec![0.0, 0.5, 1.0], vec![0.0; 3]);
        base.add_zone(zone);
        file.add_base(base);
        let path = tmp_path("coords_rt");
        file.write_text(&path).unwrap();
        let restored = CgnsFile::read_text(&path).unwrap();
        let rz = &restored.bases[0].zones[0];
        assert_eq!(rz.x.len(), 3);
        assert!((rz.x[2] - 2.0).abs() < 1e-10);
        assert!((rz.y[1] - 0.5).abs() < 1e-10);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_file_write_read_zone_type() {
        let mut file = CgnsFile::new();
        let mut base = CgnsBase::new("B", 3, 3);
        let mut zone = CgnsZone::new("Unstr", ZoneType::Unstructured);
        zone.elements = vec![1, 2, 3];
        base.add_zone(zone);
        file.add_base(base);
        let path = tmp_path("zone_type");
        file.write_text(&path).unwrap();
        let restored = CgnsFile::read_text(&path).unwrap();
        assert_eq!(restored.bases[0].zones[0].zone_type, ZoneType::Unstructured);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_file_write_read_elements() {
        let mut file = CgnsFile::new();
        let mut base = CgnsBase::new("B", 3, 3);
        let mut zone = CgnsZone::new("Z", ZoneType::Unstructured);
        zone.elements = vec![1, 2, 3, 4, 5, 6];
        base.add_zone(zone);
        file.add_base(base);
        let path = tmp_path("elements");
        file.write_text(&path).unwrap();
        let restored = CgnsFile::read_text(&path).unwrap();
        assert_eq!(restored.bases[0].zones[0].elements, vec![1, 2, 3, 4, 5, 6]);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_file_write_read_multi_zone() {
        let mut file = CgnsFile::new();
        let mut base = CgnsBase::new("B", 3, 3);
        base.add_zone(CgnsZone::new("Zone1", ZoneType::Structured));
        base.add_zone(CgnsZone::new("Zone2", ZoneType::Unstructured));
        file.add_base(base);
        let path = tmp_path("multi_zone");
        file.write_text(&path).unwrap();
        let restored = CgnsFile::read_text(&path).unwrap();
        assert_eq!(restored.bases[0].zones.len(), 2);
        assert_eq!(restored.bases[0].zones[1].zone_name, "Zone2");
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_file_write_read_multi_base() {
        let mut file = CgnsFile::new();
        file.add_base(CgnsBase::new("Base1", 3, 3));
        file.add_base(CgnsBase::new("Base2", 2, 2));
        let path = tmp_path("multi_base");
        file.write_text(&path).unwrap();
        let restored = CgnsFile::read_text(&path).unwrap();
        assert_eq!(restored.bases.len(), 2);
        assert_eq!(restored.bases[1].base_name, "Base2");
        assert_eq!(restored.bases[1].cell_dim, 2);
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_file_invalid_header() {
        let result = CgnsReader::parse("NOT_CGNS\n");
        assert!(result.is_err());
    }

    #[test]
    fn test_file_default() {
        let f = CgnsFile::default();
        assert!(f.bases.is_empty());
    }

    // ── CgnsWriter tests ─────────────────────────────────────────────────

    #[test]
    fn test_writer_structured_grid_no_solution() {
        let mut writer = CgnsWriter::new();
        let path = tmp_path("writer_noSol");
        writer
            .write_structured_grid(
                "B",
                "Z",
                vec![0.0, 1.0],
                vec![0.0, 0.0],
                vec![0.0, 0.0],
                None,
                &path,
            )
            .unwrap();
        let content = fs::read_to_string(&path).unwrap();
        assert!(content.contains("CGNS_TEXT_V1"));
        assert!(content.contains("Structured"));
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn test_writer_structured_grid_with_solution() {
        let mut writer = CgnsWriter::new();
        let mut sol = FlowSolution::new("Sol", SolutionLocation::Vertex);
        sol.add_field("Pressure", vec![101325.0, 101300.0]);
        let path = tmp_path("writer_sol");
        writer
            .write_structured_grid(
                "B",
                "Z",
                vec![0.0, 1.0],
                vec![0.0, 0.0],
                vec![0.0, 0.0],
                Some(sol),
                &path,
            )
            .unwrap();
        let sol_path = format!("{path}.sol");
        let sol_content = fs::read_to_string(&sol_path).unwrap();
        assert!(sol_content.contains("Pressure"));
        let _ = fs::remove_file(&path);
        let _ = fs::remove_file(&sol_path);
    }

    #[test]
    fn test_writer_default() {
        let w = CgnsWriter::default();
        assert!(w.file.bases.is_empty());
    }

    // ── CgnsReader tests ─────────────────────────────────────────────────

    #[test]
    fn test_reader_parse_minimal() {
        let text = "CGNS_TEXT_V1\nN_BASES 0\nEND_FILE\n";
        let file = CgnsReader::parse(text).unwrap();
        assert_eq!(file.bases.len(), 0);
    }

    #[test]
    fn test_reader_parse_with_zone_coords() {
        let text = concat!(
            "CGNS_TEXT_V1\n",
            "N_BASES 1\n",
            "BASE MyBase 3 3\n",
            "N_ZONES 1\n",
            "ZONE Z1 Structured\n",
            "N_POINTS 2\n",
            "X 0 1\n",
            "Y 0 0\n",
            "Z 0 0\n",
            "ELEMENTS 0 \n",
            "END_ZONE\n",
            "END_BASE\n",
            "END_FILE\n",
        );
        let file = CgnsReader::parse(text).unwrap();
        assert_eq!(file.bases[0].zones[0].x.len(), 2);
    }
}
