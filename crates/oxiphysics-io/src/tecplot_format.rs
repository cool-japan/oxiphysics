// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tecplot ASCII format I/O.
//!
//! Supports reading and writing Tecplot `.dat` / `.plt` ASCII files with
//! ORDERED and FE (finite element) zone types.  The implementation covers the
//! minimal subset of the Tecplot ASCII specification needed for typical physics
//! simulation output.
//!
//! Finite-element zones (`FETRIANGLE`, `FEQUADRILATERAL`, `FETETRAHEDRON`,
//! `FEBRICK`) carry their element connectivity: the writer emits one line of
//! 1-based node numbers per element after the nodal data and the reader
//! returns it, 0-based, in [`TecplotZone::connectivity`].
//!
//! # Quick start
//!
//! ```
//! use oxiphysics_io::tecplot_format::{
//!     TecplotDataset, TecplotVariable, TecplotZone, TecplotZoneType,
//! };
//!
//! let mut ds = TecplotDataset::new("My Dataset");
//! ds.variables = vec!["X".to_string(), "Y".to_string(), "P".to_string()];
//! let mut zone = TecplotZone::new_ordered("Zone 1", 3, 1, 1);
//! zone.variables.push(TecplotVariable::new("X", vec![0.0, 1.0, 2.0]));
//! zone.variables.push(TecplotVariable::new("Y", vec![0.0, 0.0, 0.0]));
//! zone.variables.push(TecplotVariable::new("P", vec![1.0, 2.0, 3.0]));
//! ds.zones.push(zone);
//!
//! // A finite-element zone: two triangles sharing an edge.
//! let mut fe = TecplotZone::new_fe(
//!     "Mesh",
//!     TecplotZoneType::FETriangle,
//!     4,
//!     2,
//!     vec![0, 1, 2, 1, 3, 2],
//! );
//! fe.variables.push(TecplotVariable::new("X", vec![0.0, 1.0, 0.0, 1.0]));
//! fe.variables.push(TecplotVariable::new("Y", vec![0.0, 0.0, 1.0, 1.0]));
//! fe.variables.push(TecplotVariable::new("P", vec![1.0, 2.0, 3.0, 4.0]));
//! ds.zones.push(fe);
//!
//! let path = std::env::temp_dir().join("oxiphysics_tecplot_quick_start.dat");
//! let path = path.to_str().unwrap();
//! ds.write(path).unwrap();
//! let back = TecplotDataset::read(path).unwrap();
//! assert_eq!(back.zones[1].connectivity, vec![0, 1, 2, 1, 3, 2]);
//! std::fs::remove_file(path).ok();
//! ```

use std::fmt::Write as FmtWrite;
use std::fs;
use std::io::{self, BufRead};

use crate::Error as IoError;

// ---------------------------------------------------------------------------
// TecplotZoneType
// ---------------------------------------------------------------------------

/// Tecplot zone topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TecplotZoneType {
    /// Structured (IJK-ordered) zone.
    Ordered,
    /// Unstructured zone — triangular elements.
    FETriangle,
    /// Unstructured zone — quadrilateral elements.
    FEQuad,
    /// Unstructured zone — tetrahedral elements.
    FETetra,
    /// Unstructured zone — hexahedral (brick) elements.
    FEBrick,
}

impl TecplotZoneType {
    /// Return the Tecplot ASCII keyword for this zone type.
    pub fn as_tecplot_str(&self) -> &'static str {
        match self {
            TecplotZoneType::Ordered => "ORDERED",
            TecplotZoneType::FETriangle => "FETRIANGLE",
            TecplotZoneType::FEQuad => "FEQUADRILATERAL",
            TecplotZoneType::FETetra => "FETETRAHEDRON",
            TecplotZoneType::FEBrick => "FEBRICK",
        }
    }

    /// Parse from a Tecplot ASCII keyword (case-insensitive).
    pub fn from_keyword(s: &str) -> Self {
        Self::from(s)
    }

    /// Number of nodes that make up one element of this zone type: 3 for
    /// `FETriangle`, 4 for `FEQuad` and `FETetra`, 8 for `FEBrick`, and 0 for
    /// `Ordered` zones, which have no element connectivity.
    pub fn nodes_per_element(&self) -> usize {
        match self {
            TecplotZoneType::Ordered => 0,
            TecplotZoneType::FETriangle => 3,
            TecplotZoneType::FEQuad => 4,
            TecplotZoneType::FETetra => 4,
            TecplotZoneType::FEBrick => 8,
        }
    }

    /// Parse the element type keyword (`ET=`) of a legacy
    /// `F=FEPOINT` / `F=FEBLOCK` zone header (case-insensitive).
    ///
    /// Returns `None` for keywords that name no supported element type.
    pub fn from_element_type(s: &str) -> Option<Self> {
        match s.trim().to_uppercase().as_str() {
            "TRIANGLE" => Some(TecplotZoneType::FETriangle),
            "QUADRILATERAL" | "QUAD" => Some(TecplotZoneType::FEQuad),
            "TETRAHEDRON" | "TETRA" => Some(TecplotZoneType::FETetra),
            "BRICK" => Some(TecplotZoneType::FEBrick),
            _ => None,
        }
    }
}

impl From<&str> for TecplotZoneType {
    fn from(s: &str) -> Self {
        match s.trim().to_uppercase().as_str() {
            "ORDERED" => TecplotZoneType::Ordered,
            "FETRIANGLE" | "FE_TRIANGLE" => TecplotZoneType::FETriangle,
            "FEQUADRILATERAL" | "FEQUAD" | "FE_QUAD" => TecplotZoneType::FEQuad,
            "FETETRAHEDRON" | "FETETRA" | "FE_TETRA" => TecplotZoneType::FETetra,
            "FEBRICK" | "FE_BRICK" => TecplotZoneType::FEBrick,
            _ => TecplotZoneType::Ordered,
        }
    }
}

// ---------------------------------------------------------------------------
// TecplotVariable
// ---------------------------------------------------------------------------

/// A single named variable (field) in a Tecplot zone.
#[derive(Debug, Clone, PartialEq)]
pub struct TecplotVariable {
    /// Variable name (e.g. `"X"`, `"Pressure"`, `"Velocity_X"`).
    pub name: String,
    /// Per-node/per-cell data values.
    pub data: Vec<f64>,
}

impl TecplotVariable {
    /// Create a new variable with the given name and data.
    pub fn new(name: impl Into<String>, data: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            data,
        }
    }

    /// Return the number of data values stored.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Return `true` when no data values are stored.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

// ---------------------------------------------------------------------------
// TecplotZone
// ---------------------------------------------------------------------------

/// A single Tecplot zone containing one or more variables.
#[derive(Debug, Clone)]
pub struct TecplotZone {
    /// Zone title string.
    pub title: String,
    /// Zone topology.
    pub zone_type: TecplotZoneType,
    /// Number of nodes in the I direction (ordered zones).
    pub i_dim: usize,
    /// Number of nodes in the J direction (ordered zones).
    pub j_dim: usize,
    /// Number of nodes in the K direction (ordered zones).
    pub k_dim: usize,
    /// Variable data for this zone.
    pub variables: Vec<TecplotVariable>,
    /// Number of nodes (FE zones).
    pub n_nodes: usize,
    /// Number of elements (FE zones).
    pub n_elements: usize,
    /// Element connectivity of an FE zone: the 0-based node indices of every
    /// element, concatenated.  Element `e` occupies
    /// `connectivity[e * npe..(e + 1) * npe]` with
    /// `npe = zone_type.nodes_per_element()`, so the length must equal
    /// `n_elements * npe`.  Empty (and ignored) for ordered zones.
    pub connectivity: Vec<usize>,
}

impl TecplotZone {
    /// Create a new ordered zone with the given IJK dimensions.
    pub fn new_ordered(title: impl Into<String>, i_dim: usize, j_dim: usize, k_dim: usize) -> Self {
        Self {
            title: title.into(),
            zone_type: TecplotZoneType::Ordered,
            i_dim,
            j_dim,
            k_dim,
            variables: Vec::new(),
            n_nodes: i_dim.saturating_mul(j_dim).saturating_mul(k_dim),
            n_elements: 0,
            connectivity: Vec::new(),
        }
    }

    /// Create a new finite-element zone.
    ///
    /// `connectivity` lists the 0-based node indices of every element,
    /// concatenated (see [`TecplotZone::connectivity`]); it is checked against
    /// `n_elements` and `n_nodes` when the zone is written.
    pub fn new_fe(
        title: impl Into<String>,
        zone_type: TecplotZoneType,
        n_nodes: usize,
        n_elements: usize,
        connectivity: Vec<usize>,
    ) -> Self {
        Self {
            title: title.into(),
            zone_type,
            i_dim: 0,
            j_dim: 0,
            k_dim: 0,
            variables: Vec::new(),
            n_nodes,
            n_elements,
            connectivity,
        }
    }

    /// The 0-based node indices of element `element`, or `None` when the zone
    /// has no such element (ordered zone, index out of range, or a connectivity
    /// array that is shorter than `n_elements` elements).
    pub fn element_nodes(&self, element: usize) -> Option<&[usize]> {
        let npe = self.zone_type.nodes_per_element();
        if npe == 0 || element >= self.n_elements {
            return None;
        }
        let start = element.checked_mul(npe)?;
        self.connectivity.get(start..start.checked_add(npe)?)
    }

    /// Total number of data points expected in this zone.
    ///
    /// The product `I * J * K` of an ordered zone saturates at `usize::MAX`, so
    /// the absurd dimensions of a corrupt file header cannot overflow.
    pub fn point_count(&self) -> usize {
        match self.zone_type {
            TecplotZoneType::Ordered => self
                .i_dim
                .saturating_mul(self.j_dim)
                .saturating_mul(self.k_dim),
            _ => self.n_nodes,
        }
    }
}

// ---------------------------------------------------------------------------
// TecplotDataset
// ---------------------------------------------------------------------------

/// A complete Tecplot dataset (title + variable list + zones).
#[derive(Debug, Clone, Default)]
pub struct TecplotDataset {
    /// Dataset title.
    pub title: String,
    /// Ordered list of variable names.
    pub variables: Vec<String>,
    /// Zones contained in the dataset.
    pub zones: Vec<TecplotZone>,
}

impl TecplotDataset {
    /// Create an empty dataset with a title.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            variables: Vec::new(),
            zones: Vec::new(),
        }
    }

    /// Write the dataset to `path` in Tecplot ASCII format.
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn write(&self, path: &str) -> Result<(), IoError> {
        let writer = TecplotWriter::new();
        writer.write(self, path)
    }

    /// Read a Tecplot ASCII file from `path` and return the parsed dataset.
    ///
    /// Returns an error if the file cannot be opened or parsed.
    pub fn read(path: &str) -> Result<Self, IoError> {
        TecplotReader::new().parse(path)
    }
}

// ---------------------------------------------------------------------------
// TecplotWriter
// ---------------------------------------------------------------------------

/// Writes a `TecplotDataset` to a Tecplot ASCII file.
#[derive(Debug, Clone, Default)]
pub struct TecplotWriter;

impl TecplotWriter {
    /// Create a new writer.
    pub fn new() -> Self {
        Self
    }

    /// Serialize `dataset` to `path` using point data packing (one record per
    /// node, all variables interleaved).
    ///
    /// Finite-element zones are followed by their element connectivity, one
    /// line of 1-based node numbers per element, as the Tecplot ASCII format
    /// requires.
    ///
    /// The whole file is assembled in memory first, so nothing is written when
    /// a zone is invalid; the error is [`IoError::General`].  A zone is
    /// invalid when
    ///
    /// * one of its variables does not hold exactly one value per point
    ///   (`point_count()`): the writer neither pads missing values nor drops
    ///   surplus ones, and a zone whose declared size is far larger than its
    ///   data cannot make it produce an unbounded amount of output, or
    /// * (FE zones) its `connectivity` does not hold exactly
    ///   `n_elements * nodes_per_element` entries or refers to a node index
    ///   `>= n_nodes`.
    pub fn write(&self, dataset: &TecplotDataset, path: &str) -> Result<(), IoError> {
        let mut buf = String::new();

        // --- File header ---
        let _ = writeln!(buf, "TITLE = \"{}\"", dataset.title);

        if !dataset.variables.is_empty() {
            let vars: Vec<String> = dataset
                .variables
                .iter()
                .map(|v| format!("\"{v}\""))
                .collect();
            let _ = writeln!(buf, "VARIABLES = {}", vars.join(", "));
        }

        // --- Zones ---
        for zone in &dataset.zones {
            self.write_zone(&mut buf, zone)?;
        }

        fs::write(path, buf).map_err(IoError::Io)
    }

    /// Append one zone to `buf`.
    fn write_zone(&self, buf: &mut String, zone: &TecplotZone) -> Result<(), IoError> {
        match zone.zone_type {
            TecplotZoneType::Ordered => {
                Self::validate_variable_lengths(zone)?;
                writeln!(
                    buf,
                    "ZONE T=\"{}\", I={}, J={}, K={}, DATAPACKING=POINT",
                    zone.title, zone.i_dim, zone.j_dim, zone.k_dim
                )
                .map_err(Self::format_error)?;
                self.write_point_data(buf, zone)?;
            }
            _ => {
                Self::validate_variable_lengths(zone)?;
                Self::validate_connectivity(zone)?;
                writeln!(
                    buf,
                    "ZONE T=\"{}\", N={}, E={}, ZONETYPE={}, DATAPACKING=POINT",
                    zone.title,
                    zone.n_nodes,
                    zone.n_elements,
                    zone.zone_type.as_tecplot_str()
                )
                .map_err(Self::format_error)?;
                self.write_point_data(buf, zone)?;
                self.write_connectivity(buf, zone)?;
            }
        }
        Ok(())
    }

    /// Convert a formatting failure into the crate error type.
    fn format_error(err: std::fmt::Error) -> IoError {
        IoError::General(format!("failed to format Tecplot output: {err}"))
    }

    /// Check that every variable of the zone holds exactly one value per point.
    fn validate_variable_lengths(zone: &TecplotZone) -> Result<(), IoError> {
        match zone
            .variables
            .iter()
            .find(|variable| variable.data.len() != zone.point_count())
        {
            Some(variable) => Err(Self::length_mismatch(zone, variable)),
            None => Ok(()),
        }
    }

    /// The error for a variable that does not hold one value per point.
    fn length_mismatch(zone: &TecplotZone, variable: &TecplotVariable) -> IoError {
        IoError::General(format!(
            "Tecplot zone \"{}\": variable \"{}\" holds {} values but the zone has {} points",
            zone.title,
            variable.name,
            variable.data.len(),
            zone.point_count()
        ))
    }

    /// Check that the connectivity of an FE zone matches its element count and
    /// only refers to existing nodes.
    fn validate_connectivity(zone: &TecplotZone) -> Result<(), IoError> {
        let npe = zone.zone_type.nodes_per_element();
        let expected = zone.n_elements.checked_mul(npe).ok_or_else(|| {
            IoError::General(format!(
                "Tecplot zone \"{}\": element count {} is too large",
                zone.title, zone.n_elements
            ))
        })?;
        if zone.connectivity.len() != expected {
            return Err(IoError::General(format!(
                "Tecplot zone \"{}\": {} {} elements need {} connectivity entries \
                 ({} nodes each) but {} were supplied",
                zone.title,
                zone.n_elements,
                zone.zone_type.as_tecplot_str(),
                expected,
                npe,
                zone.connectivity.len()
            )));
        }
        if let Some(&node) = zone.connectivity.iter().find(|&&n| n >= zone.n_nodes) {
            return Err(IoError::General(format!(
                "Tecplot zone \"{}\": connectivity refers to node {} but the zone has only {} nodes",
                zone.title, node, zone.n_nodes
            )));
        }
        Ok(())
    }

    /// Write variable data in POINT packing (all vars for node 0, then node 1, …).
    ///
    /// The variables must already have been checked with
    /// [`Self::validate_variable_lengths`]; every value is formatted straight
    /// into `buf` without a temporary string per value or per row.
    fn write_point_data(&self, buf: &mut String, zone: &TecplotZone) -> Result<(), IoError> {
        if zone.variables.is_empty() {
            return Ok(());
        }
        for point in 0..zone.point_count() {
            for (column, variable) in zone.variables.iter().enumerate() {
                if column > 0 {
                    buf.push(' ');
                }
                let Some(&value) = variable.data.get(point) else {
                    return Err(Self::length_mismatch(zone, variable));
                };
                write!(buf, "{value:.15e}").map_err(Self::format_error)?;
            }
            buf.push('\n');
        }
        Ok(())
    }

    /// Write the element connectivity of an FE zone: one line per element with
    /// 1-based node numbers.
    fn write_connectivity(&self, buf: &mut String, zone: &TecplotZone) -> Result<(), IoError> {
        let npe = zone.zone_type.nodes_per_element();
        if npe == 0 {
            return Ok(());
        }
        for element in zone.connectivity.chunks_exact(npe) {
            for (column, &node) in element.iter().enumerate() {
                if column > 0 {
                    buf.push(' ');
                }
                write!(buf, "{}", node.saturating_add(1)).map_err(Self::format_error)?;
            }
            buf.push('\n');
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// TecplotReader
// ---------------------------------------------------------------------------

/// Parses a subset of Tecplot ASCII files.
///
/// Handles `TITLE`, `VARIABLES`, and `ZONE` header lines followed by
/// whitespace-delimited data blocks: for ordered zones the floating-point
/// point data, for finite-element zones the nodal data followed by the
/// element connectivity (1-based node numbers in the file, returned 0-based).
///
/// FE zones may be declared with `ZONETYPE=FETriangle` (and the other
/// `FE…` keywords) or with the legacy `F=FEPOINT, ET=TRIANGLE` form, and the
/// node/element counts may be spelled `N=`/`E=` or `Nodes=`/`Elements=`.
/// Both `DATAPACKING=POINT` (the variables of one node side by side) and
/// `DATAPACKING=BLOCK` (all values of one variable after another; also the
/// legacy `F=BLOCK` / `F=FEBLOCK`) are understood.  Cell-centred variables
/// (`VARLOCATION=(…=CELLCENTERED)`) are rejected with [`IoError::Parse`]
/// rather than being read with the wrong layout.
#[derive(Debug, Clone, Default)]
pub struct TecplotReader;

/// How the variables of a zone's data section are interleaved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DataPacking {
    /// All variables of node 0, then all variables of node 1, ….
    Point,
    /// All nodal values of variable 0, then all of variable 1, ….
    Block,
}

/// The data records of a zone whose header has been read.
///
/// Records are consumed token by token, so the line structure of the file is
/// irrelevant: the first `point_count * variables` numbers are nodal values,
/// the next `n_elements * nodes_per_element` integers are connectivity.
struct ZoneAccumulator {
    zone: TecplotZone,
    packing: DataPacking,
    /// Number of data values expected before the connectivity block.
    expected_values: usize,
    /// Number of connectivity entries expected after the data values.
    expected_connectivity: usize,
    values: Vec<f64>,
    /// Connectivity as read from the file (1-based).
    connectivity: Vec<usize>,
}

impl ZoneAccumulator {
    fn new(zone: TecplotZone, packing: DataPacking, n_variables: usize) -> Self {
        let expected_values = zone.point_count().saturating_mul(n_variables);
        let expected_connectivity = zone
            .n_elements
            .saturating_mul(zone.zone_type.nodes_per_element());
        Self {
            zone,
            packing,
            expected_values,
            expected_connectivity,
            values: Vec::new(),
            connectivity: Vec::new(),
        }
    }

    /// Consume one data line.
    ///
    /// A line that does not start with a number (`AUXDATA`, `TEXT`, … records)
    /// carries no zone data and is skipped.  Tokens beyond what the zone
    /// expects are ignored.
    fn push_line(&mut self, line: &str) -> Result<(), IoError> {
        let mut tokens = line.split_whitespace().peekable();
        match tokens.peek() {
            Some(first) if first.parse::<f64>().is_ok() => {}
            _ => return Ok(()),
        }
        for token in tokens {
            if self.values.len() < self.expected_values {
                if let Ok(value) = token.parse::<f64>() {
                    self.values.push(value);
                }
            } else if self.connectivity.len() < self.expected_connectivity {
                let node: usize = token.parse().map_err(|_| {
                    IoError::Parse(format!(
                        "Tecplot zone \"{}\": invalid connectivity entry \"{}\"",
                        self.zone.title, token
                    ))
                })?;
                self.connectivity.push(node);
            } else {
                break;
            }
        }
        Ok(())
    }

    /// Distribute the collected values into the zone's variables according to
    /// the data packing and convert/validate the connectivity.
    fn finish(self) -> Result<TecplotZone, IoError> {
        let ZoneAccumulator {
            mut zone,
            packing,
            expected_connectivity,
            values,
            connectivity,
            ..
        } = self;

        let n_variables = zone.variables.len();
        let points_per_variable = zone.point_count();
        if n_variables > 0 {
            for (idx, &value) in values.iter().enumerate() {
                let variable = match packing {
                    DataPacking::Point => idx % n_variables,
                    DataPacking::Block => idx / points_per_variable.max(1),
                };
                if let Some(target) = zone.variables.get_mut(variable) {
                    target.data.push(value);
                }
            }
        }

        if connectivity.len() < expected_connectivity {
            return Err(IoError::Parse(format!(
                "Tecplot zone \"{}\": expected {} connectivity entries ({} elements x {} nodes) \
                 but the file provides {}",
                zone.title,
                expected_connectivity,
                zone.n_elements,
                zone.zone_type.nodes_per_element(),
                connectivity.len()
            )));
        }
        zone.connectivity.reserve(connectivity.len());
        for node in connectivity {
            if node == 0 || node > zone.n_nodes {
                return Err(IoError::Parse(format!(
                    "Tecplot zone \"{}\": connectivity node {} is outside 1..={}",
                    zone.title, node, zone.n_nodes
                )));
            }
            zone.connectivity.push(node - 1);
        }
        Ok(zone)
    }
}

impl TecplotReader {
    /// Create a new reader.
    pub fn new() -> Self {
        Self
    }

    /// Parse the Tecplot ASCII file at `path` and return a `TecplotDataset`.
    ///
    /// Returns [`IoError::Parse`] when an FE zone's connectivity is truncated,
    /// is not made of integers, or refers to a node outside `1..=N`.
    pub fn parse(&self, path: &str) -> Result<TecplotDataset, IoError> {
        let file = fs::File::open(path).map_err(IoError::Io)?;
        let reader = io::BufReader::new(file);

        let mut dataset = TecplotDataset::new("");
        let mut current: Option<ZoneAccumulator> = None;

        for line_res in reader.lines() {
            let line = line_res.map_err(IoError::Io)?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // ── TITLE ──────────────────────────────────────────────────────
            if Self::starts_with_keyword(trimmed, "TITLE") {
                dataset.title =
                    Self::extract_quoted_value(trimmed).unwrap_or_else(|| trimmed.to_string());
                continue;
            }

            // ── VARIABLES ──────────────────────────────────────────────────
            if Self::starts_with_keyword(trimmed, "VARIABLES") {
                dataset.variables = Self::parse_variables_line(trimmed);
                continue;
            }

            // ── ZONE ───────────────────────────────────────────────────────
            if Self::starts_with_keyword(trimmed, "ZONE") {
                // Flush previous zone
                if let Some(previous) = current.take() {
                    dataset.zones.push(previous.finish()?);
                }
                let (zone, packing) = Self::parse_zone_header(trimmed, &dataset.variables)?;
                current = Some(ZoneAccumulator::new(zone, packing, dataset.variables.len()));
                continue;
            }

            // ── Data line ──────────────────────────────────────────────────
            if let Some(accumulator) = current.as_mut() {
                accumulator.push_line(trimmed)?;
            }
        }

        // Flush last zone
        if let Some(last) = current.take() {
            dataset.zones.push(last.finish()?);
        }

        Ok(dataset)
    }

    /// Parse a ZONE header line and return a skeleton `TecplotZone` together
    /// with the data packing of its data section.
    fn parse_zone_header(
        line: &str,
        var_names: &[String],
    ) -> Result<(TecplotZone, DataPacking), IoError> {
        let params = Self::parse_header_pairs(line);
        let get = |keys: &[&str]| -> Option<&str> {
            params
                .iter()
                .find(|(key, _)| keys.iter().any(|wanted| key.eq_ignore_ascii_case(wanted)))
                .map(|(_, value)| value.as_str())
        };
        let get_count = |keys: &[&str], default: usize| -> usize {
            get(keys).and_then(|v| v.parse().ok()).unwrap_or(default)
        };

        let title = get(&["T", "TITLE"]).unwrap_or_default().to_string();

        let i_dim = get_count(&["I"], 1);
        let j_dim = get_count(&["J"], 1);
        let k_dim = get_count(&["K"], 1);

        let n_nodes = get_count(&["N", "NODES"], 0);
        let n_elements = get_count(&["E", "ELEMENTS"], 0);

        let zone_type = if let Some(keyword) = get(&["ZONETYPE"]) {
            TecplotZoneType::from_keyword(keyword)
        } else {
            match get(&["F"]).map(str::to_uppercase).as_deref() {
                // Legacy finite-element header: F=FEPOINT/FEBLOCK, ET=<element>.
                // Tecplot's default element type is TRIANGLE.
                Some("FEPOINT") | Some("FEBLOCK") => match get(&["ET"]) {
                    Some(element_type) => TecplotZoneType::from_element_type(element_type)
                        .unwrap_or(TecplotZoneType::Ordered),
                    None => TecplotZoneType::FETriangle,
                },
                _ => TecplotZoneType::Ordered,
            }
        };

        let legacy_format = get(&["F"]).map(str::to_uppercase);
        let packing = match get(&["DATAPACKING"]).map(str::to_uppercase).as_deref() {
            Some("BLOCK") => DataPacking::Block,
            Some("POINT") => DataPacking::Point,
            _ => match legacy_format.as_deref() {
                Some("BLOCK") | Some("FEBLOCK") => DataPacking::Block,
                _ => DataPacking::Point,
            },
        };
        if let Some(location) = get(&["VARLOCATION"])
            && location.to_uppercase().contains("CELLCENTERED")
        {
            return Err(IoError::Parse(format!(
                "Tecplot zone \"{title}\": cell-centered variables (VARLOCATION) are not supported"
            )));
        }

        let mut zone = if zone_type == TecplotZoneType::Ordered {
            TecplotZone::new_ordered(title, i_dim, j_dim, k_dim)
        } else {
            TecplotZone::new_fe(title, zone_type, n_nodes, n_elements, Vec::new())
        };

        // Pre-populate empty variable slots
        for name in var_names {
            zone.variables
                .push(TecplotVariable::new(name.clone(), Vec::new()));
        }

        Ok((zone, packing))
    }

    /// Split a keyword line into its `KEY=value` pairs.
    ///
    /// Keys are bare words directly followed by `=`.  A value is a
    /// double-quoted string (which may contain spaces, commas and `=`), a
    /// parenthesised group such as `DT=(SINGLE SINGLE)`, or an unquoted word
    /// ending at white space or a comma.  Words that are not followed by `=`
    /// (for example the leading `ZONE`) are skipped.
    fn parse_header_pairs(line: &str) -> Vec<(String, String)> {
        let chars: Vec<char> = line.chars().collect();
        let n = chars.len();
        let mut pairs = Vec::new();
        let mut i = 0;
        while i < n {
            if chars[i].is_whitespace() || chars[i] == ',' {
                i += 1;
                continue;
            }
            let key_start = i;
            while i < n && !chars[i].is_whitespace() && !matches!(chars[i], ',' | '=' | '"') {
                i += 1;
            }
            if i == key_start {
                // Stray '=' or '"' that does not belong to a pair.
                i += 1;
                continue;
            }
            let key: String = chars[key_start..i].iter().collect();

            // Optional white space between the key and '='.
            let mut j = i;
            while j < n && chars[j].is_whitespace() {
                j += 1;
            }
            if j >= n || chars[j] != '=' {
                continue;
            }
            j += 1;
            while j < n && chars[j].is_whitespace() {
                j += 1;
            }

            let value: String;
            if j < n && chars[j] == '"' {
                let start = j + 1;
                let mut end = start;
                while end < n && chars[end] != '"' {
                    end += 1;
                }
                value = chars[start..end].iter().collect();
                j = if end < n { end + 1 } else { end };
            } else if j < n && chars[j] == '(' {
                let start = j;
                let mut depth = 0usize;
                let mut end = j;
                while end < n {
                    match chars[end] {
                        '(' => depth += 1,
                        ')' => {
                            depth = depth.saturating_sub(1);
                            if depth == 0 {
                                end += 1;
                                break;
                            }
                        }
                        _ => {}
                    }
                    end += 1;
                }
                value = chars[start..end].iter().collect();
                j = end;
            } else {
                let start = j;
                let mut end = j;
                while end < n && !chars[end].is_whitespace() && chars[end] != ',' {
                    end += 1;
                }
                value = chars[start..end].iter().collect();
                j = end;
            }
            pairs.push((key, value));
            i = j;
        }
        pairs
    }

    /// Whether `line` begins with `keyword`, ignoring ASCII case.  Unlike an
    /// upper-cased copy of the line this allocates nothing, which matters
    /// because it is evaluated for every record of the (possibly huge) file.
    fn starts_with_keyword(line: &str, keyword: &str) -> bool {
        line.as_bytes()
            .get(..keyword.len())
            .is_some_and(|prefix| prefix.eq_ignore_ascii_case(keyword.as_bytes()))
    }

    /// Extract a top-level quoted value: `KEY = "value"`.
    fn extract_quoted_value(line: &str) -> Option<String> {
        let start = line.find('"')?;
        let rest = &line[start + 1..];
        let end = rest.find('"').unwrap_or(rest.len());
        Some(rest[..end].to_string())
    }

    /// Parse a `VARIABLES = "X" "Y" "P"` line into a list of variable names.
    fn parse_variables_line(line: &str) -> Vec<String> {
        // Strip up to and including '='
        let after_eq = if let Some(pos) = line.find('=') {
            &line[pos + 1..]
        } else {
            line
        };
        let mut names = Vec::new();
        let mut chars = after_eq.chars().peekable();
        while let Some(&c) = chars.peek() {
            if c == '"' {
                chars.next(); // consume opening quote
                let mut name = String::new();
                for ch in chars.by_ref() {
                    if ch == '"' {
                        break;
                    }
                    name.push(ch);
                }
                if !name.is_empty() {
                    names.push(name);
                }
            } else {
                chars.next();
            }
        }
        names
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests;
