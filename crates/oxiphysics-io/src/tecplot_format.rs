#![allow(clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tecplot ASCII format I/O.
//!
//! Supports reading and writing Tecplot `.dat` / `.plt` ASCII files with
//! ORDERED and FE (finite element) zone types.  The implementation covers the
//! minimal subset of the Tecplot ASCII specification needed for typical physics
//! simulation output.
//!
//! # Quick start
//!
//! ```no_run
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
//! ds.write("/tmp/quick_start_test.dat").unwrap();
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
    pub fn from_str(s: &str) -> Self {
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
            n_nodes: i_dim * j_dim * k_dim,
            n_elements: 0,
        }
    }

    /// Create a new finite-element zone.
    pub fn new_fe(
        title: impl Into<String>,
        zone_type: TecplotZoneType,
        n_nodes: usize,
        n_elements: usize,
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
        }
    }

    /// Total number of data points expected in this zone.
    pub fn point_count(&self) -> usize {
        match self.zone_type {
            TecplotZoneType::Ordered => self.i_dim * self.j_dim * self.k_dim,
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
            self.write_zone(&mut buf, zone);
        }

        fs::write(path, buf).map_err(IoError::Io)
    }

    /// Append one zone to `buf`.
    fn write_zone(&self, buf: &mut String, zone: &TecplotZone) {
        match zone.zone_type {
            TecplotZoneType::Ordered => {
                writeln!(
                    buf,
                    "ZONE T=\"{}\", I={}, J={}, K={}, DATAPACKING=POINT",
                    zone.title, zone.i_dim, zone.j_dim, zone.k_dim
                )
                .expect("operation should succeed");
                self.write_point_data(buf, zone);
            }
            _ => {
                writeln!(
                    buf,
                    "ZONE T=\"{}\", N={}, E={}, ZONETYPE={}, DATAPACKING=POINT",
                    zone.title,
                    zone.n_nodes,
                    zone.n_elements,
                    zone.zone_type.as_tecplot_str()
                )
                .expect("operation should succeed");
                self.write_point_data(buf, zone);
            }
        }
    }

    /// Write variable data in POINT packing (all vars for node 0, then node 1, …).
    fn write_point_data(&self, buf: &mut String, zone: &TecplotZone) {
        if zone.variables.is_empty() {
            return;
        }
        let n = zone.point_count();
        for idx in 0..n {
            let row: Vec<String> = zone
                .variables
                .iter()
                .map(|v| {
                    if idx < v.data.len() {
                        format!("{:.15e}", v.data[idx])
                    } else {
                        "0.000000000000000e0".to_string()
                    }
                })
                .collect();
            let _ = writeln!(buf, "{}", row.join(" "));
        }
    }
}

// ---------------------------------------------------------------------------
// TecplotReader
// ---------------------------------------------------------------------------

/// Parses a subset of Tecplot ASCII files.
///
/// Handles `TITLE`, `VARIABLES`, and `ZONE` header lines followed by
/// whitespace-delimited floating-point data blocks.
#[derive(Debug, Clone, Default)]
pub struct TecplotReader;

impl TecplotReader {
    /// Create a new reader.
    pub fn new() -> Self {
        Self
    }

    /// Parse the Tecplot ASCII file at `path` and return a `TecplotDataset`.
    pub fn parse(&self, path: &str) -> Result<TecplotDataset, IoError> {
        let file = fs::File::open(path).map_err(IoError::Io)?;
        let reader = io::BufReader::new(file);

        let mut dataset = TecplotDataset::new("");
        let mut current_zone: Option<TecplotZone> = None;
        // How many data points remain for the current zone
        let mut remaining_points: usize = 0;
        // Collected raw floats for the current zone before distributing
        let mut zone_floats: Vec<f64> = Vec::new();

        for line_res in reader.lines() {
            let line = line_res.map_err(IoError::Io)?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let upper = trimmed.to_uppercase();

            // ── TITLE ──────────────────────────────────────────────────────
            if upper.starts_with("TITLE") {
                dataset.title =
                    Self::extract_quoted_value(trimmed).unwrap_or_else(|| trimmed.to_string());
                continue;
            }

            // ── VARIABLES ──────────────────────────────────────────────────
            if upper.starts_with("VARIABLES") {
                dataset.variables = Self::parse_variables_line(trimmed);
                continue;
            }

            // ── ZONE ───────────────────────────────────────────────────────
            if upper.starts_with("ZONE") {
                // Flush previous zone
                if let Some(mut z) = current_zone.take() {
                    Self::distribute_floats(&mut z, &zone_floats);
                    dataset.zones.push(z);
                }
                zone_floats.clear();

                let zone = Self::parse_zone_header(trimmed, &dataset.variables);
                remaining_points = zone.point_count();
                current_zone = Some(zone);
                continue;
            }

            // ── Data line ──────────────────────────────────────────────────
            if current_zone.is_some() && remaining_points > 0 {
                for token in trimmed.split_whitespace() {
                    if let Ok(v) = token.parse::<f64>() {
                        zone_floats.push(v);
                    }
                }
                // Count rows written (POINT packing: one row = one point)
                let nvars = dataset.variables.len().max(1);
                let rows_so_far = zone_floats.len() / nvars;
                if rows_so_far >= remaining_points {
                    remaining_points = 0;
                }
            }
        }

        // Flush last zone
        if let Some(mut z) = current_zone.take() {
            Self::distribute_floats(&mut z, &zone_floats);
            dataset.zones.push(z);
        }

        Ok(dataset)
    }

    /// Distribute a flat float buffer into per-variable data (POINT packing).
    fn distribute_floats(zone: &mut TecplotZone, floats: &[f64]) {
        let nvars = zone.variables.len();
        if nvars == 0 {
            return;
        }
        for (idx, &v) in floats.iter().enumerate() {
            let var_idx = idx % nvars;
            if var_idx < zone.variables.len() {
                zone.variables[var_idx].data.push(v);
            }
        }
    }

    /// Parse a ZONE header line and return a skeleton `TecplotZone`.
    fn parse_zone_header(line: &str, var_names: &[String]) -> TecplotZone {
        let title = Self::extract_param_value(line, "T")
            .or_else(|| Self::extract_param_value(line, "TITLE"))
            .unwrap_or_default();

        let i_dim = Self::extract_param_value(line, "I")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);
        let j_dim = Self::extract_param_value(line, "J")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);
        let k_dim = Self::extract_param_value(line, "K")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);

        let n_nodes = Self::extract_param_value(line, "N")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let n_elements = Self::extract_param_value(line, "E")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let zone_type_str = Self::extract_param_value(line, "ZONETYPE")
            .or_else(|| Self::extract_param_value(line, "F"))
            .unwrap_or_else(|| "ORDERED".to_string());
        let zone_type = TecplotZoneType::from_str(&zone_type_str);

        let mut zone = if zone_type == TecplotZoneType::Ordered {
            TecplotZone::new_ordered(title, i_dim, j_dim, k_dim)
        } else {
            TecplotZone::new_fe(title, zone_type, n_nodes, n_elements)
        };

        // Pre-populate empty variable slots
        for name in var_names {
            zone.variables
                .push(TecplotVariable::new(name.clone(), Vec::new()));
        }

        zone
    }

    /// Extract the value of `KEY=value` from a keyword line (case-insensitive).
    ///
    /// Handles quoted values, e.g. `T="My Zone"`.
    fn extract_param_value(line: &str, key: &str) -> Option<String> {
        let upper = line.to_uppercase();
        let key_eq = format!("{key}=");
        // Search case-insensitively
        let pos = upper.find(&key_eq.to_uppercase())?;
        let rest = &line[pos + key_eq.len()..];
        let rest_trimmed = rest.trim_start();
        if let Some(inner) = rest_trimmed.strip_prefix('"') {
            // Quoted value
            let end = inner.find('"').unwrap_or(inner.len());
            Some(inner[..end].to_string())
        } else {
            // Unquoted: ends at next comma, space, or end of string
            let end = rest_trimmed.find([',', ' ']).unwrap_or(rest_trimmed.len());
            Some(rest_trimmed[..end].trim().to_string())
        }
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
mod tests {
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
            TecplotZoneType::from_str("ORDERED"),
            TecplotZoneType::Ordered
        );
    }

    #[test]
    fn test_zone_type_from_str_fe_triangle() {
        assert_eq!(
            TecplotZoneType::from_str("FETRIANGLE"),
            TecplotZoneType::FETriangle
        );
    }

    #[test]
    fn test_zone_type_from_str_fe_tetra() {
        assert_eq!(
            TecplotZoneType::from_str("FETETRAHEDRON"),
            TecplotZoneType::FETetra
        );
    }

    #[test]
    fn test_zone_type_from_str_case_insensitive() {
        assert_eq!(
            TecplotZoneType::from_str("ordered"),
            TecplotZoneType::Ordered
        );
    }

    #[test]
    fn test_zone_type_from_str_unknown_defaults_to_ordered() {
        assert_eq!(
            TecplotZoneType::from_str("UNKNOWN"),
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
        let z = TecplotZone::new_fe("FE Zone", TecplotZoneType::FETetra, 8, 2);
        assert_eq!(z.zone_type, TecplotZoneType::FETetra);
        assert_eq!(z.n_nodes, 8);
        assert_eq!(z.n_elements, 2);
        assert_eq!(z.point_count(), 8);
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
        let path = "/tmp/oxiphysics_tecplot_write_test.dat";
        let ds = make_1d_dataset(path, 5);
        TecplotWriter::new().write(&ds, path).expect("write failed");
        assert!(std::path::Path::new(path).exists());
    }

    #[test]
    fn test_write_contains_title() {
        let path = "/tmp/oxiphysics_tecplot_title.dat";
        let ds = make_1d_dataset(path, 3);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("Test Dataset"), "no title in output");
    }

    #[test]
    fn test_write_contains_variables() {
        let path = "/tmp/oxiphysics_tecplot_vars.dat";
        let ds = make_1d_dataset(path, 3);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("VARIABLES"), "no VARIABLES header");
        assert!(content.contains('"'), "variables not quoted");
    }

    #[test]
    fn test_write_contains_zone_keyword() {
        let path = "/tmp/oxiphysics_tecplot_zone_kw.dat";
        let ds = make_1d_dataset(path, 3);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("ZONE"), "no ZONE header");
    }

    #[test]
    fn test_write_ordered_has_ijk() {
        let path = "/tmp/oxiphysics_tecplot_ijk.dat";
        let ds = make_1d_dataset(path, 3);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("I="), "no I= in ZONE header");
        assert!(content.contains("J="), "no J= in ZONE header");
        assert!(content.contains("K="), "no K= in ZONE header");
    }

    #[test]
    fn test_write_data_lines_count() {
        // A 1D ordered zone with I=4 should produce exactly 4 data rows
        let path = "/tmp/oxiphysics_tecplot_datarows.dat";
        let ds = make_1d_dataset(path, 4);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
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
        let path = "/tmp/oxiphysics_tecplot_fe_zone.dat";
        let mut ds = TecplotDataset::new("FE Dataset");
        ds.variables = vec!["X".to_string()];
        let mut zone = TecplotZone::new_fe("FE Zone", TecplotZoneType::FETetra, 4, 1);
        zone.variables
            .push(TecplotVariable::new("X", vec![0.0, 1.0, 0.0, 0.0]));
        ds.zones.push(zone);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("FETETRAHEDRON"));
    }

    // ── Write/Read roundtrip ────────────────────────────────────────────────

    #[test]
    fn test_roundtrip_title() {
        let path = "/tmp/oxiphysics_tecplot_rt_title.dat";
        let ds = make_1d_dataset(path, 3);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        assert_eq!(parsed.title, "Test Dataset");
    }

    #[test]
    fn test_roundtrip_variable_names() {
        let path = "/tmp/oxiphysics_tecplot_rt_varnames.dat";
        let ds = make_1d_dataset(path, 3);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        assert_eq!(parsed.variables, vec!["X", "P"]);
    }

    #[test]
    fn test_roundtrip_zone_count() {
        let path = "/tmp/oxiphysics_tecplot_rt_nzones.dat";
        let ds = make_1d_dataset(path, 3);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        assert_eq!(parsed.zones.len(), 1, "expected 1 zone");
    }

    #[test]
    fn test_roundtrip_zone_type_ordered() {
        let path = "/tmp/oxiphysics_tecplot_rt_ztype.dat";
        let ds = make_1d_dataset(path, 3);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        assert_eq!(parsed.zones[0].zone_type, TecplotZoneType::Ordered);
    }

    #[test]
    fn test_roundtrip_zone_i_dim() {
        let path = "/tmp/oxiphysics_tecplot_rt_idim.dat";
        let ds = make_1d_dataset(path, 5);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        assert_eq!(parsed.zones[0].i_dim, 5);
    }

    #[test]
    fn test_roundtrip_data_values() {
        let path = "/tmp/oxiphysics_tecplot_rt_data.dat";
        let ds = make_1d_dataset(path, 3);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        let zone = &parsed.zones[0];
        // X variable should be [0,1,2]
        let x_var = zone.variables.iter().find(|v| v.name == "X").unwrap();
        assert!((x_var.data[0]).abs() < 1e-6);
        assert!((x_var.data[1] - 1.0).abs() < 1e-6);
        assert!((x_var.data[2] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_roundtrip_second_variable() {
        let path = "/tmp/oxiphysics_tecplot_rt_p.dat";
        let ds = make_1d_dataset(path, 3);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        let zone = &parsed.zones[0];
        let p_var = zone.variables.iter().find(|v| v.name == "P").unwrap();
        assert!((p_var.data[0]).abs() < 1e-6, "P[0]={}", p_var.data[0]);
        assert!((p_var.data[1] - 2.0).abs() < 1e-6, "P[1]={}", p_var.data[1]);
        assert!((p_var.data[2] - 4.0).abs() < 1e-6, "P[2]={}", p_var.data[2]);
    }

    #[test]
    fn test_roundtrip_multi_zone() {
        let path = "/tmp/oxiphysics_tecplot_rt_mz.dat";
        let mut ds = TecplotDataset::new("Multi");
        ds.variables = vec!["X".to_string()];
        for i in 0..3usize {
            let mut z = TecplotZone::new_ordered(format!("Zone {i}"), 2, 1, 1);
            z.variables
                .push(TecplotVariable::new("X", vec![i as f64, i as f64 + 1.0]));
            ds.zones.push(z);
        }
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        assert_eq!(parsed.zones.len(), 3, "expected 3 zones");
    }

    #[test]
    fn test_read_missing_file_returns_error() {
        let result = TecplotDataset::read("/tmp/no_such_file_oxiphysics_tec.dat");
        assert!(result.is_err());
    }

    #[test]
    fn test_roundtrip_large_ordered_zone() {
        let path = "/tmp/oxiphysics_tecplot_large.dat";
        let n = 50;
        let ds = make_1d_dataset(path, n);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        let zone = &parsed.zones[0];
        assert_eq!(zone.i_dim, n);
        let x_var = zone.variables.iter().find(|v| v.name == "X").unwrap();
        assert_eq!(x_var.data.len(), n);
    }

    #[test]
    fn test_roundtrip_3d_zone() {
        let path = "/tmp/oxiphysics_tecplot_3d.dat";
        let mut ds = TecplotDataset::new("3D");
        ds.variables = vec!["X".to_string()];
        let n = 2usize;
        let total = n * n * n;
        let mut zone = TecplotZone::new_ordered("3D Zone", n, n, n);
        let x: Vec<f64> = (0..total).map(|i| i as f64).collect();
        zone.variables.push(TecplotVariable::new("X", x));
        ds.zones.push(zone);
        ds.write(path).unwrap();
        let parsed = TecplotDataset::read(path).unwrap();
        assert_eq!(parsed.zones[0].i_dim, n);
        assert_eq!(parsed.zones[0].j_dim, n);
        assert_eq!(parsed.zones[0].k_dim, n);
    }

    // ── TecplotWriter helpers ────────────────────────────────────────────────

    #[test]
    fn test_writer_datapacking_point_keyword() {
        let path = "/tmp/oxiphysics_tecplot_dpkw.dat";
        let ds = make_1d_dataset(path, 2);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
        assert!(content.contains("DATAPACKING=POINT"));
    }

    #[test]
    fn test_writer_zone_title_in_output() {
        let path = "/tmp/oxiphysics_tecplot_ztitle.dat";
        let ds = make_1d_dataset(path, 2);
        TecplotWriter::new().write(&ds, path).unwrap();
        let content = fs::read_to_string(path).unwrap();
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
}
