#![allow(clippy::manual_strip)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! EnSight Gold format reader/writer for physics simulation data.
//!
//! Supports ASCII EnSight Gold case files, geometry files, and per-node
//! scalar and vector variable files.  Covers unstructured element types:
//! point, tria3 (triangle), tetra4 (tetrahedron), and hexa8 (hexahedron).

use std::fmt;
#[cfg(test)]
use std::io::BufReader;
use std::io::{self, BufRead, Write};
use std::path::Path;

// ---------------------------------------------------------------------------
// NodeIdMode
// ---------------------------------------------------------------------------

/// EnSight node ID mode, written in the geometry file header.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NodeIdMode {
    /// No node IDs stored — EnSight assigns them internally.
    #[default]
    Off,
    /// Node IDs are given explicitly in the file.
    Given,
    /// EnSight assigns node IDs automatically.
    Assign,
}

impl fmt::Display for NodeIdMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NodeIdMode::Off => write!(f, "off"),
            NodeIdMode::Given => write!(f, "given"),
            NodeIdMode::Assign => write!(f, "assign"),
        }
    }
}

// ---------------------------------------------------------------------------
// EnsightPart
// ---------------------------------------------------------------------------

/// A named part within an EnSight geometry: holds element type + connectivity.
///
/// Connectivity is stored as a flat `Vec`u32` where each group of
/// `nodes_per_element()` values represents one element.
#[derive(Debug, Clone)]
pub struct EnsightPart {
    /// Part number (1-based in EnSight convention).
    pub part_id: u32,
    /// Human-readable part description.
    pub description: String,
    /// Element type keyword (e.g. `"point"`, `"tria3"`, `"tetra4"`, `"hexa8"`).
    pub element_type: String,
    /// Flat connectivity list: `nodes_per_element * n_elements` entries (0-based node indices).
    pub connectivity: Vec<u32>,
}

impl EnsightPart {
    /// Create a new part.
    pub fn new(
        part_id: u32,
        description: impl Into<String>,
        element_type: impl Into<String>,
    ) -> Self {
        Self {
            part_id,
            description: description.into(),
            element_type: element_type.into(),
            connectivity: Vec::new(),
        }
    }

    /// Number of nodes per element for this element type.
    pub fn nodes_per_element(&self) -> usize {
        match self.element_type.as_str() {
            "point" => 1,
            "bar2" => 2,
            "tria3" => 3,
            "quad4" => 4,
            "tetra4" => 4,
            "pyramid5" => 5,
            "penta6" => 6,
            "hexa8" => 8,
            _ => 0,
        }
    }

    /// Number of elements stored in this part.
    pub fn n_elements(&self) -> usize {
        let npe = self.nodes_per_element();
        if npe == 0 {
            return 0;
        }
        self.connectivity.len() / npe
    }
}

// ---------------------------------------------------------------------------
// EnsightGeometry
// ---------------------------------------------------------------------------

/// EnSight Gold geometry: node coordinates + one or more parts.
#[derive(Debug, Clone)]
pub struct EnsightGeometry {
    /// Description line 1 (appears in .geo header).
    pub description1: String,
    /// Description line 2 (appears in .geo header).
    pub description2: String,
    /// Node ID mode.
    pub node_id_mode: NodeIdMode,
    /// Element ID mode (`"off"` | `"given"` | `"assign"`).
    pub element_id_mode: NodeIdMode,
    /// X coordinates of all nodes (0-based).
    pub x: Vec<f32>,
    /// Y coordinates of all nodes (0-based).
    pub y: Vec<f32>,
    /// Z coordinates of all nodes (0-based).
    pub z: Vec<f32>,
    /// Parts in this geometry.
    pub parts: Vec<EnsightPart>,
}

impl EnsightGeometry {
    /// Create an empty geometry with default metadata.
    pub fn new() -> Self {
        Self {
            description1: String::from("EnSight Gold Geometry"),
            description2: String::from("Created by OxiPhysics"),
            node_id_mode: NodeIdMode::Off,
            element_id_mode: NodeIdMode::Off,
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            parts: Vec::new(),
        }
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.x.len()
    }

    /// Add a node, returning its 0-based index.
    pub fn add_node(&mut self, xi: f32, yi: f32, zi: f32) -> u32 {
        let idx = self.x.len() as u32;
        self.x.push(xi);
        self.y.push(yi);
        self.z.push(zi);
        idx
    }

    /// Add a part.
    pub fn add_part(&mut self, part: EnsightPart) {
        self.parts.push(part);
    }
}

impl Default for EnsightGeometry {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// EnsightScalarVar
// ---------------------------------------------------------------------------

/// Per-node scalar variable (e.g. pressure, temperature).
#[derive(Debug, Clone)]
pub struct EnsightScalarVar {
    /// Variable name (used in .case file and as filename stem).
    pub name: String,
    /// Per-node scalar values, length == number of nodes.
    pub values: Vec<f32>,
}

impl EnsightScalarVar {
    /// Create a new scalar variable.
    pub fn new(name: impl Into<String>, values: Vec<f32>) -> Self {
        Self {
            name: name.into(),
            values,
        }
    }
}

// ---------------------------------------------------------------------------
// EnsightVectorVar
// ---------------------------------------------------------------------------

/// Per-node vector variable (e.g. velocity, displacement) with 3 components.
#[derive(Debug, Clone)]
pub struct EnsightVectorVar {
    /// Variable name (used in .case file and as filename stem).
    pub name: String,
    /// X components, length == number of nodes.
    pub vx: Vec<f32>,
    /// Y components, length == number of nodes.
    pub vy: Vec<f32>,
    /// Z components, length == number of nodes.
    pub vz: Vec<f32>,
}

impl EnsightVectorVar {
    /// Create a new vector variable.
    pub fn new(name: impl Into<String>, vx: Vec<f32>, vy: Vec<f32>, vz: Vec<f32>) -> Self {
        Self {
            name: name.into(),
            vx,
            vy,
            vz,
        }
    }

    /// Number of nodes this variable covers.
    pub fn n_nodes(&self) -> usize {
        self.vx.len()
    }
}

// ---------------------------------------------------------------------------
// EnsightCase
// ---------------------------------------------------------------------------

/// EnSight Gold case file descriptor.
#[derive(Debug, Clone)]
pub struct EnsightCase {
    /// Path to the geometry file (relative to the case file directory).
    pub geometry_file: String,
    /// Scalar variable file paths (parallel to variable names).
    pub scalar_files: Vec<(String, String)>,
    /// Vector variable file paths (parallel to variable names).
    pub vector_files: Vec<(String, String)>,
}

impl EnsightCase {
    /// Create an empty case.
    pub fn new(geometry_file: impl Into<String>) -> Self {
        Self {
            geometry_file: geometry_file.into(),
            scalar_files: Vec::new(),
            vector_files: Vec::new(),
        }
    }

    /// Register a scalar variable file.
    pub fn add_scalar(&mut self, name: impl Into<String>, file: impl Into<String>) {
        self.scalar_files.push((name.into(), file.into()));
    }

    /// Register a vector variable file.
    pub fn add_vector(&mut self, name: impl Into<String>, file: impl Into<String>) {
        self.vector_files.push((name.into(), file.into()));
    }

    /// Write the `.case` file contents to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "FORMAT")?;
        writeln!(writer, "type: ensight gold")?;
        writeln!(writer)?;
        writeln!(writer, "GEOMETRY")?;
        writeln!(writer, "model: {}", self.geometry_file)?;
        if !self.scalar_files.is_empty() || !self.vector_files.is_empty() {
            writeln!(writer)?;
            writeln!(writer, "VARIABLE")?;
            for (name, file) in &self.scalar_files {
                writeln!(writer, "scalar per node: {name} {file}")?;
            }
            for (name, file) in &self.vector_files {
                writeln!(writer, "vector per node: {name} {file}")?;
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// EnsightTimeSeries
// ---------------------------------------------------------------------------

/// Time-varying EnSight case: a sequence of geometry/variable file sets.
#[derive(Debug, Clone)]
pub struct EnsightTimeSeries {
    /// Time values for each step.
    pub times: Vec<f64>,
    /// Geometry files per time step.
    pub geo_files: Vec<String>,
}

impl EnsightTimeSeries {
    /// Create an empty time series.
    pub fn new() -> Self {
        Self {
            times: Vec::new(),
            geo_files: Vec::new(),
        }
    }

    /// Add a time step.
    pub fn push(&mut self, time: f64, geo_file: impl Into<String>) {
        self.times.push(time);
        self.geo_files.push(geo_file.into());
    }

    /// Number of time steps.
    pub fn n_steps(&self) -> usize {
        self.times.len()
    }

    /// Write a transient `.case` file section to a writer.
    pub fn write_transient_case<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writeln!(writer, "FORMAT")?;
        writeln!(writer, "type: ensight gold")?;
        writeln!(writer)?;
        writeln!(writer, "GEOMETRY")?;
        writeln!(
            writer,
            "model: 1 {}",
            if self.geo_files.is_empty() {
                "geometry"
            } else {
                &self.geo_files[0]
            }
        )?;
        writeln!(writer)?;
        writeln!(writer, "TIME")?;
        writeln!(writer, "time set: 1")?;
        writeln!(writer, "number of steps: {}", self.times.len())?;
        writeln!(writer, "time values:")?;
        for &t in &self.times {
            writeln!(writer, "  {t:.6e}")?;
        }
        Ok(())
    }
}

impl Default for EnsightTimeSeries {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// EnsightWriter
// ---------------------------------------------------------------------------

/// ASCII EnSight Gold writer: writes `.case`, `.geo`, and variable files.
pub struct EnsightWriter;

impl EnsightWriter {
    /// Write an EnSight Gold geometry file (ASCII) to a writer.
    pub fn write_geo<W: Write>(writer: &mut W, geo: &EnsightGeometry) -> io::Result<()> {
        writeln!(writer, "{}", geo.description1)?;
        writeln!(writer, "{}", geo.description2)?;
        writeln!(writer, "node id {}", geo.node_id_mode)?;
        writeln!(writer, "element id {}", geo.element_id_mode)?;

        for part in &geo.parts {
            writeln!(writer, "part")?;
            writeln!(writer, "{:10}", part.part_id)?;
            writeln!(writer, "{}", part.description)?;
            writeln!(writer, "coordinates")?;
            writeln!(writer, "{:10}", geo.n_nodes())?;
            for &v in &geo.x {
                writeln!(writer, "{v:12.5e}")?;
            }
            for &v in &geo.y {
                writeln!(writer, "{v:12.5e}")?;
            }
            for &v in &geo.z {
                writeln!(writer, "{v:12.5e}")?;
            }
            writeln!(writer, "{}", part.element_type)?;
            writeln!(writer, "{:10}", part.n_elements())?;
            let npe = part.nodes_per_element();
            for chunk in part.connectivity.chunks(npe) {
                let s: Vec<String> = chunk.iter().map(|&c| format!("{:10}", c + 1)).collect();
                writeln!(writer, "{}", s.join(""))?;
            }
        }
        Ok(())
    }

    /// Write a scalar variable file (ASCII EnSight Gold).
    pub fn write_scalar<W: Write>(writer: &mut W, var: &EnsightScalarVar) -> io::Result<()> {
        writeln!(writer, "{}", var.name)?;
        for part_id in 1u32..=1 {
            writeln!(writer, "part")?;
            writeln!(writer, "{part_id:10}")?;
            writeln!(writer, "coordinates")?;
            for &v in &var.values {
                writeln!(writer, "{v:12.5e}")?;
            }
        }
        Ok(())
    }

    /// Write a vector variable file (ASCII EnSight Gold).
    pub fn write_vector<W: Write>(writer: &mut W, var: &EnsightVectorVar) -> io::Result<()> {
        writeln!(writer, "{}", var.name)?;
        for part_id in 1u32..=1 {
            writeln!(writer, "part")?;
            writeln!(writer, "{part_id:10}")?;
            writeln!(writer, "coordinates")?;
            for i in 0..var.vx.len() {
                writeln!(writer, "{:12.5e}", var.vx[i])?;
            }
            for i in 0..var.vy.len() {
                writeln!(writer, "{:12.5e}", var.vy[i])?;
            }
            for i in 0..var.vz.len() {
                writeln!(writer, "{:12.5e}", var.vz[i])?;
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// EnsightReader
// ---------------------------------------------------------------------------

/// ASCII EnSight Gold reader.
pub struct EnsightReader;

impl EnsightReader {
    /// Parse an EnSight Gold case file from a reader.
    ///
    /// Returns an [`EnsightCase`] with paths filled in but no data loaded yet.
    pub fn read_case<R: BufRead>(reader: R) -> io::Result<EnsightCase> {
        let mut geo_file = String::new();
        let mut case = EnsightCase::new("");

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.starts_with("model:") {
                let rest = trimmed["model:".len()..].trim();
                // May be "1 filename" (transient) or just "filename"
                let parts: Vec<&str> = rest.splitn(2, ' ').collect();
                geo_file = parts.last().unwrap_or(&"").to_string();
            } else if trimmed.starts_with("scalar per node:") {
                let rest = trimmed["scalar per node:".len()..].trim();
                let parts: Vec<&str> = rest.splitn(2, ' ').collect();
                if parts.len() == 2 {
                    case.scalar_files
                        .push((parts[0].to_string(), parts[1].to_string()));
                }
            } else if trimmed.starts_with("vector per node:") {
                let rest = trimmed["vector per node:".len()..].trim();
                let parts: Vec<&str> = rest.splitn(2, ' ').collect();
                if parts.len() == 2 {
                    case.vector_files
                        .push((parts[0].to_string(), parts[1].to_string()));
                }
            }
        }
        case.geometry_file = geo_file;
        Ok(case)
    }

    /// Parse an EnSight Gold geometry file from a reader.
    ///
    /// Returns an [`EnsightGeometry`] with nodes and parts populated.
    pub fn read_geo<R: BufRead>(reader: R) -> io::Result<EnsightGeometry> {
        let mut geo = EnsightGeometry::new();
        let mut lines = reader.lines();

        // Header lines
        if let Some(l) = lines.next() {
            geo.description1 = l?.trim().to_string();
        }
        if let Some(l) = lines.next() {
            geo.description2 = l?.trim().to_string();
        }
        if let Some(l) = lines.next() {
            let s = l?;
            if s.contains("given") {
                geo.node_id_mode = NodeIdMode::Given;
            } else if s.contains("assign") {
                geo.node_id_mode = NodeIdMode::Assign;
            }
        }
        if let Some(l) = lines.next() {
            let s = l?;
            if s.contains("given") {
                geo.element_id_mode = NodeIdMode::Given;
            } else if s.contains("assign") {
                geo.element_id_mode = NodeIdMode::Assign;
            }
        }

        let mut line_buf: Vec<String> = lines.map(|l| l.unwrap_or_default()).collect();
        let mut pos = 0;

        while pos < line_buf.len() {
            let line = line_buf[pos].trim().to_string();
            pos += 1;

            if line == "part" {
                let part_id: u32 = line_buf
                    .get(pos)
                    .map(|l| l.trim().parse().unwrap_or(1))
                    .unwrap_or(1);
                pos += 1;
                let description = line_buf
                    .get(pos)
                    .map(|l| l.trim().to_string())
                    .unwrap_or_default();
                pos += 1;
                let keyword = line_buf
                    .get(pos)
                    .map(|l| l.trim().to_string())
                    .unwrap_or_default();
                pos += 1;

                if keyword == "coordinates" {
                    let n_nodes: usize = line_buf
                        .get(pos)
                        .map(|l| l.trim().parse().unwrap_or(0))
                        .unwrap_or(0);
                    pos += 1;
                    let mut xs = Vec::with_capacity(n_nodes);
                    let mut ys = Vec::with_capacity(n_nodes);
                    let mut zs = Vec::with_capacity(n_nodes);
                    for _ in 0..n_nodes {
                        let v: f32 = line_buf
                            .get(pos)
                            .map(|l| l.trim().parse().unwrap_or(0.0))
                            .unwrap_or(0.0);
                        xs.push(v);
                        pos += 1;
                    }
                    for _ in 0..n_nodes {
                        let v: f32 = line_buf
                            .get(pos)
                            .map(|l| l.trim().parse().unwrap_or(0.0))
                            .unwrap_or(0.0);
                        ys.push(v);
                        pos += 1;
                    }
                    for _ in 0..n_nodes {
                        let v: f32 = line_buf
                            .get(pos)
                            .map(|l| l.trim().parse().unwrap_or(0.0))
                            .unwrap_or(0.0);
                        zs.push(v);
                        pos += 1;
                    }
                    geo.x = xs;
                    geo.y = ys;
                    geo.z = zs;

                    // Next: element type
                    let elem_type = line_buf
                        .get(pos)
                        .map(|l| l.trim().to_string())
                        .unwrap_or_default();
                    pos += 1;
                    let n_elements: usize = line_buf
                        .get(pos)
                        .map(|l| l.trim().parse().unwrap_or(0))
                        .unwrap_or(0);
                    pos += 1;

                    let mut part = EnsightPart::new(part_id, description, &elem_type);
                    let npe = part.nodes_per_element();
                    for _ in 0..n_elements {
                        let row = line_buf
                            .get(pos)
                            .map(|l| l.trim().to_string())
                            .unwrap_or_default();
                        pos += 1;
                        // Parse space-separated or fixed-width (10-char) node ids
                        let tokens: Vec<u32> = row
                            .split_whitespace()
                            .filter_map(|t| t.parse::<u32>().ok())
                            .map(|v| v - 1) // convert to 0-based
                            .collect();
                        for k in 0..npe {
                            part.connectivity.push(*tokens.get(k).unwrap_or(&0));
                        }
                    }
                    geo.parts.push(part);
                }
            }
        }
        line_buf.clear();
        Ok(geo)
    }
}

// ---------------------------------------------------------------------------
// Convenience write function
// ---------------------------------------------------------------------------

/// Write a complete EnSight Gold dataset to disk (ASCII).
///
/// Creates `path`.case`, `path`.geo`, and per-variable files.
/// All files are placed in the directory of `path` (the stem is used as a
/// base name).
pub fn write_ensight_case(
    path: &str,
    geo: &EnsightGeometry,
    scalars: &[EnsightScalarVar],
    vectors: &[EnsightVectorVar],
) -> io::Result<()> {
    let stem = Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");
    let dir = Path::new(path).parent().unwrap_or(Path::new("."));

    let geo_name = format!("{stem}.geo");
    let case_name = format!("{stem}.case");

    // Write geometry file
    let geo_path = dir.join(&geo_name);
    let mut geo_file = std::fs::File::create(&geo_path)?;
    EnsightWriter::write_geo(&mut geo_file, geo)?;

    // Build case
    let mut case = EnsightCase::new(&geo_name);
    let mut scalar_bufs: Vec<Vec<u8>> = Vec::new();
    let mut vector_bufs: Vec<Vec<u8>> = Vec::new();

    for sv in scalars {
        let fname = format!("{stem}_{}.escl", sv.name);
        let fpath = dir.join(&fname);
        let mut buf = Vec::new();
        EnsightWriter::write_scalar(&mut buf, sv)?;
        scalar_bufs.push(buf.clone());
        std::fs::write(&fpath, &buf)?;
        case.add_scalar(&sv.name, &fname);
    }

    for vv in vectors {
        let fname = format!("{stem}_{}.evec", vv.name);
        let fpath = dir.join(&fname);
        let mut buf = Vec::new();
        EnsightWriter::write_vector(&mut buf, vv)?;
        vector_bufs.push(buf.clone());
        std::fs::write(&fpath, &buf)?;
        case.add_vector(&vv.name, &fname);
    }

    // Write case file
    let case_path = dir.join(&case_name);
    let mut case_file = std::fs::File::create(&case_path)?;
    case.write_to(&mut case_file)?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── NodeIdMode ────────────────────────────────────────────────────────

    #[test]
    fn test_node_id_mode_display_off() {
        assert_eq!(NodeIdMode::Off.to_string(), "off");
    }

    #[test]
    fn test_node_id_mode_display_given() {
        assert_eq!(NodeIdMode::Given.to_string(), "given");
    }

    #[test]
    fn test_node_id_mode_display_assign() {
        assert_eq!(NodeIdMode::Assign.to_string(), "assign");
    }

    #[test]
    fn test_node_id_mode_default() {
        assert_eq!(NodeIdMode::default(), NodeIdMode::Off);
    }

    #[test]
    fn test_node_id_mode_equality() {
        assert_eq!(NodeIdMode::Given, NodeIdMode::Given);
        assert_ne!(NodeIdMode::Given, NodeIdMode::Off);
    }

    // ── EnsightPart ───────────────────────────────────────────────────────

    #[test]
    fn test_part_nodes_per_element_point() {
        let p = EnsightPart::new(1, "test", "point");
        assert_eq!(p.nodes_per_element(), 1);
    }

    #[test]
    fn test_part_nodes_per_element_tria3() {
        let p = EnsightPart::new(1, "test", "tria3");
        assert_eq!(p.nodes_per_element(), 3);
    }

    #[test]
    fn test_part_nodes_per_element_tetra4() {
        let p = EnsightPart::new(1, "test", "tetra4");
        assert_eq!(p.nodes_per_element(), 4);
    }

    #[test]
    fn test_part_nodes_per_element_hexa8() {
        let p = EnsightPart::new(1, "test", "hexa8");
        assert_eq!(p.nodes_per_element(), 8);
    }

    #[test]
    fn test_part_n_elements_tria3() {
        let mut p = EnsightPart::new(1, "mesh", "tria3");
        p.connectivity = vec![0, 1, 2, 1, 2, 3];
        assert_eq!(p.n_elements(), 2);
    }

    #[test]
    fn test_part_n_elements_hexa8() {
        let mut p = EnsightPart::new(1, "vol", "hexa8");
        p.connectivity = vec![0; 8];
        assert_eq!(p.n_elements(), 1);
    }

    #[test]
    fn test_part_n_elements_empty() {
        let p = EnsightPart::new(1, "empty", "tetra4");
        assert_eq!(p.n_elements(), 0);
    }

    // ── EnsightGeometry ───────────────────────────────────────────────────

    #[test]
    fn test_geo_new_empty() {
        let geo = EnsightGeometry::new();
        assert_eq!(geo.n_nodes(), 0);
        assert!(geo.parts.is_empty());
    }

    #[test]
    fn test_geo_add_node() {
        let mut geo = EnsightGeometry::new();
        let idx = geo.add_node(1.0, 2.0, 3.0);
        assert_eq!(idx, 0);
        assert_eq!(geo.n_nodes(), 1);
        assert!((geo.x[0] - 1.0).abs() < 1e-7);
        assert!((geo.y[0] - 2.0).abs() < 1e-7);
        assert!((geo.z[0] - 3.0).abs() < 1e-7);
    }

    #[test]
    fn test_geo_add_multiple_nodes() {
        let mut geo = EnsightGeometry::new();
        geo.add_node(0.0, 0.0, 0.0);
        geo.add_node(1.0, 0.0, 0.0);
        geo.add_node(0.0, 1.0, 0.0);
        assert_eq!(geo.n_nodes(), 3);
    }

    #[test]
    fn test_geo_add_part() {
        let mut geo = EnsightGeometry::new();
        let part = EnsightPart::new(1, "surf", "tria3");
        geo.add_part(part);
        assert_eq!(geo.parts.len(), 1);
    }

    #[test]
    fn test_geo_default() {
        let geo = EnsightGeometry::default();
        assert_eq!(geo.n_nodes(), 0);
    }

    // ── EnsightScalarVar ──────────────────────────────────────────────────

    #[test]
    fn test_scalar_var_new() {
        let sv = EnsightScalarVar::new("pressure", vec![1.0, 2.0, 3.0]);
        assert_eq!(sv.name, "pressure");
        assert_eq!(sv.values.len(), 3);
    }

    #[test]
    fn test_scalar_var_values() {
        let sv = EnsightScalarVar::new("temp", vec![300.0, 350.0]);
        assert!((sv.values[0] - 300.0).abs() < 1e-6);
        assert!((sv.values[1] - 350.0).abs() < 1e-6);
    }

    // ── EnsightVectorVar ──────────────────────────────────────────────────

    #[test]
    fn test_vector_var_new() {
        let vv = EnsightVectorVar::new("velocity", vec![1.0, 2.0], vec![0.0, 0.0], vec![0.0, 0.0]);
        assert_eq!(vv.name, "velocity");
        assert_eq!(vv.n_nodes(), 2);
    }

    #[test]
    fn test_vector_var_n_nodes() {
        let vv = EnsightVectorVar::new("disp", vec![0.1; 5], vec![0.2; 5], vec![0.3; 5]);
        assert_eq!(vv.n_nodes(), 5);
    }

    // ── EnsightCase ───────────────────────────────────────────────────────

    #[test]
    fn test_case_new() {
        let c = EnsightCase::new("out.geo");
        assert_eq!(c.geometry_file, "out.geo");
        assert!(c.scalar_files.is_empty());
        assert!(c.vector_files.is_empty());
    }

    #[test]
    fn test_case_add_scalar() {
        let mut c = EnsightCase::new("out.geo");
        c.add_scalar("pressure", "out_pressure.escl");
        assert_eq!(c.scalar_files.len(), 1);
        assert_eq!(c.scalar_files[0].0, "pressure");
    }

    #[test]
    fn test_case_add_vector() {
        let mut c = EnsightCase::new("out.geo");
        c.add_vector("velocity", "out_velocity.evec");
        assert_eq!(c.vector_files.len(), 1);
    }

    #[test]
    fn test_case_write_to_contains_format() {
        let c = EnsightCase::new("geom.geo");
        let mut buf = Vec::new();
        c.write_to(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("FORMAT"), "missing FORMAT section");
        assert!(s.contains("type: ensight gold"));
        assert!(s.contains("geom.geo"));
    }

    #[test]
    fn test_case_write_to_with_scalar() {
        let mut c = EnsightCase::new("g.geo");
        c.add_scalar("p", "p.escl");
        let mut buf = Vec::new();
        c.write_to(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("scalar per node: p p.escl"), "output: {s}");
    }

    #[test]
    fn test_case_write_to_with_vector() {
        let mut c = EnsightCase::new("g.geo");
        c.add_vector("vel", "vel.evec");
        let mut buf = Vec::new();
        c.write_to(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("vector per node: vel vel.evec"), "output: {s}");
    }

    // ── EnsightWriter geo ─────────────────────────────────────────────────

    #[test]
    fn test_writer_geo_single_point_part() {
        let mut geo = EnsightGeometry::new();
        geo.add_node(0.0, 0.0, 0.0);
        let mut part = EnsightPart::new(1, "pts", "point");
        part.connectivity = vec![0];
        geo.add_part(part);

        let mut buf = Vec::new();
        EnsightWriter::write_geo(&mut buf, &geo).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("part"), "missing 'part' keyword");
        assert!(s.contains("coordinates"));
        assert!(s.contains("point"));
    }

    #[test]
    fn test_writer_geo_tria3_part() {
        let mut geo = EnsightGeometry::new();
        geo.add_node(0.0, 0.0, 0.0);
        geo.add_node(1.0, 0.0, 0.0);
        geo.add_node(0.0, 1.0, 0.0);
        let mut part = EnsightPart::new(1, "surf", "tria3");
        part.connectivity = vec![0, 1, 2];
        geo.add_part(part);

        let mut buf = Vec::new();
        EnsightWriter::write_geo(&mut buf, &geo).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("tria3"));
    }

    // ── EnsightWriter scalar ──────────────────────────────────────────────

    #[test]
    fn test_writer_scalar_output() {
        let sv = EnsightScalarVar::new("pressure", vec![1.5, 2.5, 3.5]);
        let mut buf = Vec::new();
        EnsightWriter::write_scalar(&mut buf, &sv).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("pressure"));
        assert!(s.contains("coordinates"));
    }

    // ── EnsightWriter vector ──────────────────────────────────────────────

    #[test]
    fn test_writer_vector_output() {
        let vv = EnsightVectorVar::new("vel", vec![1.0, 2.0], vec![0.0, 0.0], vec![0.0, 0.0]);
        let mut buf = Vec::new();
        EnsightWriter::write_vector(&mut buf, &vv).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("vel"));
        assert!(s.contains("coordinates"));
    }

    // ── EnsightReader case ────────────────────────────────────────────────

    #[test]
    fn test_reader_case_minimal() {
        let input = "FORMAT\ntype: ensight gold\n\nGEOMETRY\nmodel: geom.geo\n";
        let reader = BufReader::new(input.as_bytes());
        let case = EnsightReader::read_case(reader).unwrap();
        assert_eq!(case.geometry_file, "geom.geo");
    }

    #[test]
    fn test_reader_case_with_scalar() {
        let input = "FORMAT\ntype: ensight gold\n\nGEOMETRY\nmodel: g.geo\n\nVARIABLE\nscalar per node: pressure p.escl\n";
        let reader = BufReader::new(input.as_bytes());
        let case = EnsightReader::read_case(reader).unwrap();
        assert_eq!(case.scalar_files.len(), 1);
        assert_eq!(case.scalar_files[0].0, "pressure");
    }

    #[test]
    fn test_reader_case_with_vector() {
        let input = "FORMAT\ntype: ensight gold\n\nGEOMETRY\nmodel: g.geo\n\nVARIABLE\nvector per node: velocity v.evec\n";
        let reader = BufReader::new(input.as_bytes());
        let case = EnsightReader::read_case(reader).unwrap();
        assert_eq!(case.vector_files.len(), 1);
        assert_eq!(case.vector_files[0].0, "velocity");
    }

    #[test]
    fn test_reader_case_roundtrip() {
        let mut case_written = EnsightCase::new("out.geo");
        case_written.add_scalar("p", "out_p.escl");
        case_written.add_vector("v", "out_v.evec");
        let mut buf = Vec::new();
        case_written.write_to(&mut buf).unwrap();
        let reader = BufReader::new(buf.as_slice());
        let case_read = EnsightReader::read_case(reader).unwrap();
        assert_eq!(case_read.geometry_file, "out.geo");
        assert_eq!(case_read.scalar_files.len(), 1);
        assert_eq!(case_read.vector_files.len(), 1);
    }

    // ── EnsightTimeSeries ─────────────────────────────────────────────────

    #[test]
    fn test_time_series_new_empty() {
        let ts = EnsightTimeSeries::new();
        assert_eq!(ts.n_steps(), 0);
    }

    #[test]
    fn test_time_series_push() {
        let mut ts = EnsightTimeSeries::new();
        ts.push(0.0, "geo_0.geo");
        ts.push(1.0, "geo_1.geo");
        assert_eq!(ts.n_steps(), 2);
        assert!((ts.times[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_time_series_write_transient() {
        let mut ts = EnsightTimeSeries::new();
        ts.push(0.0, "g0.geo");
        ts.push(0.5, "g1.geo");
        let mut buf = Vec::new();
        ts.write_transient_case(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("number of steps: 2"));
        assert!(s.contains("TIME"));
    }

    #[test]
    fn test_time_series_default() {
        let ts = EnsightTimeSeries::default();
        assert_eq!(ts.n_steps(), 0);
    }

    // ── write_ensight_case integration ────────────────────────────────────

    #[test]
    fn test_write_ensight_case_creates_files() {
        let tmp_dir = std::env::temp_dir();
        let base = tmp_dir.join("test_ensight_case_output");
        let path = base.to_str().unwrap();

        let mut geo = EnsightGeometry::new();
        geo.add_node(0.0, 0.0, 0.0);
        geo.add_node(1.0, 0.0, 0.0);
        geo.add_node(0.0, 1.0, 0.0);
        let mut part = EnsightPart::new(1, "surf", "tria3");
        part.connectivity = vec![0, 1, 2];
        geo.add_part(part);

        let scalars = vec![EnsightScalarVar::new("pressure", vec![1.0, 2.0, 3.0])];
        let vectors = vec![EnsightVectorVar::new(
            "velocity",
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        )];

        write_ensight_case(path, &geo, &scalars, &vectors).unwrap();

        assert!(tmp_dir.join("test_ensight_case_output.case").exists());
        assert!(tmp_dir.join("test_ensight_case_output.geo").exists());
    }

    #[test]
    fn test_write_ensight_case_no_vars() {
        let tmp_dir = std::env::temp_dir();
        let base = tmp_dir.join("test_ensight_novar");
        let path = base.to_str().unwrap();

        let mut geo = EnsightGeometry::new();
        geo.add_node(0.0, 0.0, 0.0);
        let mut part = EnsightPart::new(1, "pt", "point");
        part.connectivity = vec![0];
        geo.add_part(part);

        write_ensight_case(path, &geo, &[], &[]).unwrap();
        assert!(tmp_dir.join("test_ensight_novar.case").exists());
        assert!(tmp_dir.join("test_ensight_novar.geo").exists());
    }

    // ── Reader geo roundtrip ──────────────────────────────────────────────

    #[test]
    fn test_reader_geo_roundtrip_nodes() {
        let mut geo = EnsightGeometry::new();
        geo.add_node(1.0, 2.0, 3.0);
        geo.add_node(4.0, 5.0, 6.0);
        let mut part = EnsightPart::new(1, "pts", "point");
        part.connectivity = vec![0, 1];
        geo.add_part(part);

        let mut buf = Vec::new();
        EnsightWriter::write_geo(&mut buf, &geo).unwrap();

        let reader = BufReader::new(buf.as_slice());
        let geo2 = EnsightReader::read_geo(reader).unwrap();
        assert_eq!(geo2.n_nodes(), 2);
        assert!((geo2.x[0] - 1.0).abs() < 1e-4, "x0={}", geo2.x[0]);
        assert!((geo2.y[0] - 2.0).abs() < 1e-4, "y0={}", geo2.y[0]);
        assert!((geo2.z[0] - 3.0).abs() < 1e-4, "z0={}", geo2.z[0]);
    }
}
