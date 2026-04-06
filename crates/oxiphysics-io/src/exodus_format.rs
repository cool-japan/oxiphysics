#![allow(clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exodus II mesh format reader and writer (Sierra/SEACAS ASCII subset).
//!
//! Implements a readable ASCII representation of the Exodus II format,
//! including element blocks, node sets, side sets and scalar field variables.

use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, Write as IoWrite};

use crate::Error;
use crate::Result;

// ── Element types ─────────────────────────────────────────────────────────────

/// Element types supported by the Exodus II format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExodusElementType {
    /// 3-node linear triangle.
    Tri3,
    /// 6-node quadratic triangle.
    Tri6,
    /// 4-node bilinear quadrilateral.
    Quad4,
    /// 8-node serendipity quadrilateral.
    Quad8,
    /// 4-node linear tetrahedron.
    Tet4,
    /// 10-node quadratic tetrahedron.
    Tet10,
    /// 8-node linear hexahedron.
    Hex8,
    /// 20-node serendipity hexahedron.
    Hex20,
    /// 2-node bar (beam / truss).
    Bar2,
}

impl ExodusElementType {
    /// Return the canonical Exodus II string name for this element type.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Tri3 => "TRI3",
            Self::Tri6 => "TRI6",
            Self::Quad4 => "QUAD4",
            Self::Quad8 => "QUAD8",
            Self::Tet4 => "TET4",
            Self::Tet10 => "TET10",
            Self::Hex8 => "HEX8",
            Self::Hex20 => "HEX20",
            Self::Bar2 => "BAR2",
        }
    }

    /// Number of nodes per element for this type.
    pub fn nodes_per_element(self) -> usize {
        match self {
            Self::Tri3 => 3,
            Self::Tri6 => 6,
            Self::Quad4 => 4,
            Self::Quad8 => 8,
            Self::Tet4 => 4,
            Self::Tet10 => 10,
            Self::Hex8 => 8,
            Self::Hex20 => 20,
            Self::Bar2 => 2,
        }
    }

    /// Parse from the canonical Exodus II string name (case-insensitive).
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_ascii_uppercase().as_str() {
            "TRI3" => Some(Self::Tri3),
            "TRI6" => Some(Self::Tri6),
            "QUAD4" => Some(Self::Quad4),
            "QUAD8" => Some(Self::Quad8),
            "TET4" => Some(Self::Tet4),
            "TET10" => Some(Self::Tet10),
            "HEX8" => Some(Self::Hex8),
            "HEX20" => Some(Self::Hex20),
            "BAR2" => Some(Self::Bar2),
            _ => None,
        }
    }
}

// ── Node set ─────────────────────────────────────────────────────────────────

/// A named set of nodes (boundary condition set).
#[derive(Debug, Clone)]
pub struct ExodusNodeSet {
    /// Unique identifier for this node set.
    pub id: usize,
    /// Human-readable name.
    pub name: String,
    /// Node IDs (0-based) in this set.
    pub node_ids: Vec<usize>,
}

impl ExodusNodeSet {
    /// Create a new node set.
    pub fn new(id: usize, name: impl Into<String>, node_ids: Vec<usize>) -> Self {
        Self {
            id,
            name: name.into(),
            node_ids,
        }
    }
}

// ── Side set ─────────────────────────────────────────────────────────────────

/// A named set of element sides (face boundary conditions).
#[derive(Debug, Clone)]
pub struct ExodusSideSet {
    /// Unique identifier for this side set.
    pub id: usize,
    /// Human-readable name.
    pub name: String,
    /// Element IDs (0-based) for each entry.
    pub elem_ids: Vec<usize>,
    /// Local side IDs within each element.
    pub side_ids: Vec<usize>,
}

impl ExodusSideSet {
    /// Create a new side set.
    pub fn new(
        id: usize,
        name: impl Into<String>,
        elem_ids: Vec<usize>,
        side_ids: Vec<usize>,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            elem_ids,
            side_ids,
        }
    }
}

// ── Element block ─────────────────────────────────────────────────────────────

/// An element block (a group of elements sharing the same type).
#[derive(Debug, Clone)]
pub struct ExodusBlock {
    /// Block identifier.
    pub id: usize,
    /// Element type for all elements in this block.
    pub element_type: ExodusElementType,
    /// Connectivity: one row per element, listing node IDs (0-based).
    pub connectivity: Vec<Vec<usize>>,
}

impl ExodusBlock {
    /// Create a new element block.
    pub fn new(id: usize, element_type: ExodusElementType, connectivity: Vec<Vec<usize>>) -> Self {
        Self {
            id,
            element_type,
            connectivity,
        }
    }

    /// Number of elements in this block.
    pub fn num_elements(&self) -> usize {
        self.connectivity.len()
    }
}

// ── Exodus mesh ───────────────────────────────────────────────────────────────

/// A complete Exodus II mesh.
#[derive(Debug, Clone)]
pub struct ExodusMesh {
    /// Node coordinates: one `[x, y, z]` per node.
    pub nodes: Vec<[f64; 3]>,
    /// Element blocks.
    pub blocks: Vec<ExodusBlock>,
    /// Named node sets.
    pub node_sets: Vec<ExodusNodeSet>,
    /// Named side sets.
    pub side_sets: Vec<ExodusSideSet>,
    /// Named scalar field variables (one value per node).
    pub variables: HashMap<String, Vec<f64>>,
}

impl ExodusMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            blocks: Vec::new(),
            node_sets: Vec::new(),
            side_sets: Vec::new(),
            variables: HashMap::new(),
        }
    }

    /// Total number of elements across all blocks.
    pub fn total_elements(&self) -> usize {
        self.blocks.iter().map(|b| b.num_elements()).sum()
    }
}

impl Default for ExodusMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ── Writer ────────────────────────────────────────────────────────────────────

/// Writes an `ExodusMesh` to a readable ASCII-like text file.
pub struct ExodusWriter;

impl ExodusWriter {
    /// Create a new `ExodusWriter`.
    pub fn new() -> Self {
        Self
    }

    /// Write `mesh` to `path` in ASCII Exodus-like format.
    ///
    /// The file uses section headers (`NODES`, `BLOCK`, `NODE_SET`,
    /// `SIDE_SET`, `VARIABLE`) that `ExodusReader::parse` can read back.
    pub fn write_text(&self, mesh: &ExodusMesh, path: &str) -> Result<()> {
        let mut f = fs::File::create(path).map_err(Error::Io)?;

        // Header
        writeln!(f, "# Exodus II ASCII (OxiPhysics)")?;
        writeln!(f, "NUM_NODES {}", mesh.nodes.len())?;
        writeln!(f, "NUM_BLOCKS {}", mesh.blocks.len())?;
        writeln!(f, "NUM_NODE_SETS {}", mesh.node_sets.len())?;
        writeln!(f, "NUM_SIDE_SETS {}", mesh.side_sets.len())?;
        writeln!(f, "NUM_VARIABLES {}", mesh.variables.len())?;
        writeln!(f)?;

        // Nodes
        writeln!(f, "NODES")?;
        for (idx, n) in mesh.nodes.iter().enumerate() {
            writeln!(f, "{} {} {} {}", idx, n[0], n[1], n[2])?;
        }
        writeln!(f)?;

        // Blocks
        for blk in &mesh.blocks {
            writeln!(
                f,
                "BLOCK {} {} {}",
                blk.id,
                blk.element_type.as_str(),
                blk.connectivity.len()
            )?;
            for row in &blk.connectivity {
                let ids: Vec<String> = row.iter().map(|v| v.to_string()).collect();
                writeln!(f, "{}", ids.join(" "))?;
            }
            writeln!(f)?;
        }

        // Node sets
        for ns in &mesh.node_sets {
            writeln!(f, "NODE_SET {} {} {}", ns.id, ns.name, ns.node_ids.len())?;
            let ids: Vec<String> = ns.node_ids.iter().map(|v| v.to_string()).collect();
            writeln!(f, "{}", ids.join(" "))?;
            writeln!(f)?;
        }

        // Side sets
        for ss in &mesh.side_sets {
            writeln!(f, "SIDE_SET {} {} {}", ss.id, ss.name, ss.elem_ids.len())?;
            for (e, s) in ss.elem_ids.iter().zip(ss.side_ids.iter()) {
                writeln!(f, "{} {}", e, s)?;
            }
            writeln!(f)?;
        }

        // Variables
        let mut var_names: Vec<&String> = mesh.variables.keys().collect();
        var_names.sort();
        for name in var_names {
            let vals = &mesh.variables[name];
            writeln!(f, "VARIABLE {} {}", name, vals.len())?;
            for v in vals {
                writeln!(f, "{}", v)?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

impl Default for ExodusWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ── Reader ────────────────────────────────────────────────────────────────────

/// Reads an `ExodusMesh` from an ASCII file written by `ExodusWriter`.
pub struct ExodusReader;

impl ExodusReader {
    /// Create a new `ExodusReader`.
    pub fn new() -> Self {
        Self
    }

    /// Parse an Exodus ASCII file at `path` and return an `ExodusMesh`.
    pub fn parse(&self, path: &str) -> Result<ExodusMesh> {
        let file = fs::File::open(path).map_err(Error::Io)?;
        let reader = BufReader::new(file);
        let mut lines: Vec<String> = Vec::new();
        for line in reader.lines() {
            let l = line.map_err(Error::Io)?;
            let trimmed = l.trim().to_string();
            if !trimmed.starts_with('#') && !trimmed.is_empty() {
                lines.push(trimmed);
            }
        }

        let mut mesh = ExodusMesh::new();
        let mut pos = 0usize;

        // Skip header counts (NUM_* lines)
        while pos < lines.len() {
            let l = &lines[pos];
            if l.starts_with("NUM_") {
                pos += 1;
            } else {
                break;
            }
        }

        while pos < lines.len() {
            let l = lines[pos].clone();
            let parts: Vec<&str> = l.split_whitespace().collect();
            if parts.is_empty() {
                pos += 1;
                continue;
            }
            match parts[0] {
                "NODES" => {
                    pos += 1;
                    while pos < lines.len() {
                        let row: Vec<&str> = lines[pos].split_whitespace().collect();
                        if row.len() < 4 {
                            break;
                        }
                        // row[0] = index (integer), row[1..4] = xyz (floats)
                        // If row[0] cannot be parsed as usize this is not a node line.
                        if row[0].parse::<usize>().is_err() {
                            break;
                        }
                        let x = row[1]
                            .parse::<f64>()
                            .map_err(|e| Error::Parse(e.to_string()))?;
                        let y = row[2]
                            .parse::<f64>()
                            .map_err(|e| Error::Parse(e.to_string()))?;
                        let z = row[3]
                            .parse::<f64>()
                            .map_err(|e| Error::Parse(e.to_string()))?;
                        mesh.nodes.push([x, y, z]);
                        pos += 1;
                    }
                }
                "BLOCK" => {
                    if parts.len() < 4 {
                        return Err(Error::Parse("malformed BLOCK header".into()));
                    }
                    let id: usize = parts[1]
                        .parse()
                        .map_err(|e: std::num::ParseIntError| Error::Parse(e.to_string()))?;
                    let etype = ExodusElementType::from_str(parts[2]).ok_or_else(|| {
                        Error::Parse(format!("unknown element type: {}", parts[2]))
                    })?;
                    let count: usize = parts[3]
                        .parse()
                        .map_err(|e: std::num::ParseIntError| Error::Parse(e.to_string()))?;
                    pos += 1;
                    let mut connectivity = Vec::with_capacity(count);
                    for _ in 0..count {
                        if pos >= lines.len() {
                            return Err(Error::Parse("unexpected end in BLOCK".into()));
                        }
                        let row: Vec<usize> = lines[pos]
                            .split_whitespace()
                            .map(|s| s.parse::<usize>().map_err(|e| Error::Parse(e.to_string())))
                            .collect::<Result<Vec<_>>>()?;
                        connectivity.push(row);
                        pos += 1;
                    }
                    mesh.blocks.push(ExodusBlock::new(id, etype, connectivity));
                }
                "NODE_SET" => {
                    if parts.len() < 4 {
                        return Err(Error::Parse("malformed NODE_SET header".into()));
                    }
                    let id: usize = parts[1]
                        .parse()
                        .map_err(|e: std::num::ParseIntError| Error::Parse(e.to_string()))?;
                    let name = parts[2].to_string();
                    let count: usize = parts[3]
                        .parse()
                        .map_err(|e: std::num::ParseIntError| Error::Parse(e.to_string()))?;
                    pos += 1;
                    let mut node_ids = Vec::with_capacity(count);
                    if pos < lines.len() && count > 0 {
                        let row: Vec<usize> = lines[pos]
                            .split_whitespace()
                            .map(|s| s.parse::<usize>().map_err(|e| Error::Parse(e.to_string())))
                            .collect::<Result<Vec<_>>>()?;
                        node_ids = row;
                        pos += 1;
                    } else if count == 0 {
                        // empty set — nothing to read
                    }
                    mesh.node_sets.push(ExodusNodeSet::new(id, name, node_ids));
                }
                "SIDE_SET" => {
                    if parts.len() < 4 {
                        return Err(Error::Parse("malformed SIDE_SET header".into()));
                    }
                    let id: usize = parts[1]
                        .parse()
                        .map_err(|e: std::num::ParseIntError| Error::Parse(e.to_string()))?;
                    let name = parts[2].to_string();
                    let count: usize = parts[3]
                        .parse()
                        .map_err(|e: std::num::ParseIntError| Error::Parse(e.to_string()))?;
                    pos += 1;
                    let mut elem_ids = Vec::with_capacity(count);
                    let mut side_ids = Vec::with_capacity(count);
                    for _ in 0..count {
                        if pos >= lines.len() {
                            return Err(Error::Parse("unexpected end in SIDE_SET".into()));
                        }
                        let row: Vec<&str> = lines[pos].split_whitespace().collect();
                        if row.len() < 2 {
                            return Err(Error::Parse("malformed SIDE_SET entry".into()));
                        }
                        elem_ids.push(
                            row[0]
                                .parse::<usize>()
                                .map_err(|e| Error::Parse(e.to_string()))?,
                        );
                        side_ids.push(
                            row[1]
                                .parse::<usize>()
                                .map_err(|e| Error::Parse(e.to_string()))?,
                        );
                        pos += 1;
                    }
                    mesh.side_sets
                        .push(ExodusSideSet::new(id, name, elem_ids, side_ids));
                }
                "VARIABLE" => {
                    if parts.len() < 3 {
                        return Err(Error::Parse("malformed VARIABLE header".into()));
                    }
                    let name = parts[1].to_string();
                    let count: usize = parts[2]
                        .parse()
                        .map_err(|e: std::num::ParseIntError| Error::Parse(e.to_string()))?;
                    pos += 1;
                    let mut vals = Vec::with_capacity(count);
                    for _ in 0..count {
                        if pos >= lines.len() {
                            return Err(Error::Parse("unexpected end in VARIABLE".into()));
                        }
                        let v: f64 = lines[pos]
                            .trim()
                            .parse()
                            .map_err(|e: std::num::ParseFloatError| Error::Parse(e.to_string()))?;
                        vals.push(v);
                        pos += 1;
                    }
                    mesh.variables.insert(name, vals);
                }
                _ => {
                    pos += 1;
                }
            }
        }

        Ok(mesh)
    }
}

impl Default for ExodusReader {
    fn default() -> Self {
        Self::new()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_tri_mesh() -> ExodusMesh {
        let mut m = ExodusMesh::new();
        m.nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.5, 1.0, 0.0],
        ];
        m.blocks.push(ExodusBlock::new(
            1,
            ExodusElementType::Tri3,
            vec![vec![0, 1, 2], vec![1, 3, 2]],
        ));
        m
    }

    fn tmp_path(name: &str) -> String {
        format!("/tmp/oxiphysics_exodus_test_{name}.exo")
    }

    // ── ExodusElementType ────────────────────────────────────────────────────

    #[test]
    fn element_type_as_str() {
        assert_eq!(ExodusElementType::Tri3.as_str(), "TRI3");
        assert_eq!(ExodusElementType::Hex20.as_str(), "HEX20");
        assert_eq!(ExodusElementType::Bar2.as_str(), "BAR2");
    }

    #[test]
    fn element_type_from_str_roundtrip() {
        for et in [
            ExodusElementType::Tri3,
            ExodusElementType::Tri6,
            ExodusElementType::Quad4,
            ExodusElementType::Quad8,
            ExodusElementType::Tet4,
            ExodusElementType::Tet10,
            ExodusElementType::Hex8,
            ExodusElementType::Hex20,
            ExodusElementType::Bar2,
        ] {
            let parsed = ExodusElementType::from_str(et.as_str());
            assert_eq!(parsed, Some(et), "round-trip failed for {:?}", et);
        }
    }

    #[test]
    fn element_type_from_str_unknown() {
        assert!(ExodusElementType::from_str("UNKNOWN").is_none());
        assert!(ExodusElementType::from_str("").is_none());
    }

    #[test]
    fn element_type_from_str_case_insensitive() {
        assert_eq!(
            ExodusElementType::from_str("tri3"),
            Some(ExodusElementType::Tri3)
        );
        assert_eq!(
            ExodusElementType::from_str("Hex8"),
            Some(ExodusElementType::Hex8)
        );
    }

    #[test]
    fn element_type_nodes_per_element() {
        assert_eq!(ExodusElementType::Tri3.nodes_per_element(), 3);
        assert_eq!(ExodusElementType::Hex20.nodes_per_element(), 20);
        assert_eq!(ExodusElementType::Bar2.nodes_per_element(), 2);
        assert_eq!(ExodusElementType::Tet10.nodes_per_element(), 10);
    }

    // ── ExodusBlock ──────────────────────────────────────────────────────────

    #[test]
    fn block_num_elements() {
        let blk = ExodusBlock::new(
            1,
            ExodusElementType::Quad4,
            vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]],
        );
        assert_eq!(blk.num_elements(), 2);
    }

    #[test]
    fn block_empty_connectivity() {
        let blk = ExodusBlock::new(2, ExodusElementType::Tet4, vec![]);
        assert_eq!(blk.num_elements(), 0);
    }

    #[test]
    fn block_stores_element_type() {
        let blk = ExodusBlock::new(
            3,
            ExodusElementType::Hex8,
            vec![vec![0, 1, 2, 3, 4, 5, 6, 7]],
        );
        assert_eq!(blk.element_type, ExodusElementType::Hex8);
    }

    // ── ExodusNodeSet ────────────────────────────────────────────────────────

    #[test]
    fn node_set_stores_data() {
        let ns = ExodusNodeSet::new(10, "inlet", vec![0, 1, 2]);
        assert_eq!(ns.id, 10);
        assert_eq!(ns.name, "inlet");
        assert_eq!(ns.node_ids, vec![0, 1, 2]);
    }

    #[test]
    fn node_set_empty() {
        let ns = ExodusNodeSet::new(1, "empty", vec![]);
        assert!(ns.node_ids.is_empty());
    }

    // ── ExodusSideSet ────────────────────────────────────────────────────────

    #[test]
    fn side_set_stores_data() {
        let ss = ExodusSideSet::new(5, "wall", vec![0, 1], vec![2, 3]);
        assert_eq!(ss.id, 5);
        assert_eq!(ss.elem_ids, vec![0, 1]);
        assert_eq!(ss.side_ids, vec![2, 3]);
    }

    #[test]
    fn side_set_mismatched_lengths_allowed() {
        // The struct stores them independently — no validation required
        let ss = ExodusSideSet::new(1, "ss", vec![0], vec![0, 1]);
        assert_eq!(ss.elem_ids.len(), 1);
        assert_eq!(ss.side_ids.len(), 2);
    }

    // ── ExodusMesh ───────────────────────────────────────────────────────────

    #[test]
    fn mesh_default_is_empty() {
        let m = ExodusMesh::default();
        assert!(m.nodes.is_empty());
        assert!(m.blocks.is_empty());
        assert!(m.node_sets.is_empty());
        assert!(m.side_sets.is_empty());
        assert!(m.variables.is_empty());
    }

    #[test]
    fn mesh_total_elements() {
        let mut m = ExodusMesh::new();
        m.blocks.push(ExodusBlock::new(
            1,
            ExodusElementType::Tri3,
            vec![vec![0, 1, 2]; 3],
        ));
        m.blocks.push(ExodusBlock::new(
            2,
            ExodusElementType::Quad4,
            vec![vec![0, 1, 2, 3]; 5],
        ));
        assert_eq!(m.total_elements(), 8);
    }

    #[test]
    fn mesh_variable_storage() {
        let mut m = ExodusMesh::new();
        m.variables
            .insert("temperature".into(), vec![1.0, 2.0, 3.0]);
        assert_eq!(m.variables["temperature"], vec![1.0, 2.0, 3.0]);
    }

    // ── Write / read roundtrip ───────────────────────────────────────────────

    #[test]
    fn write_read_nodes_roundtrip() {
        let m = simple_tri_mesh();
        let path = tmp_path("nodes_rtrip");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.nodes.len(), m.nodes.len());
        for (a, b) in m.nodes.iter().zip(m2.nodes.iter()) {
            for k in 0..3 {
                assert!((a[k] - b[k]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn write_read_blocks_roundtrip() {
        let m = simple_tri_mesh();
        let path = tmp_path("blocks_rtrip");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.blocks.len(), 1);
        assert_eq!(m2.blocks[0].element_type, ExodusElementType::Tri3);
        assert_eq!(m2.blocks[0].connectivity.len(), 2);
        assert_eq!(m2.blocks[0].connectivity[0], vec![0, 1, 2]);
    }

    #[test]
    fn write_read_node_set_roundtrip() {
        let mut m = simple_tri_mesh();
        m.node_sets.push(ExodusNodeSet::new(1, "inlet", vec![0, 1]));
        let path = tmp_path("nset_rtrip");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.node_sets.len(), 1);
        assert_eq!(m2.node_sets[0].name, "inlet");
        assert_eq!(m2.node_sets[0].node_ids, vec![0, 1]);
    }

    #[test]
    fn write_read_side_set_roundtrip() {
        let mut m = simple_tri_mesh();
        m.side_sets
            .push(ExodusSideSet::new(1, "wall", vec![0, 1], vec![2, 3]));
        let path = tmp_path("sset_rtrip");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.side_sets.len(), 1);
        assert_eq!(m2.side_sets[0].name, "wall");
        assert_eq!(m2.side_sets[0].elem_ids, vec![0, 1]);
        assert_eq!(m2.side_sets[0].side_ids, vec![2, 3]);
    }

    #[test]
    fn write_read_variable_roundtrip() {
        let mut m = simple_tri_mesh();
        m.variables
            .insert("pressure".into(), vec![1.0, 2.0, 3.0, 4.0]);
        let path = tmp_path("var_rtrip");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert!(m2.variables.contains_key("pressure"));
        let v = &m2.variables["pressure"];
        assert_eq!(v.len(), 4);
        assert!((v[0] - 1.0).abs() < 1e-12);
        assert!((v[3] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn write_read_multiple_blocks() {
        let mut m = ExodusMesh::new();
        m.nodes = vec![[0.0; 3]; 8];
        m.blocks.push(ExodusBlock::new(
            1,
            ExodusElementType::Tri3,
            vec![vec![0, 1, 2]],
        ));
        m.blocks.push(ExodusBlock::new(
            2,
            ExodusElementType::Quad4,
            vec![vec![0, 1, 2, 3]],
        ));
        m.blocks.push(ExodusBlock::new(
            3,
            ExodusElementType::Bar2,
            vec![vec![4, 5]],
        ));
        let path = tmp_path("multi_blk");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.blocks.len(), 3);
        assert_eq!(m2.blocks[1].element_type, ExodusElementType::Quad4);
        assert_eq!(m2.blocks[2].element_type, ExodusElementType::Bar2);
    }

    #[test]
    fn write_read_multiple_variables() {
        let mut m = simple_tri_mesh();
        let mut vars: HashMap<String, Vec<f64>> = HashMap::new();
        vars.insert("vel_x".into(), vec![1.0, 2.0, 3.0, 4.0]);
        vars.insert("vel_y".into(), vec![0.1, 0.2, 0.3, 0.4]);
        m.variables = vars;
        let path = tmp_path("multi_var");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.variables.len(), 2);
        assert!(m2.variables.contains_key("vel_x"));
        assert!(m2.variables.contains_key("vel_y"));
    }

    #[test]
    fn write_creates_file() {
        let m = simple_tri_mesh();
        let path = tmp_path("create_check");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        assert!(std::path::Path::new(&path).exists());
    }

    #[test]
    fn read_nonexistent_file_returns_error() {
        let result = ExodusReader::new().parse("/tmp/nonexistent_oxiphysics_exodus.exo");
        assert!(result.is_err());
    }

    #[test]
    fn write_read_tet_mesh() {
        let mut m = ExodusMesh::new();
        m.nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        m.blocks.push(ExodusBlock::new(
            1,
            ExodusElementType::Tet4,
            vec![vec![0, 1, 2, 3]],
        ));
        let path = tmp_path("tet_mesh");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.blocks[0].element_type, ExodusElementType::Tet4);
        assert_eq!(m2.blocks[0].connectivity[0], vec![0, 1, 2, 3]);
    }

    #[test]
    fn write_read_hex_mesh() {
        let mut m = ExodusMesh::new();
        m.nodes = vec![[0.0; 3]; 8];
        m.blocks.push(ExodusBlock::new(
            1,
            ExodusElementType::Hex8,
            vec![vec![0, 1, 2, 3, 4, 5, 6, 7]],
        ));
        let path = tmp_path("hex_mesh");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.blocks[0].element_type, ExodusElementType::Hex8);
        assert_eq!(m2.blocks[0].connectivity[0].len(), 8);
    }

    #[test]
    fn writer_default_works() {
        let w = ExodusWriter;
        let m = ExodusMesh::new();
        let path = tmp_path("writer_default");
        w.write_text(&m, &path).unwrap();
    }

    #[test]
    fn reader_default_works() {
        let _ = ExodusReader;
    }

    #[test]
    fn mesh_clone() {
        let m = simple_tri_mesh();
        let m2 = m.clone();
        assert_eq!(m2.nodes.len(), m.nodes.len());
        assert_eq!(m2.blocks.len(), m.blocks.len());
    }

    #[test]
    fn multiple_node_sets_roundtrip() {
        let mut m = simple_tri_mesh();
        m.node_sets.push(ExodusNodeSet::new(1, "inlet", vec![0]));
        m.node_sets
            .push(ExodusNodeSet::new(2, "outlet", vec![1, 2]));
        let path = tmp_path("multi_nset");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.node_sets.len(), 2);
    }

    #[test]
    fn multiple_side_sets_roundtrip() {
        let mut m = simple_tri_mesh();
        m.side_sets
            .push(ExodusSideSet::new(1, "top", vec![0], vec![1]));
        m.side_sets
            .push(ExodusSideSet::new(2, "bottom", vec![1], vec![0]));
        let path = tmp_path("multi_sset");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.side_sets.len(), 2);
    }

    #[test]
    fn empty_mesh_roundtrip() {
        let m = ExodusMesh::new();
        let path = tmp_path("empty_mesh");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert!(m2.nodes.is_empty());
        assert!(m2.blocks.is_empty());
    }

    #[test]
    fn block_id_preserved() {
        let mut m = ExodusMesh::new();
        m.nodes = vec![[0.0; 3]; 3];
        m.blocks.push(ExodusBlock::new(
            42,
            ExodusElementType::Tri3,
            vec![vec![0, 1, 2]],
        ));
        let path = tmp_path("block_id");
        ExodusWriter::new().write_text(&m, &path).unwrap();
        let m2 = ExodusReader::new().parse(&path).unwrap();
        assert_eq!(m2.blocks[0].id, 42);
    }
}
