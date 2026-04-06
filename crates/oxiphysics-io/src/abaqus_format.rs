#![allow(clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Abaqus INP file format I/O.
//!
//! Provides readers and writers for the Abaqus `.inp` finite-element input
//! format, covering node/element definitions, material properties, section
//! assignments, and boundary conditions.

use std::fmt::Write as FmtWrite;
use std::fs;
use std::io::{self, BufRead};

use crate::Error as IoError;

// ---------------------------------------------------------------------------
// AbaqusNode
// ---------------------------------------------------------------------------

/// A single node in an Abaqus mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct AbaqusNode {
    /// 1-based node identifier.
    pub id: usize,
    /// Node coordinates `[x, y, z]`.
    pub coordinates: [f64; 3],
}

impl AbaqusNode {
    /// Create a new node.
    pub fn new(id: usize, coordinates: [f64; 3]) -> Self {
        Self { id, coordinates }
    }
}

// ---------------------------------------------------------------------------
// ElementType
// ---------------------------------------------------------------------------

/// Abaqus element type identifier.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ElementType {
    /// 4-node 3-D tetrahedral element.
    C3D4,
    /// 8-node 3-D hexahedral element.
    C3D8,
    /// 4-node shell element.
    S4,
    /// 2-node truss element.
    T3D2,
    /// Unknown / unsupported element type.
    Unknown(String),
}

impl ElementType {
    /// Canonical Abaqus keyword string.
    pub fn as_str(&self) -> &str {
        match self {
            ElementType::C3D4 => "C3D4",
            ElementType::C3D8 => "C3D8",
            ElementType::S4 => "S4",
            ElementType::T3D2 => "T3D2",
            ElementType::Unknown(s) => s.as_str(),
        }
    }

    /// Parse from string (case-insensitive).
    pub fn from_str(s: &str) -> Self {
        match s.trim().to_uppercase().as_str() {
            "C3D4" => ElementType::C3D4,
            "C3D8" => ElementType::C3D8,
            "S4" => ElementType::S4,
            "T3D2" => ElementType::T3D2,
            other => ElementType::Unknown(other.to_string()),
        }
    }
}

// ---------------------------------------------------------------------------
// AbaqusElement
// ---------------------------------------------------------------------------

/// A single finite element in an Abaqus mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct AbaqusElement {
    /// 1-based element identifier.
    pub id: usize,
    /// Element topology type.
    pub element_type: ElementType,
    /// Node connectivity (1-based node ids).
    pub node_ids: Vec<usize>,
}

impl AbaqusElement {
    /// Create a new element.
    pub fn new(id: usize, element_type: ElementType, node_ids: Vec<usize>) -> Self {
        Self {
            id,
            element_type,
            node_ids,
        }
    }
}

// ---------------------------------------------------------------------------
// AbaqusSection
// ---------------------------------------------------------------------------

/// An Abaqus section assignment mapping elements to a material.
#[derive(Debug, Clone)]
pub struct AbaqusSection {
    /// Section name.
    pub name: String,
    /// Name of the associated material.
    pub material_name: String,
    /// Element set indices assigned to this section.
    pub elements: Vec<usize>,
}

impl AbaqusSection {
    /// Create a new section.
    pub fn new(
        name: impl Into<String>,
        material_name: impl Into<String>,
        elements: Vec<usize>,
    ) -> Self {
        Self {
            name: name.into(),
            material_name: material_name.into(),
            elements,
        }
    }
}

// ---------------------------------------------------------------------------
// AbaqusMaterial
// ---------------------------------------------------------------------------

/// Elastic material properties used by Abaqus.
#[derive(Debug, Clone)]
pub struct ElasticProps {
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
}

/// Simple isotropic hardening plasticity.
#[derive(Debug, Clone)]
pub struct PlasticProps {
    /// Initial yield stress (Pa).
    pub yield_stress: f64,
    /// Linear hardening modulus (Pa).
    pub hardening_modulus: f64,
}

/// Material definition for an Abaqus model.
#[derive(Debug, Clone)]
pub struct AbaqusMaterial {
    /// Material name.
    pub name: String,
    /// Elastic properties.
    pub elastic: ElasticProps,
    /// Material density (kg/m³).
    pub density: f64,
    /// Optional plastic properties.
    pub plastic: Option<PlasticProps>,
}

impl AbaqusMaterial {
    /// Create a purely elastic material.
    pub fn new_elastic(
        name: impl Into<String>,
        young_modulus: f64,
        poisson_ratio: f64,
        density: f64,
    ) -> Self {
        Self {
            name: name.into(),
            elastic: ElasticProps {
                young_modulus,
                poisson_ratio,
            },
            density,
            plastic: None,
        }
    }

    /// Create an elastoplastic material.
    pub fn new_plastic(
        name: impl Into<String>,
        young_modulus: f64,
        poisson_ratio: f64,
        density: f64,
        yield_stress: f64,
        hardening_modulus: f64,
    ) -> Self {
        Self {
            name: name.into(),
            elastic: ElasticProps {
                young_modulus,
                poisson_ratio,
            },
            density,
            plastic: Some(PlasticProps {
                yield_stress,
                hardening_modulus,
            }),
        }
    }
}

// ---------------------------------------------------------------------------
// BoundaryCondition
// ---------------------------------------------------------------------------

/// Abaqus boundary condition type.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryCondition {
    /// All 6 DOFs fixed (clamped).
    Encastre {
        /// Node-set name the BC applies to.
        node_set: String,
    },
    /// Translation DOFs fixed, rotations free.
    Pinned {
        /// Node-set name the BC applies to.
        node_set: String,
    },
    /// Symmetry plane constraint.
    SymmetryPlane {
        /// Node-set name the BC applies to.
        node_set: String,
        /// Axis normal to the symmetry plane: 1=X, 2=Y, 3=Z.
        axis: u8,
    },
}

impl BoundaryCondition {
    /// Return the Abaqus keyword string for this BC type.
    pub fn keyword(&self) -> &'static str {
        match self {
            BoundaryCondition::Encastre { .. } => "ENCASTRE",
            BoundaryCondition::Pinned { .. } => "PINNED",
            BoundaryCondition::SymmetryPlane { .. } => "SYMMETRY",
        }
    }
}

// ---------------------------------------------------------------------------
// AbaqusMesh
// ---------------------------------------------------------------------------

/// Complete Abaqus mesh model.
#[derive(Debug, Clone, Default)]
pub struct AbaqusMesh {
    /// All nodes in the model.
    pub nodes: Vec<AbaqusNode>,
    /// All elements in the model.
    pub elements: Vec<AbaqusElement>,
    /// Section assignments.
    pub sections: Vec<AbaqusSection>,
    /// Material definitions.
    pub materials: Vec<AbaqusMaterial>,
    /// Boundary conditions.
    pub boundary_conditions: Vec<BoundaryCondition>,
}

impl AbaqusMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }
}

// ---------------------------------------------------------------------------
// AbaqusWriter
// ---------------------------------------------------------------------------

/// Writes an `AbaqusMesh` to an Abaqus `.inp` text file.
#[derive(Debug, Clone, Default)]
pub struct AbaqusWriter;

impl AbaqusWriter {
    /// Create a new writer.
    pub fn new() -> Self {
        Self
    }

    /// Serialize `mesh` to `path` in Abaqus INP format.
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn write(&self, mesh: &AbaqusMesh, path: &str) -> Result<(), IoError> {
        let mut buf = String::new();

        let _ = writeln!(buf, "** Generated by OxiPhysics AbaqusWriter");
        let _ = writeln!(buf, "*Heading");
        let _ = writeln!(buf, "OxiPhysics model");

        // --- Nodes ---
        let _ = writeln!(buf, "*Node");
        for node in &mesh.nodes {
            writeln!(
                buf,
                "{}, {:.15e}, {:.15e}, {:.15e}",
                node.id, node.coordinates[0], node.coordinates[1], node.coordinates[2]
            )
            .expect("operation should succeed");
        }

        // --- Elements (grouped by type) ---
        // Collect unique types
        let mut types_seen: Vec<String> = Vec::new();
        for el in &mesh.elements {
            let t = el.element_type.as_str().to_string();
            if !types_seen.contains(&t) {
                types_seen.push(t);
            }
        }

        for etype in &types_seen {
            let _ = writeln!(buf, "*Element, type={etype}");
            for el in mesh
                .elements
                .iter()
                .filter(|e| e.element_type.as_str() == etype)
            {
                let ids: Vec<String> = el.node_ids.iter().map(|n| n.to_string()).collect();
                let _ = writeln!(buf, "{}, {}", el.id, ids.join(", "));
            }
        }

        // --- Materials ---
        for mat in &mesh.materials {
            let _ = writeln!(buf, "*Material, name={}", mat.name);
            let _ = writeln!(buf, "*Density");
            let _ = writeln!(buf, "{:.15e}", mat.density);
            let _ = writeln!(buf, "*Elastic");
            writeln!(
                buf,
                "{:.15e}, {:.15e}",
                mat.elastic.young_modulus, mat.elastic.poisson_ratio
            )
            .expect("operation should succeed");
            if let Some(p) = &mat.plastic {
                let _ = writeln!(buf, "*Plastic");
                let _ = writeln!(buf, "{:.15e}, {:.15e}", p.yield_stress, p.hardening_modulus);
            }
        }

        // --- Sections ---
        for sec in &mesh.sections {
            let el_set: Vec<String> = sec.elements.iter().map(|e| e.to_string()).collect();
            writeln!(
                buf,
                "*Solid Section, elset={}, material={}",
                sec.name, sec.material_name
            )
            .expect("operation should succeed");
            if !el_set.is_empty() {
                let _ = writeln!(buf, "{}", el_set.join(", "));
            }
        }

        // --- Boundary Conditions ---
        for bc in &mesh.boundary_conditions {
            match bc {
                BoundaryCondition::Encastre { node_set } => {
                    let _ = writeln!(buf, "*Boundary");
                    let _ = writeln!(buf, "{node_set}, ENCASTRE");
                }
                BoundaryCondition::Pinned { node_set } => {
                    let _ = writeln!(buf, "*Boundary");
                    let _ = writeln!(buf, "{node_set}, PINNED");
                }
                BoundaryCondition::SymmetryPlane { node_set, axis } => {
                    let _ = writeln!(buf, "*Boundary");
                    let sym = match axis {
                        1 => "XSYMM",
                        2 => "YSYMM",
                        3 => "ZSYMM",
                        _ => "XSYMM",
                    };
                    let _ = writeln!(buf, "{node_set}, {sym}");
                }
            }
        }

        let _ = writeln!(buf, "*End Part");

        fs::write(path, buf).map_err(IoError::Io)
    }
}

// ---------------------------------------------------------------------------
// AbaqusReader
// ---------------------------------------------------------------------------

/// Parses a subset of Abaqus `.inp` files: `*NODE` and `*ELEMENT` sections.
#[derive(Debug, Clone, Default)]
pub struct AbaqusReader;

impl AbaqusReader {
    /// Create a new reader.
    pub fn new() -> Self {
        Self
    }

    /// Parse `path` and return an `AbaqusMesh` with nodes and elements
    /// populated.
    ///
    /// Only `*Node` and `*Element` keyword blocks are processed; all other
    /// blocks are silently skipped.
    pub fn parse(&self, path: &str) -> Result<AbaqusMesh, IoError> {
        let file = fs::File::open(path).map_err(IoError::Io)?;
        let reader = io::BufReader::new(file);

        let mut mesh = AbaqusMesh::new();
        let mut current_block = Block::None;

        for line_res in reader.lines() {
            let line = line_res.map_err(IoError::Io)?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with("**") {
                continue;
            }

            if trimmed.starts_with('*') {
                // New keyword
                let upper = trimmed.to_uppercase();
                if upper.starts_with("*NODE") && !upper.starts_with("*NSET") {
                    current_block = Block::Node;
                } else if upper.starts_with("*ELEMENT") {
                    // Extract type= parameter
                    let etype = Self::extract_param(trimmed, "TYPE")
                        .unwrap_or_else(|| "UNKNOWN".to_string());
                    current_block = Block::Element(ElementType::from_str(&etype));
                } else {
                    current_block = Block::None;
                }
                continue;
            }

            match &current_block {
                Block::Node => {
                    if let Some(node) = Self::parse_node_line(trimmed) {
                        mesh.nodes.push(node);
                    }
                }
                Block::Element(etype) => {
                    if let Some(el) = Self::parse_element_line(trimmed, etype.clone()) {
                        mesh.elements.push(el);
                    }
                }
                Block::None => {}
            }
        }

        Ok(mesh)
    }

    /// Extract value of `key=value` from a keyword line (case-insensitive).
    fn extract_param(line: &str, key: &str) -> Option<String> {
        let upper = line.to_uppercase();
        let key_eq = format!("{key}=");
        let pos = upper.find(&key_eq)?;
        let rest = &line[pos + key_eq.len()..];
        // Value ends at next comma or end of line
        let end = rest.find(',').unwrap_or(rest.len());
        Some(rest[..end].trim().to_string())
    }

    /// Parse a node data line: `id, x, y, z`
    fn parse_node_line(line: &str) -> Option<AbaqusNode> {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            return None;
        }
        let id: usize = parts[0].trim().parse().ok()?;
        let x: f64 = parts[1].trim().parse().ok()?;
        let y: f64 = parts[2].trim().parse().ok()?;
        let z: f64 = parts[3].trim().parse().ok()?;
        Some(AbaqusNode::new(id, [x, y, z]))
    }

    /// Parse an element data line: `id, n1, n2, ...`
    fn parse_element_line(line: &str, etype: ElementType) -> Option<AbaqusElement> {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 2 {
            return None;
        }
        let id: usize = parts[0].trim().parse().ok()?;
        let node_ids: Vec<usize> = parts[1..]
            .iter()
            .filter_map(|s| s.trim().parse().ok())
            .collect();
        if node_ids.is_empty() {
            return None;
        }
        Some(AbaqusElement::new(id, etype, node_ids))
    }
}

/// Internal parsing state machine.
#[derive(Debug, Clone)]
enum Block {
    /// Outside any recognized block.
    None,
    /// Inside a `*Node` block.
    Node,
    /// Inside an `*Element` block with the given type.
    Element(ElementType),
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- ElementType ---

    #[test]
    fn test_element_type_from_str_c3d4() {
        assert_eq!(ElementType::from_str("C3D4"), ElementType::C3D4);
    }

    #[test]
    fn test_element_type_from_str_case_insensitive() {
        assert_eq!(ElementType::from_str("c3d8"), ElementType::C3D8);
    }

    #[test]
    fn test_element_type_from_str_s4() {
        assert_eq!(ElementType::from_str("S4"), ElementType::S4);
    }

    #[test]
    fn test_element_type_from_str_t3d2() {
        assert_eq!(ElementType::from_str("T3D2"), ElementType::T3D2);
    }

    #[test]
    fn test_element_type_from_str_unknown() {
        match ElementType::from_str("FOOBAR") {
            ElementType::Unknown(s) => assert_eq!(s, "FOOBAR"),
            _ => panic!("expected Unknown"),
        }
    }

    #[test]
    fn test_element_type_as_str() {
        assert_eq!(ElementType::C3D4.as_str(), "C3D4");
        assert_eq!(ElementType::C3D8.as_str(), "C3D8");
        assert_eq!(ElementType::S4.as_str(), "S4");
        assert_eq!(ElementType::T3D2.as_str(), "T3D2");
    }

    // --- AbaqusNode ---

    #[test]
    fn test_node_new() {
        let n = AbaqusNode::new(1, [1.0, 2.0, 3.0]);
        assert_eq!(n.id, 1);
        assert_eq!(n.coordinates, [1.0, 2.0, 3.0]);
    }

    // --- AbaqusMaterial ---

    #[test]
    fn test_material_elastic() {
        let m = AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0);
        assert_eq!(m.name, "Steel");
        assert!((m.elastic.young_modulus - 210e9).abs() < 1.0);
        assert!(m.plastic.is_none());
    }

    #[test]
    fn test_material_plastic() {
        let m = AbaqusMaterial::new_plastic("Steel", 210e9, 0.3, 7800.0, 250e6, 1e9);
        assert!(m.plastic.is_some());
        let p = m.plastic.unwrap();
        assert!((p.yield_stress - 250e6).abs() < 1.0);
    }

    // --- BoundaryCondition ---

    #[test]
    fn test_bc_keyword_encastre() {
        let bc = BoundaryCondition::Encastre {
            node_set: "FIXED".to_string(),
        };
        assert_eq!(bc.keyword(), "ENCASTRE");
    }

    #[test]
    fn test_bc_keyword_pinned() {
        let bc = BoundaryCondition::Pinned {
            node_set: "PIN".to_string(),
        };
        assert_eq!(bc.keyword(), "PINNED");
    }

    #[test]
    fn test_bc_keyword_symmetry() {
        let bc = BoundaryCondition::SymmetryPlane {
            node_set: "SYM".to_string(),
            axis: 2,
        };
        assert_eq!(bc.keyword(), "SYMMETRY");
    }

    // --- AbaqusWriter / AbaqusReader roundtrip ---

    fn sample_mesh() -> AbaqusMesh {
        let mut mesh = AbaqusMesh::new();
        mesh.nodes = vec![
            AbaqusNode::new(1, [0.0, 0.0, 0.0]),
            AbaqusNode::new(2, [1.0, 0.0, 0.0]),
            AbaqusNode::new(3, [0.0, 1.0, 0.0]),
            AbaqusNode::new(4, [0.0, 0.0, 1.0]),
        ];
        mesh.elements = vec![AbaqusElement::new(1, ElementType::C3D4, vec![1, 2, 3, 4])];
        mesh
    }

    #[test]
    fn test_write_creates_file() {
        let path = "/tmp/oxiphysics_abaqus_test_write.inp";
        let mesh = sample_mesh();
        let writer = AbaqusWriter::new();
        writer.write(&mesh, path).expect("write failed");
        assert!(std::path::Path::new(path).exists());
    }

    #[test]
    fn test_roundtrip_node_count() {
        let path = "/tmp/oxiphysics_abaqus_roundtrip.inp";
        let mesh = sample_mesh();
        let writer = AbaqusWriter::new();
        writer.write(&mesh, path).expect("write failed");

        let reader = AbaqusReader::new();
        let parsed = reader.parse(path).expect("parse failed");
        assert_eq!(parsed.nodes.len(), 4);
    }

    #[test]
    fn test_roundtrip_element_count() {
        let path = "/tmp/oxiphysics_abaqus_rt_elem.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        assert_eq!(parsed.elements.len(), 1);
    }

    #[test]
    fn test_roundtrip_node_ids() {
        let path = "/tmp/oxiphysics_abaqus_rt_nodeids.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        let ids: Vec<usize> = parsed.nodes.iter().map(|n| n.id).collect();
        assert_eq!(ids, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_roundtrip_node_coordinates() {
        let path = "/tmp/oxiphysics_abaqus_rt_coords.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        let n1 = &parsed.nodes[0];
        assert!((n1.coordinates[0]).abs() < 1e-10);
        let n2 = &parsed.nodes[1];
        assert!((n2.coordinates[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_roundtrip_element_type() {
        let path = "/tmp/oxiphysics_abaqus_rt_etype.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        assert_eq!(parsed.elements[0].element_type, ElementType::C3D4);
    }

    #[test]
    fn test_roundtrip_element_nodes() {
        let path = "/tmp/oxiphysics_abaqus_rt_enodes.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        assert_eq!(parsed.elements[0].node_ids, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_roundtrip_element_id() {
        let path = "/tmp/oxiphysics_abaqus_rt_eid.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        assert_eq!(parsed.elements[0].id, 1);
    }

    #[test]
    fn test_multiple_element_types() {
        let path = "/tmp/oxiphysics_abaqus_multtype.inp";
        let mut mesh = AbaqusMesh::new();
        mesh.nodes = vec![
            AbaqusNode::new(1, [0.0, 0.0, 0.0]),
            AbaqusNode::new(2, [1.0, 0.0, 0.0]),
            AbaqusNode::new(3, [0.0, 1.0, 0.0]),
            AbaqusNode::new(4, [0.0, 0.0, 1.0]),
            AbaqusNode::new(5, [1.0, 1.0, 0.0]),
        ];
        mesh.elements = vec![
            AbaqusElement::new(1, ElementType::C3D4, vec![1, 2, 3, 4]),
            AbaqusElement::new(2, ElementType::T3D2, vec![4, 5]),
        ];
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        assert_eq!(parsed.elements.len(), 2);
    }

    #[test]
    fn test_write_contains_heading() {
        let path = "/tmp/oxiphysics_abaqus_heading.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("*Heading"), "no *Heading found");
    }

    #[test]
    fn test_write_contains_node_keyword() {
        let path = "/tmp/oxiphysics_abaqus_nkw.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("*Node"), "no *Node found");
    }

    #[test]
    fn test_write_contains_element_keyword() {
        let path = "/tmp/oxiphysics_abaqus_ekw.inp";
        let mesh = sample_mesh();
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("*Element"), "no *Element found");
    }

    #[test]
    fn test_write_material() {
        let path = "/tmp/oxiphysics_abaqus_mat.inp";
        let mut mesh = sample_mesh();
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0));
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("*Material"), "no *Material found");
        assert!(content.contains("Steel"));
    }

    #[test]
    fn test_write_plastic_material() {
        let path = "/tmp/oxiphysics_abaqus_plastic.inp";
        let mut mesh = sample_mesh();
        mesh.materials.push(AbaqusMaterial::new_plastic(
            "Steel", 210e9, 0.3, 7800.0, 250e6, 1e9,
        ));
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("*Plastic"), "no *Plastic found");
    }

    #[test]
    fn test_write_section() {
        let path = "/tmp/oxiphysics_abaqus_sec.inp";
        let mut mesh = sample_mesh();
        mesh.sections
            .push(AbaqusSection::new("SEC1", "Steel", vec![1]));
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("*Solid Section"), "no *Solid Section");
    }

    #[test]
    fn test_write_bc_encastre() {
        let path = "/tmp/oxiphysics_abaqus_bc_enc.inp";
        let mut mesh = sample_mesh();
        mesh.boundary_conditions.push(BoundaryCondition::Encastre {
            node_set: "FIXED".to_string(),
        });
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("ENCASTRE"), "no ENCASTRE");
    }

    #[test]
    fn test_write_bc_pinned() {
        let path = "/tmp/oxiphysics_abaqus_bc_pin.inp";
        let mut mesh = sample_mesh();
        mesh.boundary_conditions.push(BoundaryCondition::Pinned {
            node_set: "PINSET".to_string(),
        });
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("PINNED"), "no PINNED");
    }

    #[test]
    fn test_write_bc_symmetry() {
        let path = "/tmp/oxiphysics_abaqus_bc_sym.inp";
        let mut mesh = sample_mesh();
        mesh.boundary_conditions
            .push(BoundaryCondition::SymmetryPlane {
                node_set: "SYMSET".to_string(),
                axis: 2,
            });
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.contains("YSYMM"), "no YSYMM");
    }

    #[test]
    fn test_parse_empty_file() {
        let path = "/tmp/oxiphysics_abaqus_empty.inp";
        std::fs::write(path, "** empty\n").unwrap();
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        assert!(parsed.nodes.is_empty());
        assert!(parsed.elements.is_empty());
    }

    #[test]
    fn test_parse_missing_file() {
        let result = AbaqusReader::new().parse("/tmp/does_not_exist_oxiphysics.inp");
        assert!(result.is_err());
    }

    #[test]
    fn test_abaqus_mesh_default() {
        let m = AbaqusMesh::default();
        assert!(m.nodes.is_empty());
        assert!(m.elements.is_empty());
    }

    #[test]
    fn test_section_new() {
        let s = AbaqusSection::new("S1", "Mat1", vec![1, 2, 3]);
        assert_eq!(s.name, "S1");
        assert_eq!(s.elements, vec![1, 2, 3]);
    }

    #[test]
    fn test_abaqus_element_new() {
        let e = AbaqusElement::new(5, ElementType::S4, vec![10, 11, 12, 13]);
        assert_eq!(e.id, 5);
        assert_eq!(e.element_type, ElementType::S4);
        assert_eq!(e.node_ids.len(), 4);
    }

    #[test]
    fn test_large_mesh_roundtrip() {
        let path = "/tmp/oxiphysics_abaqus_large.inp";
        let mut mesh = AbaqusMesh::new();
        // 100 nodes
        for i in 1..=100 {
            mesh.nodes.push(AbaqusNode::new(i, [i as f64, 0.0, 0.0]));
        }
        // 24 tetrahedral elements
        for i in 0..24usize {
            mesh.elements.push(AbaqusElement::new(
                i + 1,
                ElementType::C3D4,
                vec![
                    i * 4 % 97 + 1,
                    i * 4 % 97 + 2,
                    i * 4 % 97 + 3,
                    i * 4 % 97 + 4,
                ],
            ));
        }
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        assert_eq!(parsed.nodes.len(), 100);
        assert_eq!(parsed.elements.len(), 24);
    }

    #[test]
    fn test_node_roundtrip_precision() {
        let path = "/tmp/oxiphysics_abaqus_prec.inp";
        let mut mesh = AbaqusMesh::new();
        mesh.nodes.push(AbaqusNode::new(
            1,
            [1.23456789012345, -9.87654321098765, 2.89793238462643],
        ));
        AbaqusWriter::new().write(&mesh, path).expect("write");
        let parsed = AbaqusReader::new().parse(path).expect("parse");
        let c = parsed.nodes[0].coordinates;
        assert!((c[0] - 1.23456789012345).abs() < 1e-10);
        assert!((c[1] - (-9.87654321098765)).abs() < 1e-10);
        assert!((c[2] - 2.89793238462643).abs() < 1e-10);
    }
}
