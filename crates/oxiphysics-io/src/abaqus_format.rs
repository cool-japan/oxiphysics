// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Abaqus INP file format I/O.
//!
//! Provides readers and writers for the Abaqus `.inp` finite-element input
//! format, covering node/element definitions, material properties, section
//! assignments, and boundary conditions.

use std::collections::{BTreeMap, HashMap, HashSet};
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
    pub fn from_keyword(s: &str) -> Self {
        Self::from(s)
    }
}

impl From<&str> for ElementType {
    fn from(s: &str) -> Self {
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
///
/// The writer emits an element set named after the section and a section
/// definition that refers to it: a `*Solid Section` for continuum and truss
/// elements, a `*Shell Section` (which needs [`AbaqusSection::thickness`]) for
/// [`ElementType::S4`] shell elements.
#[derive(Debug, Clone)]
pub struct AbaqusSection {
    /// Section name, also the name of the element set written for it.
    pub name: String,
    /// Name of the associated material.
    pub material_name: String,
    /// Ids of the elements assigned to this section.
    pub elements: Vec<usize>,
    /// Shell thickness: the data line of the `*Shell Section` written for a
    /// section of [`ElementType::S4`] elements.
    ///
    /// It is required (finite and positive) for such a section and must be
    /// `None` for every other section, whose data line stays empty.
    pub thickness: Option<f64>,
}

impl AbaqusSection {
    /// Create a new section without a shell thickness.
    pub fn new(
        name: impl Into<String>,
        material_name: impl Into<String>,
        elements: Vec<usize>,
    ) -> Self {
        Self {
            name: name.into(),
            material_name: material_name.into(),
            elements,
            thickness: None,
        }
    }

    /// Set the shell thickness (see [`AbaqusSection::thickness`]).
    #[must_use]
    pub fn with_thickness(mut self, thickness: f64) -> Self {
        self.thickness = Some(thickness);
        self
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
    /// Named node sets (`*Nset`): set name to 1-based node ids.
    ///
    /// The boundary conditions refer to node sets by name, so every set a
    /// boundary condition uses must be defined here.  The writer emits one
    /// `*Nset` block per entry before the `*Boundary` block.
    pub node_sets: BTreeMap<String, Vec<usize>>,
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
    /// # Errors
    ///
    /// Returns an error when the mesh cannot be written as a valid deck (see
    /// [`AbaqusWriter::write_string`]) or when the file cannot be created or
    /// written.
    pub fn write(&self, mesh: &AbaqusMesh, path: &str) -> Result<(), IoError> {
        let text = self.write_string(mesh)?;
        fs::write(path, text).map_err(IoError::Io)
    }

    /// Serialize `mesh` into Abaqus INP text.
    ///
    /// The deck lists the nodes, the elements (grouped by type, at most 16
    /// entries per data line, continuation lines end with a comma), the node
    /// sets, the materials, one `*Elset` plus section definition per section
    /// (`*Solid Section`, or `*Shell Section` for [`ElementType::S4`] shells)
    /// and the boundary conditions.  A linearly hardening material is written
    /// as the two-row `*Plastic` table `(sy, 0)`, `(sy + 0.1 H, 0.1)`.
    ///
    /// # Errors
    ///
    /// The writer refuses to produce a deck that Abaqus would reject and
    /// returns [`IoError::General`] when
    ///
    /// - a name (element type, material, section, node set, boundary node set)
    ///   is empty, has blanks at either end, starts with `*` or contains a
    ///   comma, an equals sign or a control character,
    /// - a node or element id is zero or used twice, a coordinate or a
    ///   material property is not finite,
    /// - a material, a section or a node set is defined twice (names are
    ///   compared without regard to case),
    /// - an element has no nodes or refers to a node that is not in the mesh,
    /// - a node set or a section is empty or lists an id that is not in the
    ///   mesh,
    /// - a section uses a material that is not defined, or an element belongs
    ///   to two sections,
    /// - a section mixes [`ElementType::S4`] shells with other elements, a
    ///   shell section has no valid thickness, or another section has one,
    /// - a boundary condition names a node set that is not in
    ///   [`AbaqusMesh::node_sets`], or a symmetry axis is not 1, 2 or 3.
    pub fn write_string(&self, mesh: &AbaqusMesh) -> Result<String, IoError> {
        let section_kinds = Self::validate(mesh)?;
        let mut buf = String::new();

        writeln!(buf, "** Generated by OxiPhysics AbaqusWriter").map_err(fmt_error)?;
        writeln!(buf, "*Heading").map_err(fmt_error)?;
        writeln!(buf, "OxiPhysics model").map_err(fmt_error)?;

        // --- Nodes ---
        writeln!(buf, "*Node").map_err(fmt_error)?;
        for node in &mesh.nodes {
            writeln!(
                buf,
                "{}, {:.15e}, {:.15e}, {:.15e}",
                node.id, node.coordinates[0], node.coordinates[1], node.coordinates[2]
            )
            .map_err(fmt_error)?;
        }

        // --- Elements (grouped by type) ---
        let mut types_seen: Vec<&str> = Vec::new();
        for el in &mesh.elements {
            let t = el.element_type.as_str();
            if !types_seen.contains(&t) {
                types_seen.push(t);
            }
        }
        for etype in types_seen {
            writeln!(buf, "*Element, type={etype}").map_err(fmt_error)?;
            for el in mesh
                .elements
                .iter()
                .filter(|e| e.element_type.as_str() == etype)
            {
                push_element(&mut buf, el)?;
            }
        }

        // --- Node sets (referenced by the boundary conditions) ---
        for (name, ids) in &mesh.node_sets {
            writeln!(buf, "*Nset, nset={name}").map_err(fmt_error)?;
            push_id_lines(&mut buf, ids)?;
        }

        // --- Materials ---
        for mat in &mesh.materials {
            writeln!(buf, "*Material, name={}", mat.name).map_err(fmt_error)?;
            writeln!(buf, "*Density").map_err(fmt_error)?;
            writeln!(buf, "{:.15e}", mat.density).map_err(fmt_error)?;
            writeln!(buf, "*Elastic").map_err(fmt_error)?;
            writeln!(
                buf,
                "{:.15e}, {:.15e}",
                mat.elastic.young_modulus, mat.elastic.poisson_ratio
            )
            .map_err(fmt_error)?;
            if let Some(p) = &mat.plastic {
                // Rows of (yield stress, plastic strain); the stress grows
                // linearly with the hardening modulus between them.
                writeln!(buf, "*Plastic").map_err(fmt_error)?;
                writeln!(buf, "{:.15e}, {:.15e}", p.yield_stress, 0.0).map_err(fmt_error)?;
                writeln!(
                    buf,
                    "{:.15e}, {:.15e}",
                    p.yield_stress + p.hardening_modulus * PLASTIC_TABLE_STRAIN,
                    PLASTIC_TABLE_STRAIN
                )
                .map_err(fmt_error)?;
            }
        }

        // --- Sections: an element set named after the section, then the section ---
        for (sec, kind) in mesh.sections.iter().zip(&section_kinds) {
            writeln!(buf, "*Elset, elset={}", sec.name).map_err(fmt_error)?;
            push_id_lines(&mut buf, &sec.elements)?;
            match kind {
                SectionKind::Solid => {
                    writeln!(
                        buf,
                        "*Solid Section, elset={}, material={}",
                        sec.name, sec.material_name
                    )
                    .map_err(fmt_error)?;
                    // The data line of a solid section is empty (a thickness
                    // only matters for planar sections).
                    writeln!(buf, ",").map_err(fmt_error)?;
                }
                SectionKind::Shell { thickness } => {
                    writeln!(
                        buf,
                        "*Shell Section, elset={}, material={}",
                        sec.name, sec.material_name
                    )
                    .map_err(fmt_error)?;
                    writeln!(buf, "{thickness:.15e}").map_err(fmt_error)?;
                }
            }
        }

        // --- Boundary Conditions ---
        if !mesh.boundary_conditions.is_empty() {
            writeln!(buf, "*Boundary").map_err(fmt_error)?;
        }
        for bc in &mesh.boundary_conditions {
            match bc {
                BoundaryCondition::Encastre { node_set } => {
                    writeln!(buf, "{node_set}, ENCASTRE").map_err(fmt_error)?;
                }
                BoundaryCondition::Pinned { node_set } => {
                    writeln!(buf, "{node_set}, PINNED").map_err(fmt_error)?;
                }
                BoundaryCondition::SymmetryPlane { node_set, axis } => {
                    let sym = match axis {
                        1 => "XSYMM",
                        2 => "YSYMM",
                        3 => "ZSYMM",
                        other => {
                            return Err(IoError::General(format!(
                                "Abaqus symmetry axis must be 1, 2 or 3, got {other}"
                            )));
                        }
                    };
                    writeln!(buf, "{node_set}, {sym}").map_err(fmt_error)?;
                }
            }
        }

        Ok(buf)
    }

    /// Check that `mesh` can be written as a self-consistent deck and decide
    /// which kind of section definition each of its sections needs.
    fn validate(mesh: &AbaqusMesh) -> Result<Vec<SectionKind>, IoError> {
        let mut node_ids: HashSet<usize> = HashSet::with_capacity(mesh.nodes.len());
        for node in &mesh.nodes {
            check_id("node", node.id)?;
            if let Some(coordinate) = node.coordinates.iter().find(|c| !c.is_finite()) {
                return Err(IoError::General(format!(
                    "Abaqus node {} has the non-finite coordinate {coordinate}",
                    node.id
                )));
            }
            if !node_ids.insert(node.id) {
                return Err(IoError::General(format!(
                    "Abaqus node {} is defined twice",
                    node.id
                )));
            }
        }

        let mut element_types: HashMap<usize, &ElementType> =
            HashMap::with_capacity(mesh.elements.len());
        for element in &mesh.elements {
            check_id("element", element.id)?;
            check_name("element type", element.element_type.as_str())?;
            if element_types
                .insert(element.id, &element.element_type)
                .is_some()
            {
                return Err(IoError::General(format!(
                    "Abaqus element {} is defined twice",
                    element.id
                )));
            }
            if element.node_ids.is_empty() {
                return Err(IoError::General(format!(
                    "Abaqus element {} has no nodes",
                    element.id
                )));
            }
            if let Some(missing) = element.node_ids.iter().find(|n| !node_ids.contains(*n)) {
                return Err(IoError::General(format!(
                    "Abaqus element {} refers to node {missing}, which is not in the mesh",
                    element.id
                )));
            }
        }

        // Abaqus compares names case-insensitively.
        let mut node_set_names: HashSet<String> = HashSet::with_capacity(mesh.node_sets.len());
        for (name, ids) in &mesh.node_sets {
            check_name("node set", name)?;
            if !node_set_names.insert(name.to_ascii_uppercase()) {
                return Err(IoError::General(format!(
                    "Abaqus node set {name} has the same name as another node set up to case"
                )));
            }
            if ids.is_empty() {
                return Err(IoError::General(format!("Abaqus node set {name} is empty")));
            }
            if let Some(missing) = ids.iter().find(|n| !node_ids.contains(*n)) {
                return Err(IoError::General(format!(
                    "Abaqus node set {name} lists node {missing}, which is not in the mesh"
                )));
            }
        }

        let mut material_names: HashSet<String> = HashSet::with_capacity(mesh.materials.len());
        for mat in &mesh.materials {
            check_name("material", &mat.name)?;
            check_material_numbers(mat)?;
            if !material_names.insert(mat.name.to_ascii_uppercase()) {
                return Err(IoError::General(format!(
                    "Abaqus material {} is defined twice (names are case-insensitive)",
                    mat.name
                )));
            }
        }

        let mut section_names: HashSet<String> = HashSet::with_capacity(mesh.sections.len());
        let mut section_of_element: HashMap<usize, &str> = HashMap::new();
        let mut kinds = Vec::with_capacity(mesh.sections.len());
        for sec in &mesh.sections {
            check_name("section", &sec.name)?;
            check_name("section material", &sec.material_name)?;
            if !section_names.insert(sec.name.to_ascii_uppercase()) {
                return Err(IoError::General(format!(
                    "Abaqus section {} is defined twice (names are case-insensitive)",
                    sec.name
                )));
            }
            if !material_names.contains(&sec.material_name.to_ascii_uppercase()) {
                return Err(IoError::General(format!(
                    "Abaqus section {} uses material {}, which is not defined",
                    sec.name, sec.material_name
                )));
            }
            if sec.elements.is_empty() {
                return Err(IoError::General(format!(
                    "Abaqus section {} has no elements",
                    sec.name
                )));
            }
            let mut shells = 0usize;
            for id in &sec.elements {
                let Some(element_type) = element_types.get(id) else {
                    return Err(IoError::General(format!(
                        "Abaqus section {} lists element {id}, which is not in the mesh",
                        sec.name
                    )));
                };
                if matches!(element_type, ElementType::S4) {
                    shells += 1;
                }
                // Section names are unique, so a different name is a different
                // section; the same element may repeat within one section.
                if let Some(other) = section_of_element.insert(*id, sec.name.as_str())
                    && !other.eq_ignore_ascii_case(&sec.name)
                {
                    return Err(IoError::General(format!(
                        "Abaqus element {id} is assigned to the sections {other} and {}",
                        sec.name
                    )));
                }
            }
            kinds.push(classify_section(sec, shells)?);
        }

        for bc in &mesh.boundary_conditions {
            let node_set = match bc {
                BoundaryCondition::Encastre { node_set }
                | BoundaryCondition::Pinned { node_set }
                | BoundaryCondition::SymmetryPlane { node_set, .. } => node_set,
            };
            check_name("boundary node set", node_set)?;
            if !node_set_names.contains(&node_set.to_ascii_uppercase()) {
                return Err(IoError::General(format!(
                    "Abaqus boundary condition refers to node set {node_set}, which is not in AbaqusMesh::node_sets"
                )));
            }
        }

        Ok(kinds)
    }
}

/// Section definition a section is written as.
#[derive(Debug, Clone, Copy, PartialEq)]
enum SectionKind {
    /// `*Solid Section` with an empty data line (continuum and truss
    /// elements).
    Solid,
    /// `*Shell Section` with the given thickness on its data line.
    Shell { thickness: f64 },
}

/// Decide the kind of `sec`, which holds `shells` elements of type
/// [`ElementType::S4`] (counting repeated ids).
fn classify_section(sec: &AbaqusSection, shells: usize) -> Result<SectionKind, IoError> {
    if shells == 0 {
        return match sec.thickness {
            None => Ok(SectionKind::Solid),
            Some(_) => Err(IoError::General(format!(
                "Abaqus section {} has a thickness but no S4 shell elements; only shell sections take one",
                sec.name
            ))),
        };
    }
    if shells != sec.elements.len() {
        return Err(IoError::General(format!(
            "Abaqus section {} mixes S4 shell elements with other elements; a section is either a shell section or a solid section",
            sec.name
        )));
    }
    match sec.thickness {
        Some(thickness) if thickness.is_finite() && thickness > 0.0 => {
            Ok(SectionKind::Shell { thickness })
        }
        Some(thickness) => Err(IoError::General(format!(
            "Abaqus section {} has the invalid shell thickness {thickness}",
            sec.name
        ))),
        None => Err(IoError::General(format!(
            "Abaqus section {} holds S4 shell elements, which need a thickness (AbaqusSection::with_thickness)",
            sec.name
        ))),
    }
}

/// Reject an id that Abaqus would not accept: ids are positive integers.
fn check_id(kind: &str, id: usize) -> Result<(), IoError> {
    if id == 0 {
        return Err(IoError::General(format!(
            "Abaqus {kind} ids start at 1, got 0"
        )));
    }
    Ok(())
}

/// Reject a material whose numbers cannot be written to a deck.
fn check_material_numbers(mat: &AbaqusMaterial) -> Result<(), IoError> {
    let plastic = mat
        .plastic
        .iter()
        .flat_map(|p| [p.yield_stress, p.hardening_modulus]);
    let numbers = [
        mat.elastic.young_modulus,
        mat.elastic.poisson_ratio,
        mat.density,
    ]
    .into_iter()
    .chain(plastic);
    for number in numbers {
        if !number.is_finite() {
            return Err(IoError::General(format!(
                "Abaqus material {} has the non-finite property value {number}",
                mat.name
            )));
        }
    }
    Ok(())
}

/// Maximum number of entries on one Abaqus data line.
const MAX_ENTRIES_PER_LINE: usize = 16;

/// Plastic strain of the second row of the `*Plastic` table written for a
/// linearly hardening material.
const PLASTIC_TABLE_STRAIN: f64 = 0.1;

fn fmt_error(error: std::fmt::Error) -> IoError {
    IoError::General(format!("formatting the Abaqus deck failed: {error}"))
}

/// Reject a name that would break a keyword or a data line: an empty one, one
/// with blanks at either end (Abaqus trims them), one that starts with `*`
/// (the data line `*FIXED, ENCASTRE` would read as a keyword) and one with a
/// comma, an equals sign or a control character.
fn check_name(kind: &str, name: &str) -> Result<(), IoError> {
    if name.trim().is_empty()
        || name.trim() != name
        || name.starts_with('*')
        || name.chars().any(|c| c == ',' || c == '=' || c.is_control())
    {
        return Err(IoError::General(format!(
            "Abaqus {kind} name {name:?} is empty, has blanks at either end, starts with `*` or contains a comma, an equals sign or a control character"
        )));
    }
    Ok(())
}

/// Append `ids` as data lines of at most [`MAX_ENTRIES_PER_LINE`] entries.
fn push_id_lines(buf: &mut String, ids: &[usize]) -> Result<(), IoError> {
    for chunk in ids.chunks(MAX_ENTRIES_PER_LINE) {
        let entries: Vec<String> = chunk.iter().map(usize::to_string).collect();
        writeln!(buf, "{}", entries.join(", ")).map_err(fmt_error)?;
    }
    Ok(())
}

/// Append one element definition.
///
/// The element id and its nodes are written with at most
/// [`MAX_ENTRIES_PER_LINE`] entries per line; every line except the last ends
/// with a comma, which tells Abaqus that the definition continues.
fn push_element(buf: &mut String, element: &AbaqusElement) -> Result<(), IoError> {
    let entries: Vec<String> = std::iter::once(element.id)
        .chain(element.node_ids.iter().copied())
        .map(|v| v.to_string())
        .collect();
    let line_count = entries.len().div_ceil(MAX_ENTRIES_PER_LINE);
    for (index, chunk) in entries.chunks(MAX_ENTRIES_PER_LINE).enumerate() {
        let continued = if index + 1 < line_count { "," } else { "" };
        writeln!(buf, "{}{continued}", chunk.join(", ")).map_err(fmt_error)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// AbaqusReader
// ---------------------------------------------------------------------------

/// Parses a subset of Abaqus `.inp` files: `*NODE` and `*ELEMENT` sections.
///
/// Element definitions that continue on the next data line (a line ending in
/// a comma, as written for elements with more than 15 nodes) are joined.
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
        // Data lines of an element definition that has not ended yet.
        let mut pending_element = String::new();

        for line_res in reader.lines() {
            let line = line_res.map_err(IoError::Io)?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with("**") {
                continue;
            }

            if trimmed.starts_with('*') {
                Self::flush_element(&mut pending_element, &current_block, &mut mesh);
                // New keyword
                let upper = trimmed.to_uppercase();
                if upper.starts_with("*NODE") && !upper.starts_with("*NSET") {
                    current_block = Block::Node;
                } else if upper.starts_with("*ELEMENT") {
                    // Extract type= parameter
                    let etype = Self::extract_param(trimmed, "TYPE")
                        .unwrap_or_else(|| "UNKNOWN".to_string());
                    current_block = Block::Element(ElementType::from_keyword(&etype));
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
                Block::Element(_) => {
                    if !pending_element.is_empty() {
                        pending_element.push(' ');
                    }
                    pending_element.push_str(trimmed);
                    // A trailing comma continues the definition on the next line.
                    if !trimmed.ends_with(',') {
                        Self::flush_element(&mut pending_element, &current_block, &mut mesh);
                    }
                }
                Block::None => {}
            }
        }
        Self::flush_element(&mut pending_element, &current_block, &mut mesh);

        Ok(mesh)
    }

    /// Parse the joined data lines collected in `pending` as one element of
    /// the current block and store it in `mesh`.
    fn flush_element(pending: &mut String, block: &Block, mesh: &mut AbaqusMesh) {
        if let Block::Element(etype) = block
            && let Some(el) = Self::parse_element_line(pending, etype.clone())
        {
            mesh.elements.push(el);
        }
        pending.clear();
    }

    /// Extract value of `key=value` from a keyword line (case-insensitive).
    ///
    /// Only ASCII letters are folded, so a byte offset found in the folded
    /// copy is also valid in `line`, whatever other characters it contains.
    fn extract_param(line: &str, key: &str) -> Option<String> {
        let upper = line.to_ascii_uppercase();
        let key_eq = format!("{}=", key.to_ascii_uppercase());
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
        assert_eq!(ElementType::from_keyword("C3D4"), ElementType::C3D4);
    }

    #[test]
    fn test_element_type_from_str_case_insensitive() {
        assert_eq!(ElementType::from_keyword("c3d8"), ElementType::C3D8);
    }

    #[test]
    fn test_element_type_from_str_s4() {
        assert_eq!(ElementType::from_keyword("S4"), ElementType::S4);
    }

    #[test]
    fn test_element_type_from_str_t3d2() {
        assert_eq!(ElementType::from_keyword("T3D2"), ElementType::T3D2);
    }

    #[test]
    fn test_element_type_from_str_unknown() {
        match ElementType::from_keyword("FOOBAR") {
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
        let path = std::env::temp_dir().join("oxiphysics_abaqus_test_write.inp");
        let mesh = sample_mesh();
        let writer = AbaqusWriter::new();
        writer
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write failed");
        assert!(path.exists());
    }

    #[test]
    fn test_roundtrip_node_count() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_roundtrip.inp");
        let mesh = sample_mesh();
        let writer = AbaqusWriter::new();
        writer
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write failed");

        let reader = AbaqusReader::new();
        let parsed = reader
            .parse(path.to_str().unwrap_or(""))
            .expect("parse failed");
        assert_eq!(parsed.nodes.len(), 4);
    }

    #[test]
    fn test_roundtrip_element_count() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_rt_elem.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements.len(), 1);
    }

    #[test]
    fn test_roundtrip_node_ids() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_rt_nodeids.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        let ids: Vec<usize> = parsed.nodes.iter().map(|n| n.id).collect();
        assert_eq!(ids, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_roundtrip_node_coordinates() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_rt_coords.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        let n1 = &parsed.nodes[0];
        assert!((n1.coordinates[0]).abs() < 1e-10);
        let n2 = &parsed.nodes[1];
        assert!((n2.coordinates[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_roundtrip_element_type() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_rt_etype.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements[0].element_type, ElementType::C3D4);
    }

    #[test]
    fn test_roundtrip_element_nodes() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_rt_enodes.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements[0].node_ids, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_roundtrip_element_id() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_rt_eid.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements[0].id, 1);
    }

    #[test]
    fn test_multiple_element_types() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_multtype.inp");
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
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements.len(), 2);
    }

    #[test]
    fn test_write_contains_heading() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_heading.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("*Heading"), "no *Heading found");
    }

    #[test]
    fn test_write_contains_node_keyword() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_nkw.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("*Node"), "no *Node found");
    }

    #[test]
    fn test_write_contains_element_keyword() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_ekw.inp");
        let mesh = sample_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("*Element"), "no *Element found");
    }

    #[test]
    fn test_write_material() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_mat.inp");
        let mut mesh = sample_mesh();
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0));
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("*Material"), "no *Material found");
        assert!(content.contains("Steel"));
    }

    #[test]
    fn test_write_plastic_material() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_plastic.inp");
        let mut mesh = sample_mesh();
        mesh.materials.push(AbaqusMaterial::new_plastic(
            "Steel", 210e9, 0.3, 7800.0, 250e6, 1e9,
        ));
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("*Plastic"), "no *Plastic found");
    }

    #[test]
    fn test_write_section() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_sec.inp");
        let mut mesh = sample_mesh();
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0));
        mesh.sections
            .push(AbaqusSection::new("SEC1", "Steel", vec![1]));
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("*Solid Section"), "no *Solid Section");
        assert!(content.contains("*Elset, elset=SEC1"), "no *Elset");
    }

    #[test]
    fn test_write_bc_encastre() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_bc_enc.inp");
        let mut mesh = sample_mesh();
        mesh.node_sets.insert("FIXED".to_string(), vec![1, 2]);
        mesh.boundary_conditions.push(BoundaryCondition::Encastre {
            node_set: "FIXED".to_string(),
        });
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("FIXED, ENCASTRE"), "no ENCASTRE");
    }

    #[test]
    fn test_write_bc_pinned() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_bc_pin.inp");
        let mut mesh = sample_mesh();
        mesh.node_sets.insert("PINSET".to_string(), vec![3]);
        mesh.boundary_conditions.push(BoundaryCondition::Pinned {
            node_set: "PINSET".to_string(),
        });
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("PINSET, PINNED"), "no PINNED");
    }

    #[test]
    fn test_write_bc_symmetry() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_bc_sym.inp");
        let mut mesh = sample_mesh();
        mesh.node_sets.insert("SYMSET".to_string(), vec![1, 3, 4]);
        mesh.boundary_conditions
            .push(BoundaryCondition::SymmetryPlane {
                node_set: "SYMSET".to_string(),
                axis: 2,
            });
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("SYMSET, YSYMM"), "no YSYMM");
    }

    #[test]
    fn test_parse_empty_file() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_empty.inp");
        std::fs::write(&path, "** empty\n").unwrap();
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert!(parsed.nodes.is_empty());
        assert!(parsed.elements.is_empty());
    }

    #[test]
    fn test_parse_missing_file() {
        let path = std::env::temp_dir().join("does_not_exist_oxiphysics.inp");
        let result = AbaqusReader::new().parse(path.to_str().unwrap_or(""));
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
        let path = std::env::temp_dir().join("oxiphysics_abaqus_large.inp");
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
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.nodes.len(), 100);
        assert_eq!(parsed.elements.len(), 24);
    }

    #[test]
    fn test_node_roundtrip_precision() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_prec.inp");
        let mut mesh = AbaqusMesh::new();
        mesh.nodes.push(AbaqusNode::new(
            1,
            [1.23456789012345, -9.87654321098765, 2.89793238462643],
        ));
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        let c = parsed.nodes[0].coordinates;
        assert!((c[0] - 1.23456789012345).abs() < 1e-10);
        assert!((c[1] - (-9.87654321098765)).abs() < 1e-10);
        assert!((c[2] - 2.89793238462643).abs() < 1e-10);
    }

    // --- Valid deck structure ---

    /// A mesh that uses every feature of the writer: two element types, an
    /// elastic and a plastic material, two sections, node sets and one boundary
    /// condition of each kind.
    fn full_mesh() -> AbaqusMesh {
        let mut mesh = sample_mesh();
        mesh.nodes.push(AbaqusNode::new(5, [1.0, 1.0, 0.0]));
        mesh.elements
            .push(AbaqusElement::new(2, ElementType::T3D2, vec![4, 5]));
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Rubber", 1.0e7, 0.49, 1100.0));
        mesh.materials.push(AbaqusMaterial::new_plastic(
            "Steel", 210e9, 0.3, 7800.0, 250e6, 1e9,
        ));
        mesh.sections
            .push(AbaqusSection::new("SOLIDS", "Steel", vec![1]));
        mesh.sections
            .push(AbaqusSection::new("BARS", "Rubber", vec![2]));
        mesh.node_sets.insert("FIXED".to_string(), vec![1, 2, 3]);
        mesh.node_sets.insert("MIRROR".to_string(), vec![4, 5]);
        mesh.boundary_conditions.push(BoundaryCondition::Encastre {
            node_set: "FIXED".to_string(),
        });
        mesh.boundary_conditions
            .push(BoundaryCondition::SymmetryPlane {
                node_set: "MIRROR".to_string(),
                axis: 3,
            });
        mesh
    }

    fn deck_lines(text: &str) -> Vec<&str> {
        text.lines().collect()
    }

    fn position_of(lines: &[&str], needle: &str) -> usize {
        lines
            .iter()
            .position(|l| *l == needle)
            .unwrap_or_else(|| panic!("line {needle:?} not found in\n{}", lines.join("\n")))
    }

    #[test]
    fn test_deck_defines_every_set_before_it_is_used() {
        let text = AbaqusWriter::new().write_string(&full_mesh()).unwrap();
        let lines = deck_lines(&text);
        let boundary = position_of(&lines, "*Boundary");
        for set in ["FIXED", "MIRROR"] {
            let nset = position_of(&lines, &format!("*Nset, nset={set}"));
            assert!(nset < boundary, "{set} defined after *Boundary:\n{text}");
        }
        assert_eq!(
            lines[position_of(&lines, "*Nset, nset=FIXED") + 1],
            "1, 2, 3"
        );
        assert_eq!(lines[position_of(&lines, "*Nset, nset=MIRROR") + 1], "4, 5");
        assert_eq!(lines[boundary + 1], "FIXED, ENCASTRE");
        assert_eq!(lines[boundary + 2], "MIRROR, ZSYMM");
    }

    #[test]
    fn test_deck_section_references_an_elset_and_has_an_empty_data_line() {
        let text = AbaqusWriter::new().write_string(&full_mesh()).unwrap();
        let lines = deck_lines(&text);
        let elset = position_of(&lines, "*Elset, elset=SOLIDS");
        // The element ids are the data of the *Elset, not of the section.
        assert_eq!(lines[elset + 1], "1");
        assert_eq!(
            lines[elset + 2],
            "*Solid Section, elset=SOLIDS, material=Steel"
        );
        assert_eq!(lines[elset + 3], ",");
        let bars = position_of(&lines, "*Elset, elset=BARS");
        assert_eq!(lines[bars + 1], "2");
        assert_eq!(
            lines[bars + 2],
            "*Solid Section, elset=BARS, material=Rubber"
        );
        assert_eq!(lines[bars + 3], ",");
    }

    #[test]
    fn test_deck_has_no_end_part_without_a_part() {
        let text = AbaqusWriter::new().write_string(&full_mesh()).unwrap();
        assert!(!text.to_uppercase().contains("*END PART"), "{text}");
        assert!(!text.to_uppercase().contains("*PART"), "{text}");
    }

    #[test]
    fn test_plastic_table_is_stress_versus_plastic_strain() {
        let text = AbaqusWriter::new().write_string(&full_mesh()).unwrap();
        let lines = deck_lines(&text);
        let plastic = position_of(&lines, "*Plastic");
        let row = |index: usize| -> Vec<f64> {
            lines[index]
                .split(',')
                .map(|v| v.trim().parse::<f64>().unwrap())
                .collect()
        };
        assert_eq!(row(plastic + 1), vec![250e6, 0.0]);
        let second = row(plastic + 2);
        assert!((second[0] - (250e6 + 1e9 * 0.1)).abs() < 1.0, "{second:?}");
        assert_eq!(second[1], 0.1);
        // The table has exactly the two rows: the next line is a keyword.
        assert!(
            lines[plastic + 3].starts_with('*'),
            "{}",
            lines[plastic + 3]
        );
        // The purely elastic material has no table.
        assert_eq!(text.matches("*Plastic").count(), 1);
    }

    #[test]
    fn test_reader_reparses_nodes_and_elements_of_a_full_deck() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_full_deck.inp");
        let mesh = full_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.nodes, mesh.nodes);
        // The writer groups the elements by type; the mesh was built in that
        // order already.
        assert_eq!(parsed.elements, mesh.elements);
    }

    #[test]
    fn test_no_data_line_has_more_than_sixteen_entries() {
        let mut mesh = AbaqusMesh::new();
        for id in 1..=40 {
            mesh.nodes.push(AbaqusNode::new(id, [id as f64, 0.0, 0.0]));
        }
        mesh.elements.push(AbaqusElement::new(
            1,
            ElementType::from_keyword("C3D20"),
            (1..=20).collect(),
        ));
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0));
        mesh.sections
            .push(AbaqusSection::new("ALL", "Steel", vec![1]));
        mesh.node_sets
            .insert("MANY".to_string(), (1..=40).collect());
        mesh.boundary_conditions.push(BoundaryCondition::Pinned {
            node_set: "MANY".to_string(),
        });
        let text = AbaqusWriter::new().write_string(&mesh).unwrap();
        for line in text.lines() {
            if line.starts_with('*') || line.starts_with("**") {
                continue;
            }
            let entries = line.split(',').filter(|s| !s.trim().is_empty()).count();
            assert!(entries <= 16, "{entries} entries in {line:?}");
        }
        let lines = deck_lines(&text);
        let nset = position_of(&lines, "*Nset, nset=MANY");
        assert_eq!(lines[nset + 1].split(',').count(), 16);
        assert_eq!(lines[nset + 2].split(',').count(), 16);
        assert_eq!(lines[nset + 3].split(',').count(), 8);
    }

    #[test]
    fn test_element_with_more_than_fifteen_nodes_is_continued_and_rejoined() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_c3d20.inp");
        let mut mesh = AbaqusMesh::new();
        for id in 1..=21 {
            mesh.nodes.push(AbaqusNode::new(id, [id as f64, 0.0, 0.0]));
        }
        let c3d20 = ElementType::from_keyword("C3D20");
        mesh.elements
            .push(AbaqusElement::new(7, c3d20.clone(), (1..=20).collect()));
        mesh.elements
            .push(AbaqusElement::new(8, c3d20, (2..=21).collect()));
        let text = AbaqusWriter::new().write_string(&mesh).unwrap();
        let lines = deck_lines(&text);
        let first = position_of(&lines, "*Element, type=C3D20") + 1;
        assert_eq!(
            lines[first],
            "7, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,"
        );
        assert_eq!(lines[first + 1], "16, 17, 18, 19, 20");
        assert!(lines[first + 2].starts_with("8, 2, 3,"));
        std::fs::write(&path, &text).unwrap();
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements, mesh.elements);
    }

    #[test]
    fn test_reader_joins_continued_element_lines_of_a_foreign_file() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_foreign_continuation.inp");
        let text = "*Node\n1, 0, 0, 0\n*Element, type=C3D20\n1, 1, 2, 3, 4, 5, 6, 7, 8,\n  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20\n2, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6,\n5, 4, 3, 2,\n";
        std::fs::write(&path, text).unwrap();
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements.len(), 2);
        assert_eq!(parsed.elements[0].node_ids, (1..=20).collect::<Vec<_>>());
        assert_eq!(
            parsed.elements[1].node_ids,
            (2..=20).rev().collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_node_sets_and_sections_do_not_disturb_the_reader() {
        // The data lines of *Nset, *Elset, *Solid Section and *Boundary must
        // not be taken for nodes or elements.
        let path = std::env::temp_dir().join("oxiphysics_abaqus_sets_ignored.inp");
        let mesh = full_mesh();
        AbaqusWriter::new()
            .write(&mesh, path.to_str().unwrap_or(""))
            .expect("write");
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.nodes.len(), 5);
        assert_eq!(parsed.elements.len(), 2);
    }

    fn write_error(mesh: &AbaqusMesh) -> String {
        match AbaqusWriter::new().write_string(mesh) {
            Ok(text) => panic!("the deck must be rejected:\n{text}"),
            Err(error) => error.to_string(),
        }
    }

    #[test]
    fn test_boundary_condition_on_an_undefined_node_set_is_rejected() {
        let mut mesh = sample_mesh();
        mesh.boundary_conditions.push(BoundaryCondition::Encastre {
            node_set: "NOWHERE".to_string(),
        });
        let message = write_error(&mesh);
        assert!(message.contains("NOWHERE"), "{message}");
        assert!(message.contains("node_sets"), "{message}");
        // The set name is matched case-insensitively like Abaqus does.
        mesh.node_sets.insert("nowhere".to_string(), vec![1]);
        assert!(AbaqusWriter::new().write_string(&mesh).is_ok());
    }

    #[test]
    fn test_section_with_undefined_material_or_missing_elements_is_rejected() {
        let mut mesh = sample_mesh();
        mesh.sections
            .push(AbaqusSection::new("SEC", "Missing", vec![1]));
        let message = write_error(&mesh);
        assert!(message.contains("Missing"), "{message}");

        mesh.materials
            .push(AbaqusMaterial::new_elastic("Missing", 1.0, 0.3, 1.0));
        assert!(AbaqusWriter::new().write_string(&mesh).is_ok());

        mesh.sections[0].elements = vec![99];
        assert!(write_error(&mesh).contains("element 99"));
        mesh.sections[0].elements.clear();
        assert!(write_error(&mesh).contains("no elements"));
    }

    #[test]
    fn test_invalid_symmetry_axis_is_rejected_instead_of_replaced() {
        let mut mesh = sample_mesh();
        mesh.node_sets.insert("S".to_string(), vec![1]);
        mesh.boundary_conditions
            .push(BoundaryCondition::SymmetryPlane {
                node_set: "S".to_string(),
                axis: 7,
            });
        assert!(write_error(&mesh).contains("axis"));
    }

    /// One `S4` shell on the unit square (nodes 1 to 4) with a material.
    fn shell_mesh() -> AbaqusMesh {
        let mut mesh = AbaqusMesh::new();
        for (id, x, y) in [(1, 0.0, 0.0), (2, 1.0, 0.0), (3, 1.0, 1.0), (4, 0.0, 1.0)] {
            mesh.nodes.push(AbaqusNode::new(id, [x, y, 0.0]));
        }
        mesh.elements
            .push(AbaqusElement::new(1, ElementType::S4, vec![1, 2, 3, 4]));
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0));
        mesh
    }

    #[test]
    fn test_names_that_read_as_keywords_or_lose_blanks_are_rejected() {
        for name in ["*FIX", " FIX", "FIX ", "**c", "  "] {
            let mut mesh = sample_mesh();
            mesh.node_sets.insert(name.to_string(), vec![1]);
            mesh.boundary_conditions.push(BoundaryCondition::Pinned {
                node_set: name.to_string(),
            });
            let message = write_error(&mesh);
            assert!(message.contains("node set name"), "{name:?}: {message}");
        }
        // Interior blanks and punctuation other than the separators are fine.
        let mut mesh = sample_mesh();
        mesh.node_sets.insert("left side-1_a".to_string(), vec![1]);
        mesh.boundary_conditions.push(BoundaryCondition::Pinned {
            node_set: "left side-1_a".to_string(),
        });
        assert!(AbaqusWriter::new().write_string(&mesh).is_ok());
    }

    #[test]
    fn test_shell_section_is_written_with_its_thickness() {
        let mut mesh = shell_mesh();
        mesh.sections
            .push(AbaqusSection::new("PLATE", "Steel", vec![1]).with_thickness(0.0025));
        let text = AbaqusWriter::new().write_string(&mesh).unwrap();
        let lines = deck_lines(&text);
        let elset = position_of(&lines, "*Elset, elset=PLATE");
        assert_eq!(lines[elset + 1], "1");
        assert_eq!(
            lines[elset + 2],
            "*Shell Section, elset=PLATE, material=Steel"
        );
        assert_eq!(lines[elset + 3].trim().parse::<f64>().unwrap(), 0.0025);
        assert!(!text.contains("*Solid Section"), "{text}");
    }

    #[test]
    fn test_shell_section_needs_a_valid_thickness_and_other_sections_take_none() {
        let mut mesh = shell_mesh();
        mesh.sections
            .push(AbaqusSection::new("PLATE", "Steel", vec![1]));
        let message = write_error(&mesh);
        assert!(message.contains("need a thickness"), "{message}");
        for thickness in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            mesh.sections[0].thickness = Some(thickness);
            assert!(write_error(&mesh).contains("invalid shell thickness"));
        }

        // A solid section with a thickness is refused too.
        let mut solid = sample_mesh();
        solid
            .materials
            .push(AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0));
        solid
            .sections
            .push(AbaqusSection::new("SOLID", "Steel", vec![1]).with_thickness(1.0));
        assert!(write_error(&solid).contains("no S4 shell elements"));
    }

    #[test]
    fn test_section_mixing_shells_and_solids_is_rejected() {
        let mut mesh = shell_mesh();
        mesh.nodes.push(AbaqusNode::new(5, [0.0, 0.0, 1.0]));
        mesh.elements
            .push(AbaqusElement::new(2, ElementType::T3D2, vec![1, 5]));
        mesh.sections
            .push(AbaqusSection::new("BOTH", "Steel", vec![1, 2]).with_thickness(0.1));
        assert!(write_error(&mesh).contains("mixes S4 shell elements"));
        // Split into one section per kind, the deck is fine.
        mesh.sections[0] = AbaqusSection::new("PLATE", "Steel", vec![1]).with_thickness(0.1);
        mesh.sections
            .push(AbaqusSection::new("BARS", "Steel", vec![2]));
        let text = AbaqusWriter::new().write_string(&mesh).unwrap();
        assert_eq!(text.matches("*Shell Section").count(), 1);
        assert_eq!(text.matches("*Solid Section").count(), 1);
    }

    #[test]
    fn test_duplicate_ids_and_names_are_rejected() {
        let mut mesh = sample_mesh();
        mesh.nodes.push(AbaqusNode::new(2, [5.0, 5.0, 5.0]));
        assert!(write_error(&mesh).contains("node 2 is defined twice"));
        mesh.nodes.pop();

        mesh.elements
            .push(AbaqusElement::new(1, ElementType::T3D2, vec![1, 2]));
        assert!(write_error(&mesh).contains("element 1 is defined twice"));
        mesh.elements.pop();

        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel", 1.0, 0.3, 1.0));
        mesh.materials
            .push(AbaqusMaterial::new_elastic("STEEL", 2.0, 0.3, 1.0));
        assert!(write_error(&mesh).contains("material STEEL is defined twice"));
        mesh.materials.pop();

        mesh.sections
            .push(AbaqusSection::new("A", "Steel", vec![1]));
        mesh.sections
            .push(AbaqusSection::new("a", "Steel", vec![1]));
        assert!(write_error(&mesh).contains("section a is defined twice"));
        mesh.sections.pop();

        mesh.node_sets.insert("Fixed".to_string(), vec![1]);
        mesh.node_sets.insert("FIXED".to_string(), vec![2]);
        assert!(write_error(&mesh).contains("up to case"));
    }

    #[test]
    fn test_element_assigned_to_two_sections_is_rejected() {
        let mut mesh = sample_mesh();
        mesh.nodes.push(AbaqusNode::new(5, [1.0, 1.0, 1.0]));
        mesh.elements
            .push(AbaqusElement::new(2, ElementType::T3D2, vec![4, 5]));
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel", 210e9, 0.3, 7800.0));
        mesh.sections
            .push(AbaqusSection::new("FIRST", "Steel", vec![1, 2, 2]));
        // Repeating an id within one section is harmless.
        assert!(AbaqusWriter::new().write_string(&mesh).is_ok());
        mesh.sections
            .push(AbaqusSection::new("SECOND", "Steel", vec![2]));
        let message = write_error(&mesh);
        assert!(
            message.contains("element 2 is assigned to the sections FIRST and SECOND"),
            "{message}"
        );
    }

    #[test]
    fn test_zero_ids_and_non_finite_numbers_are_rejected() {
        let mut mesh = sample_mesh();
        mesh.nodes.push(AbaqusNode::new(0, [0.0; 3]));
        assert!(write_error(&mesh).contains("ids start at 1"));
        mesh.nodes.pop();

        mesh.nodes[0].coordinates[1] = f64::NAN;
        assert!(write_error(&mesh).contains("non-finite coordinate"));
        mesh.nodes[0].coordinates[1] = 0.0;

        mesh.elements[0].id = 0;
        assert!(write_error(&mesh).contains("ids start at 1"));
        mesh.elements[0].id = 1;

        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut with_material = mesh.clone();
            with_material
                .materials
                .push(AbaqusMaterial::new_elastic("Bad", bad, 0.3, 1.0));
            assert!(write_error(&with_material).contains("non-finite property"));
            let mut with_plastic = mesh.clone();
            with_plastic
                .materials
                .push(AbaqusMaterial::new_plastic("Bad", 1.0, 0.3, 1.0, 1.0, bad));
            assert!(write_error(&with_plastic).contains("non-finite property"));
        }
    }

    #[test]
    fn test_keyword_parameters_survive_characters_that_change_length_when_uppercased() {
        // U+0149 uppercases to two characters (three bytes instead of two) and
        // U+FB01 to two ASCII letters (two bytes instead of three).
        for decoration in ["\u{149}", "\u{fb01}", "\u{149}\u{fb01}\u{149}"] {
            let line = format!("*Element, elset={decoration}, type=c3d8");
            assert_eq!(
                AbaqusReader::extract_param(&line, "TYPE"),
                Some("c3d8".to_string()),
                "{line}"
            );
        }
        // A keyword line that ends right after `type=` used to slice past the
        // end of the line.
        assert_eq!(
            AbaqusReader::extract_param("*Element, \u{149} type=", "TYPE"),
            Some(String::new())
        );
        assert_eq!(
            AbaqusReader::extract_param("*Element, \u{149}", "TYPE"),
            None
        );
    }

    #[test]
    fn test_reader_takes_the_element_type_from_a_keyword_line_with_non_ascii_text() {
        let path = std::env::temp_dir().join("oxiphysics_abaqus_non_ascii_keyword.inp");
        let text = "*Node\n1, 0, 0, 0\n2, 1, 0, 0\n*Element, elset=\u{149}\u{149}\u{149}\u{149}, type=T3D2\n1, 1, 2\n";
        std::fs::write(&path, text).unwrap();
        let parsed = AbaqusReader::new()
            .parse(path.to_str().unwrap_or(""))
            .expect("parse");
        assert_eq!(parsed.elements.len(), 1);
        assert_eq!(parsed.elements[0].element_type, ElementType::T3D2);
        assert_eq!(parsed.elements[0].node_ids, vec![1, 2]);
    }

    #[test]
    fn test_names_that_would_break_the_deck_are_rejected() {
        let mut mesh = sample_mesh();
        mesh.node_sets.insert("A,B".to_string(), vec![1]);
        assert!(write_error(&mesh).contains("node set"));
        mesh.node_sets.clear();
        mesh.materials
            .push(AbaqusMaterial::new_elastic("Steel\n*Step", 1.0, 0.3, 1.0));
        assert!(write_error(&mesh).contains("material"));
    }

    #[test]
    fn test_element_and_node_set_references_to_missing_nodes_are_rejected() {
        let mut mesh = sample_mesh();
        mesh.elements
            .push(AbaqusElement::new(2, ElementType::T3D2, vec![4, 42]));
        assert!(write_error(&mesh).contains("node 42"));
        mesh.elements.pop();
        mesh.node_sets.insert("BAD".to_string(), vec![1, 77]);
        assert!(write_error(&mesh).contains("node 77"));
        mesh.node_sets.insert("BAD".to_string(), Vec::new());
        assert!(write_error(&mesh).contains("empty"));
    }
}
