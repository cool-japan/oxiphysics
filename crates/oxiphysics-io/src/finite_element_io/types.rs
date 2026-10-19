//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;
use std::fmt::Write as _;

use super::functions::DofIndex;

/// Element type mapping between software conventions.
#[derive(Debug, Clone)]
pub struct ElementTypeMapping {
    /// Internal element type.
    pub fe_type: FeElementType,
    /// Abaqus keyword.
    pub abaqus: &'static str,
    /// NASTRAN card keyword.
    pub nastran: &'static str,
    /// Gmsh type tag.
    pub gmsh_tag: usize,
    /// VTK cell type code.
    pub vtk_code: u8,
}
/// A complete finite-element mesh with materials, BCs, and analysis steps.
#[derive(Debug, Clone, Default)]
pub struct FeMesh {
    /// All nodes in the mesh.
    pub nodes: Vec<FeNode>,
    /// All elements.
    pub elements: Vec<FeElement>,
    /// Material definitions.
    pub materials: Vec<LinearElasticMaterial>,
    /// Analysis steps.
    pub steps: Vec<AnalysisStep>,
    /// Model-level boundary conditions: `*BOUNDARY` data given outside any
    /// `*STEP` (before the first step).
    ///
    /// They hold in every step.  The Abaqus reader ([`AbaqusInpIo`]) and the
    /// CalculiX reader and writer ([`crate::calculix_format`]) keep them here,
    /// so a file without an explicit step still yields its supports and a file
    /// with several steps does not repeat them per step.  Code that wants every
    /// condition acting in step `i` must therefore consider both this vector
    /// and `steps[i].bcs`; other readers may instead attach pre-step conditions
    /// to a default step.
    pub global_bcs: Vec<DirichletBc>,
    /// Node set lookup: set name → list of node IDs.
    pub node_sets: HashMap<String, Vec<usize>>,
    /// Element set lookup: set name → list of element IDs.
    pub element_sets: HashMap<String, Vec<usize>>,
}
impl FeMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self::default()
    }
    /// Number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }
    /// Number of elements.
    pub fn num_elements(&self) -> usize {
        self.elements.len()
    }
    /// Find a node by ID (linear scan).
    pub fn find_node(&self, id: usize) -> Option<&FeNode> {
        self.nodes.iter().find(|n| n.id == id)
    }
    /// Find an element by ID (linear scan).
    pub fn find_element(&self, id: usize) -> Option<&FeElement> {
        self.elements.iter().find(|e| e.id == id)
    }
    /// Compute axis-aligned bounding box `([xmin,ymin,zmin], [xmax,ymax,zmax])`.
    /// Returns `None` if the mesh has no nodes.
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for n in &self.nodes {
            for k in 0..3 {
                if n.coords[k] < lo[k] {
                    lo[k] = n.coords[k];
                }
                if n.coords[k] > hi[k] {
                    hi[k] = n.coords[k];
                }
            }
        }
        Some((lo, hi))
    }
}
/// ANSYS MAPDL input deck: what [`MapdlDeck::parse`] reads of an APDL input
/// file or `CDWRITE` CDB file.
#[derive(Debug, Clone, Default)]
pub struct MapdlDeck {
    /// Nodes parsed from `NBLOCK` sections.
    pub nodes: Vec<FeNode>,
    /// Elements parsed from `EBLOCK` sections.
    pub elements: Vec<FeElement>,
    /// Material property table: mat_id → (E, nu, density).
    pub materials: HashMap<usize, (f64, f64, f64)>,
    /// Displacement constraints: node_id → \[(dof, value)\].
    pub disp_constraints: HashMap<usize, Vec<(DofIndex, f64)>>,
    /// Nodal forces: node_id → \[(dof, magnitude)\].
    pub nodal_forces: HashMap<usize, Vec<(DofIndex, f64)>>,
}
impl MapdlDeck {
    /// Create an empty MAPDL deck.
    pub fn new() -> Self {
        Self::default()
    }
    /// Parse the text of an ANSYS MAPDL input deck or `CDWRITE` CDB file.
    ///
    /// Reads the fixed-format `NBLOCK` and `EBLOCK` (`SOLID` layout) blocks,
    /// the `ET` element-type table and the commands `MP,EX`, `MP,NUXY`,
    /// `MP,DENS`, `D` and `F`.
    ///
    /// Each block is cut into fields by the Fortran format line that follows
    /// its command, for example `(3i9,6e21.13e3)`, so values that abut because
    /// a negative number fills its whole field are read correctly.  An
    /// `NBLOCK` row holds `node, entity, line location, x, y, z`; fields that a
    /// short row leaves out are zero.  An `EBLOCK` row holds eleven attribute
    /// fields (material, element type, ..., node count, ..., element number)
    /// followed by the node numbers, which continue on further rows when the
    /// element has more than eight nodes.  The element type follows from the
    /// deck's `ET` table and the node count: solids with 4, 8, 10 or 20 nodes
    /// become `Tet4`, `Hex8`, `Tet10` and `Hex20`, shells and beams of the
    /// well-known ANSYS families become surface and line elements, and any
    /// other element is kept as [`FeElementType::Unknown`] rather than being
    /// mistaken for a solid.
    ///
    /// Blocks and rows that cannot be interpreted are skipped and the rest of
    /// the deck is still read; use [`MapdlDeck::try_parse`] to be told about
    /// them.
    pub fn parse(source: &str) -> Self {
        Self::parse_collecting_issues(source).0
    }
    /// Like [`MapdlDeck::parse`], but fails instead of skipping what it
    /// cannot interpret.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::Parse`], naming the line, for a block whose
    /// format line is missing or unsupported, for an `EBLOCK` that is not in
    /// the `SOLID` layout, for rows with invalid numbers and for an element
    /// whose node list is cut off.
    pub fn try_parse(source: &str) -> crate::Result<Self> {
        let (deck, issues) = Self::parse_collecting_issues(source);
        if issues.is_empty() {
            Ok(deck)
        } else {
            Err(crate::Error::Parse(format!(
                "MAPDL deck: {}",
                issues.join("; ")
            )))
        }
    }
    fn parse_collecting_issues(source: &str) -> (Self, Vec<String>) {
        let lines: Vec<&str> = source.lines().collect();
        let mut deck = Self::new();
        let mut issues = Vec::new();
        let mut element_types = super::mapdl_cdb::ElementTypeTable::new();
        let mut at = 0;
        while let Some(&line) = lines.get(at) {
            at += 1;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('!') {
                continue;
            }
            let upper = trimmed.to_ascii_uppercase();
            if super::mapdl_cdb::is_command(&upper, "NBLOCK") {
                at = super::mapdl_cdb::read_nblock(&lines, at, &mut deck.nodes, &mut issues);
            } else if super::mapdl_cdb::is_command(&upper, "EBLOCK") {
                at = super::mapdl_cdb::read_eblock(
                    &lines,
                    at,
                    trimmed,
                    &element_types,
                    &mut deck.elements,
                    &mut issues,
                );
            } else if upper.starts_with("ET,") {
                super::mapdl_cdb::record_element_type(trimmed, &mut element_types);
            } else if upper.starts_with("MP,") {
                Self::apply_mp_command(trimmed, &mut deck);
            } else if upper.starts_with("D,") {
                Self::apply_d_command(trimmed, &mut deck);
            } else if upper.starts_with("F,") {
                Self::apply_f_command(trimmed, &mut deck);
            }
        }
        (deck, issues)
    }
    fn apply_mp_command(line: &str, deck: &mut MapdlDeck) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            return;
        }
        let prop = parts[1].trim().to_uppercase();
        let mat_id: usize = match parts[2].trim().parse() {
            Ok(v) => v,
            Err(_) => return,
        };
        let value: f64 = match parts[3].trim().parse() {
            Ok(v) => v,
            Err(_) => return,
        };
        let entry = deck.materials.entry(mat_id).or_insert((0.0, 0.3, 0.0));
        match prop.as_str() {
            "EX" => entry.0 = value,
            "NUXY" => entry.1 = value,
            "DENS" => entry.2 = value,
            _ => {}
        }
    }
    fn apply_d_command(line: &str, deck: &mut MapdlDeck) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 3 {
            return;
        }
        let node_id: usize = match parts[1].trim().parse() {
            Ok(v) => v,
            Err(_) => return,
        };
        let dof_str = parts[2].trim().to_uppercase();
        let dof: DofIndex = match dof_str.as_str() {
            "UX" => 1,
            "UY" => 2,
            "UZ" => 3,
            "ROTX" => 4,
            "ROTY" => 5,
            "ROTZ" => 6,
            _ => return,
        };
        let value: f64 = parts
            .get(3)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0);
        deck.disp_constraints
            .entry(node_id)
            .or_default()
            .push((dof, value));
    }
    fn apply_f_command(line: &str, deck: &mut MapdlDeck) {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            return;
        }
        let node_id: usize = match parts[1].trim().parse() {
            Ok(v) => v,
            Err(_) => return,
        };
        let dof_str = parts[2].trim().to_uppercase();
        let dof: DofIndex = match dof_str.as_str() {
            "FX" => 1,
            "FY" => 2,
            "FZ" => 3,
            "MX" => 4,
            "MY" => 5,
            "MZ" => 6,
            _ => return,
        };
        let mag: f64 = match parts[3].trim().parse() {
            Ok(v) => v,
            Err(_) => return,
        };
        deck.nodal_forces
            .entry(node_id)
            .or_default()
            .push((dof, mag));
    }
}
/// A Dirichlet (displacement) boundary condition on a node DOF.
#[derive(Debug, Clone, PartialEq)]
pub struct DirichletBc {
    /// Node ID.
    pub node_id: usize,
    /// Degree of freedom (1–6).
    pub dof: DofIndex,
    /// Prescribed value (0.0 for a fixed support).
    pub value: f64,
}
impl DirichletBc {
    /// Create a fixed-support condition (value = 0).
    pub fn fixed(node_id: usize, dof: DofIndex) -> Self {
        Self {
            node_id,
            dof,
            value: 0.0,
        }
    }
    /// Create a prescribed displacement condition.
    pub fn prescribed(node_id: usize, dof: DofIndex, value: f64) -> Self {
        Self {
            node_id,
            dof,
            value,
        }
    }
}
/// FORCE card.
#[derive(Debug, Clone, PartialEq)]
pub struct NastranForce {
    /// Set ID.
    pub sid: usize,
    /// Grid point ID.
    pub gid: usize,
    /// Scale factor.
    pub scale: f64,
    /// Force direction vector `[nx, ny, nz]`.
    pub direction: [f64; 3],
}
/// A parsed Gmsh mesh.
#[derive(Debug, Clone, Default)]
pub struct GmshMesh {
    /// Format version detected.
    pub version: Option<GmshVersion>,
    /// Nodes (ID, coordinates).
    pub nodes: Vec<FeNode>,
    /// Elements.
    pub elements: Vec<FeElement>,
    /// Physical group names: (dim, tag) → name.
    pub physical_names: HashMap<(i32, i32), String>,
}
impl GmshMesh {
    /// Create an empty Gmsh mesh.
    pub fn new() -> Self {
        Self::default()
    }
    /// Parse a Gmsh `.msh` ASCII file in format 2.x or 4.1.
    ///
    /// Nodes and elements come from `$Nodes` and `$Elements`, physical group
    /// names from `$PhysicalNames`; every other section is skipped.  A 4.1
    /// file is read through its entity blocks (`$Nodes`: `entityDim entityTag
    /// parametric numNodes`, then the node tags, then the coordinates;
    /// `$Elements`: `entityDim entityTag elementType numElements`, then
    /// `elementTag nodeTag...`).  Gmsh element type numbers without a
    /// [`FeElementType`] are kept as `Unknown("gmsh:<number>")`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::Parse`], naming the line, for binary files
    /// (file type 1), for format versions other than 2.x and 4.1, for a
    /// missing `$MeshFormat`, for sections that are cut off or whose counts
    /// disagree with their content, for malformed numbers, and for element
    /// lines that declare more tags than they hold or list the wrong number
    /// of nodes for their type.
    pub fn parse(source: &str) -> Result<Self, crate::Error> {
        let mut mesh = Self::new();
        super::gmsh_msh::parse_into(&mut mesh, source)?;
        Ok(mesh)
    }
    /// Write the mesh in Gmsh v2 ASCII format.
    pub fn write_v2(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "$MeshFormat\n2.2 0 8\n$EndMeshFormat");
        let _ = writeln!(out, "$Nodes\n{}", self.nodes.len());
        for n in &self.nodes {
            let _ = writeln!(
                out,
                "{} {:.10e} {:.10e} {:.10e}",
                n.id, n.coords[0], n.coords[1], n.coords[2]
            );
        }
        let _ = writeln!(out, "$EndNodes");
        let _ = writeln!(out, "$Elements\n{}", self.elements.len());
        for e in &self.elements {
            let type_tag = Self::fe_to_gmsh(&e.element_type);
            let ids: Vec<String> = e.connectivity.iter().map(|n| n.to_string()).collect();
            let _ = writeln!(out, "{} {} 0 {}", e.id, type_tag, ids.join(" "));
        }
        let _ = writeln!(out, "$EndElements");
        out
    }
    fn fe_to_gmsh(t: &FeElementType) -> usize {
        super::gmsh_msh::tag_of_element_type(t)
    }
}
/// A load/analysis step.
#[derive(Debug, Clone)]
pub struct AnalysisStep {
    /// Step name.
    pub name: String,
    /// Boundary conditions active in this step.
    pub bcs: Vec<DirichletBc>,
    /// Loads applied in this step.
    pub forces: Vec<NodalForce>,
    /// Step time period (for transient analyses).
    pub time_period: f64,
    /// Suggested initial increment size.
    pub initial_increment: f64,
}
impl AnalysisStep {
    /// Create a new step with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            bcs: Vec::new(),
            forces: Vec::new(),
            time_period: 1.0,
            initial_increment: 0.01,
        }
    }
}
/// Parse state for the Abaqus .inp reader.
#[derive(Debug)]
pub(super) struct AbaqusParseState {
    pub(super) node_section: bool,
    pub(super) element_section: bool,
    pub(super) current_elem_type: FeElementType,
    pub(super) node_set_section: Option<String>,
    pub(super) elem_set_section: Option<String>,
    pub(super) material_name: Option<String>,
    pub(super) in_elastic: bool,
    /// Inside the data of a `*BOUNDARY` keyword that prescribes displacements.
    pub(super) in_boundary: bool,
    /// Inside the data of a `*CLOAD` keyword.
    pub(super) in_cload: bool,
    /// The current `*NSET` / `*ELSET` keyword has the `GENERATE` parameter.
    pub(super) generate: bool,
    /// Between `*STEP` and `*END STEP`.
    pub(super) in_step: bool,
    /// Set named by `NSET=` on the current `*NODE` keyword.
    pub(super) node_collect_set: Option<String>,
    /// Set named by `ELSET=` on the current `*ELEMENT` keyword.
    pub(super) element_collect_set: Option<String>,
    /// Element record that continues on the next data line, with the number
    /// of the line it started on.
    pub(super) pending_element: Option<(usize, String)>,
}
/// Isotropic linear-elastic material.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearElasticMaterial {
    /// Unique material identifier.
    pub id: usize,
    /// Material name.
    pub name: String,
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Density \[kg/m³\].
    pub density: f64,
}
impl LinearElasticMaterial {
    /// Create a new linear-elastic material.
    pub fn new(
        id: usize,
        name: impl Into<String>,
        young_modulus: f64,
        poisson_ratio: f64,
        density: f64,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            young_modulus,
            poisson_ratio,
            density,
        }
    }
    /// Shear modulus derived from `E` and `ν`.
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
    /// Bulk modulus derived from `E` and `ν`.
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
}
/// NASTRAN element (CTETRA or CHEXA).
#[derive(Debug, Clone, PartialEq)]
pub struct NastranElement {
    /// Element ID (EID).
    pub eid: usize,
    /// Property ID (PID).
    pub pid: usize,
    /// Element type.
    pub element_type: FeElementType,
    /// Node connectivity.
    pub nodes: Vec<usize>,
}
/// Container for a parsed NASTRAN Bulk Data section.
#[derive(Debug, Clone, Default)]
pub struct NastranBulkData {
    /// Grid points.
    pub grids: Vec<NastranGrid>,
    /// Elements.
    pub elements: Vec<NastranElement>,
    /// MAT1 materials.
    pub materials: Vec<NastranMat1>,
    /// SPC constraints.
    pub spcs: Vec<NastranSpc>,
    /// FORCE loads.
    pub forces: Vec<NastranForce>,
}
impl NastranBulkData {
    /// Create empty bulk data.
    pub fn new() -> Self {
        Self::default()
    }
    /// Parse NASTRAN Bulk Data from a string.
    ///
    /// Only the section between `BEGIN BULK` and `ENDDATA` is read.  Small
    /// field (8-column), large field (`CARD*`, 16-column) and free-field
    /// (comma separated) cards are understood, and cards that continue on
    /// further lines (a marker such as `+C1` or `*K1`, a blank name field, or
    /// a leading comma) are joined before they are interpreted.
    ///
    /// The cards `GRID`, `CTETRA` (4 or 10 nodes), `CHEXA` (8 or 20 nodes),
    /// `MAT1`, `SPC` (one or two constraints per card) and `FORCE` are read;
    /// a card is recognised by its exact name, so `FORCE1` and `SPC1` are not
    /// mistaken for `FORCE` and `SPC`.
    ///
    /// Real numbers may be written the Nastran way without the exponent
    /// letter (`2.1+11` is `2.1E+11`).  Blank fields take the defaults of the
    /// card definition (`0.0` for grid coordinates and for the direction of a
    /// `FORCE`, `0` for coordinate system ids); the elastic constants of a
    /// `MAT1` card that are left blank are derived from the other two, with
    /// Poisson's ratio taken as `0.3` when only one constant is given.  A card
    /// with a field that is not a number, or one that is malformed in another
    /// way, such as an element with a missing or surplus node, is skipped
    /// instead of being read with invented values.
    pub fn parse(source: &str) -> Self {
        let mut joiner = super::nastran_bulk::CardJoiner::new();
        let mut in_bulk = false;
        for raw_line in source.lines() {
            let line = super::nastran_bulk::strip_comment(raw_line);
            if line.trim().is_empty() {
                continue;
            }
            match super::nastran_bulk::section_marker(line) {
                Some(super::nastran_bulk::SectionMarker::BeginBulk) => in_bulk = true,
                Some(super::nastran_bulk::SectionMarker::EndData) => break,
                None if in_bulk => joiner.push_line(line),
                None => {}
            }
        }
        let mut data = Self::new();
        for card in joiner.finish() {
            data.add_card(&card);
        }
        data
    }
    /// Interpret one joined card (`fields[0]` is the card name).
    fn add_card(&mut self, fields: &[String]) {
        let name = fields
            .first()
            .map(|name| name.trim_end_matches('*').to_ascii_uppercase())
            .unwrap_or_default();
        match name.as_str() {
            "GRID" => {
                if let Some(g) = Self::parse_grid(fields) {
                    self.grids.push(g);
                }
            }
            "CTETRA" => {
                if let Some(e) = Self::parse_ctetra(fields) {
                    self.elements.push(e);
                }
            }
            "CHEXA" => {
                if let Some(e) = Self::parse_chexa(fields) {
                    self.elements.push(e);
                }
            }
            "MAT1" => {
                if let Some(m) = Self::parse_mat1(fields) {
                    self.materials.push(m);
                }
            }
            "SPC" => self
                .spcs
                .extend(super::nastran_bulk::spc_constraints(fields)),
            "FORCE" => {
                if let Some(f) = Self::parse_force(fields) {
                    self.forces.push(f);
                }
            }
            _ => {}
        }
    }
    fn parse_grid(f: &[String]) -> Option<NastranGrid> {
        use super::nastran_bulk::{optional_id, optional_real, required_id};
        let id = required_id(f, 1)?;
        let cp = optional_id(f, 2, 0)?;
        let x = optional_real(f, 3, 0.0)?;
        let y = optional_real(f, 4, 0.0)?;
        let z = optional_real(f, 5, 0.0)?;
        let cd = optional_id(f, 6, 0)?;
        Some(NastranGrid {
            id,
            cp,
            coords: [x, y, z],
            cd,
        })
    }
    fn parse_ctetra(f: &[String]) -> Option<NastranElement> {
        let eid: usize = f.get(1)?.parse().ok()?;
        let pid = super::nastran_bulk::property_id(f, eid)?;
        let nodes = super::nastran_bulk::element_node_ids(f)?;
        let element_type = match nodes.len() {
            4 => FeElementType::Tet4,
            10 => FeElementType::Tet10,
            _ => return None,
        };
        Some(NastranElement {
            eid,
            pid,
            element_type,
            nodes,
        })
    }
    fn parse_chexa(f: &[String]) -> Option<NastranElement> {
        let eid: usize = f.get(1)?.parse().ok()?;
        let pid = super::nastran_bulk::property_id(f, eid)?;
        let nodes = super::nastran_bulk::element_node_ids(f)?;
        let element_type = match nodes.len() {
            8 => FeElementType::Hex8,
            20 => FeElementType::Hex20,
            _ => return None,
        };
        Some(NastranElement {
            eid,
            pid,
            element_type,
            nodes,
        })
    }
    fn parse_mat1(f: &[String]) -> Option<NastranMat1> {
        use super::nastran_bulk::{blank_or_real, complete_moduli, optional_real, required_id};
        let mid = required_id(f, 1)?;
        let (e, g, nu) = complete_moduli(
            blank_or_real(f, 2)?,
            blank_or_real(f, 3)?,
            blank_or_real(f, 4)?,
        )?;
        let rho = optional_real(f, 5, 0.0)?;
        Some(NastranMat1 { mid, e, g, nu, rho })
    }
    fn parse_force(f: &[String]) -> Option<NastranForce> {
        use super::nastran_bulk::{optional_id, optional_real, required_id, required_real};
        let sid = required_id(f, 1)?;
        let gid = required_id(f, 2)?;
        // The coordinate system id is checked but not used: forces are kept
        // as written, so a direction given in a local system is not rotated.
        optional_id(f, 3, 0)?;
        let scale = required_real(f, 4)?;
        let nx = optional_real(f, 5, 0.0)?;
        let ny = optional_real(f, 6, 0.0)?;
        let nz = optional_real(f, 7, 0.0)?;
        Some(NastranForce {
            sid,
            gid,
            scale,
            direction: [nx, ny, nz],
        })
    }
    /// Write bulk data in free-field comma-separated format.
    ///
    /// Elements are written as `CTETRA` / `CHEXA` cards whose grid points beyond
    /// the sixth continue on `+` lines, so `parse(write(data))` returns the same
    /// elements.  An element of a type that has no card here (or whose node count
    /// does not fit its type) is written as a `$ unsupported element` comment.
    pub fn write(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "$ Generated by OxiPhysics finite_element_io");
        let _ = writeln!(out, "BEGIN BULK");
        for g in &self.grids {
            let _ = writeln!(
                out,
                "GRID,{},{},{:.10e},{:.10e},{:.10e},{}",
                g.id, g.cp, g.coords[0], g.coords[1], g.coords[2], g.cd
            );
        }
        let mut markers = 0usize;
        for element in &self.elements {
            super::nastran_bulk::write_element(&mut out, element, &mut markers);
        }
        for m in &self.materials {
            let _ = writeln!(
                out,
                "MAT1,{},{:.6e},{:.6e},{:.4},{:.4}",
                m.mid, m.e, m.g, m.nu, m.rho
            );
        }
        for s in &self.spcs {
            let _ = writeln!(
                out,
                "SPC,{},{},{},{:.6e}",
                s.sid, s.gid, s.components, s.value
            );
        }
        for f in &self.forces {
            let _ = writeln!(
                out,
                "FORCE,{},{},0,{:.6e},{:.6},{:.6},{:.6}",
                f.sid, f.gid, f.scale, f.direction[0], f.direction[1], f.direction[2]
            );
        }
        let _ = writeln!(out, "ENDDATA");
        out
    }
}
/// A simplified in-memory ExodusII-like mesh.
#[derive(Debug, Clone, Default)]
pub struct ExodusLikeMesh {
    /// Title string.
    pub title: String,
    /// Node coordinates (flat: \[x0,y0,z0, x1,y1,z1, ...\]).
    pub coordinates: Vec<f64>,
    /// Number of nodes.
    pub num_nodes: usize,
    /// Element blocks: block_id → (element_type, connectivity flat list).
    pub element_blocks: HashMap<usize, (FeElementType, Vec<usize>)>,
    /// Node sets: set_id → \[node_ids\].
    pub node_sets: HashMap<usize, Vec<usize>>,
    /// Side sets: set_id → \[(element_id, face_ordinal)\].
    pub side_sets: HashMap<usize, Vec<(usize, usize)>>,
    /// Nodal variable names.
    pub nodal_variable_names: Vec<String>,
    /// Nodal variable data: var_index → \[values per node\].
    pub nodal_variables: HashMap<usize, Vec<f64>>,
}
impl ExodusLikeMesh {
    /// Create an empty mesh.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            ..Self::default()
        }
    }
    /// Access x-coordinate of node `i` (0-based).
    pub fn x(&self, i: usize) -> f64 {
        self.coordinates[3 * i]
    }
    /// Access y-coordinate of node `i` (0-based).
    pub fn y(&self, i: usize) -> f64 {
        self.coordinates[3 * i + 1]
    }
    /// Access z-coordinate of node `i` (0-based).
    pub fn z(&self, i: usize) -> f64 {
        self.coordinates[3 * i + 2]
    }
    /// Total element count across all blocks.
    pub fn total_elements(&self) -> usize {
        self.element_blocks
            .values()
            .map(|(et, conn)| {
                let npe = et.node_count().unwrap_or(1);
                conn.len() / npe.max(1)
            })
            .sum()
    }
    /// Serialize to a simplified text representation (not real ExodusII binary).
    pub fn serialize_text(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "TITLE: {}", self.title);
        let _ = writeln!(out, "NUM_NODES: {}", self.num_nodes);
        let _ = writeln!(out, "COORDINATES:");
        for i in 0..self.num_nodes {
            let _ = writeln!(
                out,
                "  {} {:.10e} {:.10e} {:.10e}",
                i + 1,
                self.x(i),
                self.y(i),
                self.z(i)
            );
        }
        let _ = writeln!(out, "ELEMENT_BLOCKS: {}", self.element_blocks.len());
        for (blk_id, (et, conn)) in &self.element_blocks {
            let npe = et.node_count().unwrap_or(1);
            let ne = conn.len() / npe.max(1);
            let _ = writeln!(out, "  BLOCK {} {:?} ne={} npe={}", blk_id, et, ne, npe);
        }
        out
    }
}
/// MAT1 (isotropic material) card data.
#[derive(Debug, Clone, PartialEq)]
pub struct NastranMat1 {
    /// Material ID.
    pub mid: usize,
    /// Young's modulus.
    pub e: f64,
    /// Shear modulus (may be derived from E and nu).
    pub g: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Density.
    pub rho: f64,
}
/// A single NASTRAN GRID card.
#[derive(Debug, Clone, PartialEq)]
pub struct NastranGrid {
    /// Grid point ID.
    pub id: usize,
    /// Coordinate system ID (0 = basic).
    pub cp: usize,
    /// Coordinates `[x, y, z]`.
    pub coords: [f64; 3],
    /// Analysis coordinate system ID.
    pub cd: usize,
}
impl NastranGrid {
    /// Create a grid point in the basic coordinate system.
    pub fn new(id: usize, coords: [f64; 3]) -> Self {
        Self {
            id,
            cp: 0,
            coords,
            cd: 0,
        }
    }
}
/// Multi-point constraint (MPC) definition.
#[derive(Debug, Clone)]
pub struct MultiPointConstraint {
    /// Dependent node ID.
    pub dependent_node: usize,
    /// Dependent DOF.
    pub dependent_dof: DofIndex,
    /// List of (independent_node_id, independent_dof, coefficient).
    pub terms: Vec<(usize, DofIndex, f64)>,
    /// Constant term on the right-hand side.
    pub rhs: f64,
}
impl MultiPointConstraint {
    /// Create a new MPC.
    pub fn new(dependent_node: usize, dependent_dof: DofIndex) -> Self {
        Self {
            dependent_node,
            dependent_dof,
            terms: Vec::new(),
            rhs: 0.0,
        }
    }
    /// Add an independent term.
    pub fn add_term(&mut self, node: usize, dof: DofIndex, coeff: f64) {
        self.terms.push((node, dof, coeff));
    }
}
/// A concentrated nodal force.
#[derive(Debug, Clone, PartialEq)]
pub struct NodalForce {
    /// Node ID.
    pub node_id: usize,
    /// Degree of freedom (1–6).
    pub dof: DofIndex,
    /// Force magnitude \[N\] (or moment \[N·m\] for rotational DOFs).
    pub magnitude: f64,
}
impl NodalForce {
    /// Create a nodal force.
    pub fn new(node_id: usize, dof: DofIndex, magnitude: f64) -> Self {
        Self {
            node_id,
            dof,
            magnitude,
        }
    }
}
/// Nodal stress data (Voigt notation: σxx, σyy, σzz, σxy, σyz, σxz).
#[derive(Debug, Clone)]
pub struct NodalStresses {
    /// Node IDs.
    pub node_ids: Vec<usize>,
    /// Stress tensors in Voigt notation per node.
    pub stresses: Vec<[f64; 6]>,
}
impl NodalStresses {
    /// Create a zero-stress field for a list of node IDs.
    pub fn zeros(node_ids: Vec<usize>) -> Self {
        let n = node_ids.len();
        Self {
            node_ids,
            stresses: vec![[0.0; 6]; n],
        }
    }
    /// Von Mises stress at node `i`.
    pub fn von_mises(&self, i: usize) -> f64 {
        let s = self.stresses[i];
        let ds = [s[0] - s[1], s[1] - s[2], s[2] - s[0]];
        let shear_sq = s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
        let vm_sq = 0.5 * (ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2]) + 3.0 * shear_sq;
        vm_sq.sqrt()
    }
}
/// Abaqus `.inp` file reader and writer.
#[derive(Debug, Default)]
pub struct AbaqusInpIo;
impl AbaqusInpIo {
    /// Parse an Abaqus `.inp` file from a string.
    ///
    /// Reads `*NODE`, `*ELEMENT`, `*NSET`, `*ELSET`, `*MATERIAL`, `*ELASTIC`
    /// (isotropic), `*DENSITY`, `*BOUNDARY`, `*CLOAD`, `*STEP` and
    /// `*END STEP`; other keywords are skipped.  Nodes may omit trailing coordinates (two-dimensional models),
    /// elements may continue over several lines, `*NSET` / `*ELSET` accept
    /// `GENERATE` ranges, and `NSET=` / `ELSET=` on `*NODE` / `*ELEMENT`
    /// collect sets.  `*BOUNDARY` (`node|set, first[, last[, value]]` or a
    /// named constraint such as `ENCASTRE`) and `*CLOAD`
    /// (`node|set, dof, magnitude`) become [`DirichletBc`]s and
    /// [`NodalForce`]s, one per node and degree of freedom: inside a `*STEP`
    /// in that step, in the model data in [`FeMesh::global_bcs`].
    ///
    /// Data lines that cannot be interpreted (malformed numbers, references to
    /// undefined sets, `*CLOAD` outside a step, elements with the wrong number
    /// of nodes) are skipped; [`AbaqusInpIo::try_parse`] reports them instead.
    pub fn parse(source: &str) -> FeMesh {
        super::abaqus_inp::parse_inp(source).0
    }
    /// Like [`AbaqusInpIo::parse`], but fails instead of skipping data lines
    /// that cannot be interpreted.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::Parse`] listing every skipped line, each with
    /// its line number.
    pub fn try_parse(source: &str) -> crate::Result<FeMesh> {
        let (mesh, issues) = super::abaqus_inp::parse_inp(source);
        if issues.is_empty() {
            Ok(mesh)
        } else {
            Err(crate::Error::Parse(format!(
                "Abaqus input: {}",
                issues.join("; ")
            )))
        }
    }
    /// Write a mesh to Abaqus `.inp` format.
    pub fn write(mesh: &FeMesh) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "** Generated by OxiPhysics finite_element_io");
        let _ = writeln!(out, "*NODE");
        for n in &mesh.nodes {
            let _ = writeln!(
                out,
                "{}, {:.10e}, {:.10e}, {:.10e}",
                n.id, n.coords[0], n.coords[1], n.coords[2]
            );
        }
        let mut by_type: HashMap<String, Vec<&FeElement>> = HashMap::new();
        for e in &mesh.elements {
            let key = format!("{:?}", e.element_type);
            by_type.entry(key).or_default().push(e);
        }
        for (type_name, elems) in &by_type {
            let abaqus_type = Self::fe_type_to_abaqus(type_name);
            let _ = writeln!(out, "*ELEMENT, TYPE={}", abaqus_type);
            for e in elems {
                let ids: Vec<String> = e.connectivity.iter().map(|i| i.to_string()).collect();
                let _ = writeln!(out, "{}, {}", e.id, ids.join(", "));
            }
        }
        for mat in &mesh.materials {
            let _ = writeln!(out, "*MATERIAL, NAME={}", mat.name);
            let _ = writeln!(out, "*ELASTIC");
            let _ = writeln!(out, "{:.6e}, {:.6}", mat.young_modulus, mat.poisson_ratio);
        }
        out
    }
    fn fe_type_to_abaqus(type_name: &str) -> &'static str {
        match type_name {
            "Tet4" => "C3D4",
            "Tet10" => "C3D10",
            "Hex8" => "C3D8",
            "Hex20" => "C3D20",
            "Tri3" => "S3",
            "Quad4" => "S4",
            "Line2" => "T3D2",
            _ => "UNKNOWN",
        }
    }
}
/// SPC (single-point constraint) card.
#[derive(Debug, Clone, PartialEq)]
pub struct NastranSpc {
    /// SPC set ID.
    pub sid: usize,
    /// Grid point ID.
    pub gid: usize,
    /// DOF component string (e.g. "123456").
    pub components: String,
    /// Displacement value.
    pub value: f64,
}
/// A 3-D node: an integer ID and spatial coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct FeNode {
    /// 1-based node identifier.
    pub id: usize,
    /// Cartesian coordinates `[x, y, z]`.
    pub coords: [f64; 3],
}
impl FeNode {
    /// Create a new node.
    pub fn new(id: usize, coords: [f64; 3]) -> Self {
        Self { id, coords }
    }
    /// Euclidean distance to another node.
    pub fn distance(&self, other: &FeNode) -> f64 {
        let dx = self.coords[0] - other.coords[0];
        let dy = self.coords[1] - other.coords[1];
        let dz = self.coords[2] - other.coords[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
/// Supported finite-element types across all formats.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FeElementType {
    /// 4-node linear tetrahedron.
    Tet4,
    /// 10-node quadratic tetrahedron.
    Tet10,
    /// 8-node linear hexahedron.
    Hex8,
    /// 20-node quadratic hexahedron.
    Hex20,
    /// 3-node linear triangle (shell / 2-D).
    Tri3,
    /// 6-node quadratic triangle.
    Tri6,
    /// 4-node linear quadrilateral.
    Quad4,
    /// 8-node quadratic quadrilateral.
    Quad8,
    /// 2-node line (truss / beam).
    Line2,
    /// Unknown element type (carries the raw type string).
    Unknown(String),
}
impl FeElementType {
    /// Number of nodes per element for built-in types; `None` for `Unknown`.
    pub fn node_count(&self) -> Option<usize> {
        match self {
            FeElementType::Tet4 => Some(4),
            FeElementType::Tet10 => Some(10),
            FeElementType::Hex8 => Some(8),
            FeElementType::Hex20 => Some(20),
            FeElementType::Tri3 => Some(3),
            FeElementType::Tri6 => Some(6),
            FeElementType::Quad4 => Some(4),
            FeElementType::Quad8 => Some(8),
            FeElementType::Line2 => Some(2),
            FeElementType::Unknown(_) => None,
        }
    }
    /// Whether this is a volumetric (3-D) element.
    pub fn is_volumetric(&self) -> bool {
        matches!(
            self,
            FeElementType::Tet4 | FeElementType::Tet10 | FeElementType::Hex8 | FeElementType::Hex20
        )
    }
}
/// Gmsh mesh format version.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GmshVersion {
    /// Version 2 ASCII format.
    V2,
    /// Version 4 ASCII format.
    V4,
}
/// Solver state checkpoint for restart capability.
#[derive(Debug, Clone)]
pub struct RestartCheckpoint {
    /// Simulation title.
    pub title: String,
    /// Step number at the time of writing.
    pub step: usize,
    /// Simulation time at checkpoint.
    pub time: f64,
    /// Nodal displacement solution vector (flat: \[ux1, uy1, uz1, ux2, ...\]).
    pub displacements: Vec<f64>,
    /// Nodal velocity vector (flat).
    pub velocities: Vec<f64>,
    /// Solver residual norm at this checkpoint.
    pub residual_norm: f64,
    /// User-defined metadata key-value pairs.
    pub metadata: HashMap<String, String>,
}
impl RestartCheckpoint {
    /// Create a new checkpoint.
    pub fn new(title: impl Into<String>, step: usize, time: f64) -> Self {
        Self {
            title: title.into(),
            step,
            time,
            displacements: Vec::new(),
            velocities: Vec::new(),
            residual_norm: 0.0,
            metadata: HashMap::new(),
        }
    }
    /// Serialize to a simple text format.
    pub fn serialize(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "RESTART_FILE");
        let _ = writeln!(out, "TITLE: {}", self.title);
        let _ = writeln!(out, "STEP: {}", self.step);
        let _ = writeln!(out, "TIME: {:.15e}", self.time);
        let _ = writeln!(out, "RESIDUAL_NORM: {:.15e}", self.residual_norm);
        let _ = writeln!(out, "METADATA:");
        for (k, v) in &self.metadata {
            let _ = writeln!(out, "  {}={}", k, v);
        }
        let _ = writeln!(out, "DISPLACEMENTS: {}", self.displacements.len());
        for v in &self.displacements {
            let _ = writeln!(out, "  {:.15e}", v);
        }
        let _ = writeln!(out, "VELOCITIES: {}", self.velocities.len());
        for v in &self.velocities {
            let _ = writeln!(out, "  {:.15e}", v);
        }
        out
    }
    /// Deserialize from the text format produced by [`RestartCheckpoint::serialize`].
    pub fn deserialize(text: &str) -> Option<Self> {
        let mut title = String::new();
        let mut step = 0usize;
        let mut time = 0.0f64;
        let mut residual_norm = 0.0f64;
        let mut displacements = Vec::new();
        let mut velocities = Vec::new();
        let mut metadata = HashMap::new();
        let mut mode = "";
        for line in text.lines() {
            let t = line.trim();
            if t.is_empty() {
                continue;
            }
            if t == "RESTART_FILE" {
                continue;
            }
            if let Some(rest) = t.strip_prefix("TITLE:") {
                title = rest.trim().to_string();
            } else if let Some(rest) = t.strip_prefix("STEP:") {
                step = rest.trim().parse().unwrap_or(0);
            } else if let Some(rest) = t.strip_prefix("TIME:") {
                time = rest.trim().parse().unwrap_or(0.0);
            } else if let Some(rest) = t.strip_prefix("RESIDUAL_NORM:") {
                residual_norm = rest.trim().parse().unwrap_or(0.0);
            } else if t.starts_with("METADATA:") {
                mode = "metadata";
            } else if t.starts_with("DISPLACEMENTS:") {
                mode = "displacements";
            } else if t.starts_with("VELOCITIES:") {
                mode = "velocities";
            } else {
                match mode {
                    "metadata" => {
                        if let Some(eq) = t.find('=') {
                            let k = t[..eq].trim().to_string();
                            let v = t[eq + 1..].trim().to_string();
                            metadata.insert(k, v);
                        }
                    }
                    "displacements" => {
                        if let Ok(v) = t.parse::<f64>() {
                            displacements.push(v);
                        }
                    }
                    "velocities" => {
                        if let Ok(v) = t.parse::<f64>() {
                            velocities.push(v);
                        }
                    }
                    _ => {}
                }
            }
        }
        Some(Self {
            title,
            step,
            time,
            displacements,
            velocities,
            residual_norm,
            metadata,
        })
    }
}
/// A single finite element.
#[derive(Debug, Clone, PartialEq)]
pub struct FeElement {
    /// 1-based element identifier.
    pub id: usize,
    /// Element topology.
    pub element_type: FeElementType,
    /// Ordered list of node IDs that form this element.
    pub connectivity: Vec<usize>,
    /// Optional material/section tag.
    pub material_id: Option<usize>,
}
impl FeElement {
    /// Create a new element.
    pub fn new(id: usize, element_type: FeElementType, connectivity: Vec<usize>) -> Self {
        Self {
            id,
            element_type,
            connectivity,
            material_id: None,
        }
    }
    /// Number of nodes in this element's connectivity list.
    pub fn num_nodes(&self) -> usize {
        self.connectivity.len()
    }
}
/// Mesh partitioning metadata for parallel FE simulations.
#[derive(Debug, Clone)]
pub struct MeshPartition {
    /// Partition (rank) index, 0-based.
    pub rank: usize,
    /// Total number of partitions.
    pub num_parts: usize,
    /// Node IDs owned by this partition.
    pub owned_nodes: Vec<usize>,
    /// Node IDs that are ghost (shared from neighbouring partitions).
    pub ghost_nodes: Vec<usize>,
    /// Element IDs owned by this partition.
    pub owned_elements: Vec<usize>,
    /// Neighbouring partition ranks.
    pub neighbours: Vec<usize>,
}
impl MeshPartition {
    /// Create a new partition descriptor.
    pub fn new(rank: usize, num_parts: usize) -> Self {
        Self {
            rank,
            num_parts,
            owned_nodes: Vec::new(),
            ghost_nodes: Vec::new(),
            owned_elements: Vec::new(),
            neighbours: Vec::new(),
        }
    }
    /// Total node count (owned + ghost).
    pub fn total_nodes(&self) -> usize {
        self.owned_nodes.len() + self.ghost_nodes.len()
    }
    /// Fraction of overall partitions that are neighbours.
    pub fn neighbour_fraction(&self) -> f64 {
        if self.num_parts <= 1 {
            0.0
        } else {
            self.neighbours.len() as f64 / (self.num_parts - 1) as f64
        }
    }
}
/// Nodal displacement data for VTK export.
#[derive(Debug, Clone)]
pub struct NodalDisplacements {
    /// Node IDs in the same order as `displacements`.
    pub node_ids: Vec<usize>,
    /// Displacement vectors `[dx, dy, dz]` per node.
    pub displacements: Vec<[f64; 3]>,
}
impl NodalDisplacements {
    /// Create a zero-displacement field for a list of node IDs.
    pub fn zeros(node_ids: Vec<usize>) -> Self {
        let n = node_ids.len();
        Self {
            node_ids,
            displacements: vec![[0.0; 3]; n],
        }
    }
    /// Magnitude of displacement at index `i`.
    pub fn magnitude(&self, i: usize) -> f64 {
        let d = self.displacements[i];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }
    /// Maximum displacement magnitude over all nodes.
    pub fn max_magnitude(&self) -> f64 {
        (0..self.node_ids.len())
            .map(|i| self.magnitude(i))
            .fold(0.0_f64, f64::max)
    }
}
/// VTK result exporter for FE data.
#[derive(Debug, Default)]
pub struct VtkResultExporter;
impl VtkResultExporter {
    /// Export mesh + displacements + stresses to VTK legacy ASCII format.
    pub fn export_legacy_ascii(
        mesh: &FeMesh,
        displacements: Option<&NodalDisplacements>,
        stresses: Option<&NodalStresses>,
    ) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "# vtk DataFile Version 3.0");
        let _ = writeln!(out, "OxiPhysics FE result");
        let _ = writeln!(out, "ASCII");
        let _ = writeln!(out, "DATASET UNSTRUCTURED_GRID");
        let _ = writeln!(out, "POINTS {} double", mesh.nodes.len());
        for n in &mesh.nodes {
            let _ = writeln!(
                out,
                "{:.10e} {:.10e} {:.10e}",
                n.coords[0], n.coords[1], n.coords[2]
            );
        }
        let node_map: HashMap<usize, usize> = mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(i, n)| (n.id, i))
            .collect();
        let total_conn: usize = mesh.elements.iter().map(|e| e.connectivity.len() + 1).sum();
        let _ = writeln!(out, "CELLS {} {}", mesh.elements.len(), total_conn);
        for e in &mesh.elements {
            let ids: Vec<String> = e
                .connectivity
                .iter()
                .map(|nid| {
                    node_map
                        .get(nid)
                        .map(|i| i.to_string())
                        .unwrap_or_else(|| "0".to_string())
                })
                .collect();
            let _ = writeln!(out, "{} {}", e.connectivity.len(), ids.join(" "));
        }
        let _ = writeln!(out, "CELL_TYPES {}", mesh.elements.len());
        for e in &mesh.elements {
            let code = Self::fe_to_vtk_code(&e.element_type) as u32;
            let _ = writeln!(out, "{}", code);
        }
        let n_pts = mesh.nodes.len();
        if displacements.is_some() || stresses.is_some() {
            let _ = writeln!(out, "POINT_DATA {}", n_pts);
        }
        if let Some(disp) = displacements {
            let _ = writeln!(out, "VECTORS displacement double");
            let disp_map: HashMap<usize, &[f64; 3]> = disp
                .node_ids
                .iter()
                .zip(disp.displacements.iter())
                .map(|(id, d)| (*id, d))
                .collect();
            for n in &mesh.nodes {
                let d = disp_map.get(&n.id).copied().unwrap_or(&[0.0; 3]);
                let _ = writeln!(out, "{:.10e} {:.10e} {:.10e}", d[0], d[1], d[2]);
            }
        }
        if let Some(stress) = stresses {
            let _ = writeln!(out, "SCALARS von_mises double 1");
            let _ = writeln!(out, "LOOKUP_TABLE default");
            let stress_map: HashMap<usize, usize> = stress
                .node_ids
                .iter()
                .enumerate()
                .map(|(i, id)| (*id, i))
                .collect();
            for n in &mesh.nodes {
                let vm = stress_map
                    .get(&n.id)
                    .map(|&i| stress.von_mises(i))
                    .unwrap_or(0.0);
                let _ = writeln!(out, "{:.10e}", vm);
            }
        }
        out
    }
    fn fe_to_vtk_code(t: &FeElementType) -> VtkCellCode {
        match t {
            FeElementType::Tet4 | FeElementType::Tet10 => VtkCellCode::Tetra,
            FeElementType::Hex8 | FeElementType::Hex20 => VtkCellCode::Hexahedron,
            FeElementType::Tri3 | FeElementType::Tri6 => VtkCellCode::Triangle,
            FeElementType::Quad4 | FeElementType::Quad8 => VtkCellCode::Quad,
            FeElementType::Line2 => VtkCellCode::Line,
            FeElementType::Unknown(_) => VtkCellCode::Tetra,
        }
    }
}
/// VTK cell type codes (legacy VTK format).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VtkCellCode {
    /// VTK_TETRA (4-node tet).
    Tetra = 10,
    /// VTK_HEXAHEDRON (8-node hex).
    Hexahedron = 12,
    /// VTK_TRIANGLE (3-node triangle).
    Triangle = 5,
    /// VTK_QUAD (4-node quad).
    Quad = 9,
    /// VTK_LINE (2-node line).
    Line = 3,
}
/// Simple round-robin mesh partitioner.
#[derive(Debug, Default)]
pub struct RoundRobinPartitioner;
impl RoundRobinPartitioner {
    /// Partition `mesh` into `num_parts` and return one descriptor per part.
    pub fn partition(mesh: &FeMesh, num_parts: usize) -> Vec<MeshPartition> {
        let np = num_parts.max(1);
        let mut parts: Vec<MeshPartition> = (0..np).map(|r| MeshPartition::new(r, np)).collect();
        for (i, n) in mesh.nodes.iter().enumerate() {
            parts[i % np].owned_nodes.push(n.id);
        }
        for (i, e) in mesh.elements.iter().enumerate() {
            parts[i % np].owned_elements.push(e.id);
        }
        parts
    }
}
