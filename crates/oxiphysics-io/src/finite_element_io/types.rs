//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::if_same_then_else, clippy::manual_strip)]
use std::collections::HashMap;
#[allow(unused_imports)]
use std::fmt::Write as _;

use super::functions::DofIndex;
#[allow(unused_imports)]
use super::functions::*;

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
/// ANSYS MAPDL command deck (simplified subset).
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
    /// Parse a simplified MAPDL text input deck.
    ///
    /// Handles: `NBLOCK`, `EBLOCK`, `MP,EX`, `MP,NUXY`, `MP,DENS`, `D`, `F`.
    pub fn parse(source: &str) -> Self {
        let mut deck = Self::new();
        let mut in_nblock = false;
        let mut in_eblock = false;
        let mut _nblock_count = 0usize;
        let mut _eblock_count = 0usize;
        for line in source.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('!') {
                continue;
            }
            let upper = trimmed.to_uppercase();
            if upper.starts_with("NBLOCK") {
                in_nblock = true;
                in_eblock = false;
                continue;
            }
            if upper.starts_with("EBLOCK") {
                in_eblock = true;
                in_nblock = false;
                continue;
            }
            if upper.starts_with("-1") {
                in_nblock = false;
                in_eblock = false;
                continue;
            }
            if in_nblock {
                if let Some(n) = Self::parse_nblock_line(trimmed) {
                    _nblock_count += 1;
                    deck.nodes.push(n);
                }
                continue;
            }
            if in_eblock {
                if let Some(e) = Self::parse_eblock_line(trimmed) {
                    _eblock_count += 1;
                    deck.elements.push(e);
                }
                continue;
            }
            if upper.starts_with("MP,") {
                Self::apply_mp_command(trimmed, &mut deck);
                continue;
            }
            if upper.starts_with("D,") {
                Self::apply_d_command(trimmed, &mut deck);
                continue;
            }
            if upper.starts_with("F,") {
                Self::apply_f_command(trimmed, &mut deck);
                continue;
            }
        }
        deck
    }
    fn parse_nblock_line(line: &str) -> Option<FeNode> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return None;
        }
        let id: usize = parts[0].parse().ok()?;
        let x: f64 = parts[1].parse().ok()?;
        let y: f64 = parts[2].parse().ok()?;
        let z: f64 = parts[3].parse().ok()?;
        Some(FeNode::new(id, [x, y, z]))
    }
    fn parse_eblock_line(line: &str) -> Option<FeElement> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 {
            return None;
        }
        let id: usize = parts[4].parse().ok()?;
        let connectivity: Vec<usize> = parts[5..]
            .iter()
            .filter_map(|s| s.parse::<usize>().ok())
            .collect();
        Some(FeElement::new(id, FeElementType::Hex8, connectivity))
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
    /// Parse a Gmsh `.msh` ASCII file (v2 or v4).
    pub fn parse(source: &str) -> Self {
        let mut mesh = Self::new();
        let lines: Vec<&str> = source.lines().collect();
        let mut i = 0;
        while i < lines.len() {
            let line = lines[i].trim();
            match line {
                "$MeshFormat" => {
                    i += 1;
                    if i < lines.len() {
                        let ver_line = lines[i].trim();
                        let parts: Vec<&str> = ver_line.split_whitespace().collect();
                        if let Some(v) = parts.first() {
                            if v.starts_with("4") {
                                mesh.version = Some(GmshVersion::V4);
                            } else {
                                mesh.version = Some(GmshVersion::V2);
                            }
                        }
                    }
                    while i < lines.len() && lines[i].trim() != "$EndMeshFormat" {
                        i += 1;
                    }
                }
                "$PhysicalNames" => {
                    i += 1;
                    if i < lines.len() {
                        let _count: usize = lines[i].trim().parse().unwrap_or(0);
                        i += 1;
                    }
                    while i < lines.len() && lines[i].trim() != "$EndPhysicalNames" {
                        let parts: Vec<&str> = lines[i].trim().splitn(3, ' ').collect();
                        if parts.len() == 3 {
                            let dim: i32 = parts[0].parse().unwrap_or(0);
                            let tag: i32 = parts[1].parse().unwrap_or(0);
                            let name = parts[2].trim_matches('"').to_string();
                            mesh.physical_names.insert((dim, tag), name);
                        }
                        i += 1;
                    }
                }
                "$Nodes" => {
                    let version = mesh.version.unwrap_or(GmshVersion::V2);
                    i += 1;
                    match version {
                        GmshVersion::V2 => {
                            while i < lines.len() && lines[i].trim() != "$EndNodes" {
                                let l = lines[i].trim();
                                let parts: Vec<&str> = l.split_whitespace().collect();
                                if parts.len() >= 4
                                    && let (Ok(id), Ok(x), Ok(y), Ok(z)) = (
                                        parts[0].parse::<usize>(),
                                        parts[1].parse::<f64>(),
                                        parts[2].parse::<f64>(),
                                        parts[3].parse::<f64>(),
                                    )
                                {
                                    mesh.nodes.push(FeNode::new(id, [x, y, z]));
                                }
                                i += 1;
                            }
                        }
                        GmshVersion::V4 => {
                            i += 1;
                            while i < lines.len() && lines[i].trim() != "$EndNodes" {
                                let l = lines[i].trim();
                                let parts: Vec<&str> = l.split_whitespace().collect();
                                if parts.len() >= 4
                                    && let (Ok(id), Ok(x), Ok(y), Ok(z)) = (
                                        parts[0].parse::<usize>(),
                                        parts[1].parse::<f64>(),
                                        parts[2].parse::<f64>(),
                                        parts[3].parse::<f64>(),
                                    )
                                {
                                    mesh.nodes.push(FeNode::new(id, [x, y, z]));
                                }
                                i += 1;
                            }
                        }
                    }
                }
                "$Elements" => {
                    i += 1;
                    while i < lines.len() && lines[i].trim() != "$EndElements" {
                        let l = lines[i].trim();
                        let parts: Vec<&str> = l.split_whitespace().collect();
                        if parts.len() >= 3
                            && let (Ok(id), Ok(type_tag)) =
                                (parts[0].parse::<usize>(), parts[1].parse::<usize>())
                        {
                            let ntags: usize = parts[2].parse().unwrap_or(0);
                            let data_start = 3 + ntags;
                            let elem_type = Self::gmsh_type(type_tag);
                            let connectivity: Vec<usize> = parts[data_start..]
                                .iter()
                                .filter_map(|s| s.parse::<usize>().ok())
                                .collect();
                            if !connectivity.is_empty() {
                                mesh.elements
                                    .push(FeElement::new(id, elem_type, connectivity));
                            }
                        }
                        i += 1;
                    }
                }
                _ => {}
            }
            i += 1;
        }
        mesh
    }
    fn gmsh_type(tag: usize) -> FeElementType {
        match tag {
            1 => FeElementType::Line2,
            2 => FeElementType::Tri3,
            3 => FeElementType::Quad4,
            4 => FeElementType::Tet4,
            5 => FeElementType::Hex8,
            9 => FeElementType::Tri6,
            10 => FeElementType::Quad8,
            11 => FeElementType::Tet10,
            _ => FeElementType::Unknown(format!("gmsh:{}", tag)),
        }
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
        match t {
            FeElementType::Line2 => 1,
            FeElementType::Tri3 => 2,
            FeElementType::Quad4 => 3,
            FeElementType::Tet4 => 4,
            FeElementType::Hex8 => 5,
            FeElementType::Tri6 => 9,
            FeElementType::Quad8 => 10,
            FeElementType::Tet10 => 11,
            _ => 0,
        }
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
    /// Parse NASTRAN free-field or fixed-field (8-char) Bulk Data from a string.
    pub fn parse(source: &str) -> Self {
        let mut data = Self::new();
        let mut in_bulk = false;
        for line in source.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('$') {
                continue;
            }
            let upper = trimmed.to_uppercase();
            if upper.starts_with("BEGIN BULK") {
                in_bulk = true;
                continue;
            }
            if upper.starts_with("ENDDATA") {
                break;
            }
            if !in_bulk
                && !upper.starts_with("NASTRAN")
                && !upper.starts_with("SOL")
                && !upper.starts_with("CEND")
                && !upper.starts_with("TITLE")
                && !upper.starts_with("SUBTITLE")
                && !upper.starts_with("LOAD")
                && !upper.starts_with("SPC")
                && !upper.starts_with("SUBCASE")
            {
                continue;
            }
            if upper.starts_with("GRID") {
                if let Some(g) = Self::parse_grid(trimmed) {
                    data.grids.push(g);
                }
            } else if upper.starts_with("CTETRA") {
                if let Some(e) = Self::parse_ctetra(trimmed) {
                    data.elements.push(e);
                }
            } else if upper.starts_with("CHEXA") {
                if let Some(e) = Self::parse_chexa(trimmed) {
                    data.elements.push(e);
                }
            } else if upper.starts_with("MAT1") {
                if let Some(m) = Self::parse_mat1(trimmed) {
                    data.materials.push(m);
                }
            } else if upper.starts_with("SPC ") || upper.starts_with("SPC,") {
                if let Some(s) = Self::parse_spc(trimmed) {
                    data.spcs.push(s);
                }
            } else if upper.starts_with("FORCE")
                && let Some(f) = Self::parse_force(trimmed)
            {
                data.forces.push(f);
            }
        }
        data
    }
    fn split_bulk(line: &str) -> Vec<String> {
        if line.contains(',') {
            return line.split(',').map(|s| s.trim().to_string()).collect();
        }
        let mut fields = Vec::new();
        let chars: Vec<char> = line.chars().collect();
        let mut col = 0;
        while col < chars.len() {
            let end = (col + 8).min(chars.len());
            let field: String = chars[col..end].iter().collect();
            fields.push(field.trim().to_string());
            col += 8;
        }
        fields
    }
    fn parse_grid(line: &str) -> Option<NastranGrid> {
        let f = Self::split_bulk(line);
        if f.len() < 6 {
            return None;
        }
        let id: usize = f[1].parse().ok()?;
        let cp: usize = f[2].parse().unwrap_or(0);
        let x: f64 = f[3].parse().ok()?;
        let y: f64 = f[4].parse().ok()?;
        let z: f64 = f[5].parse().ok()?;
        let cd: usize = f.get(6).and_then(|s| s.parse().ok()).unwrap_or(0);
        Some(NastranGrid {
            id,
            cp,
            coords: [x, y, z],
            cd,
        })
    }
    fn parse_ctetra(line: &str) -> Option<NastranElement> {
        let f = Self::split_bulk(line);
        if f.len() < 7 {
            return None;
        }
        let eid: usize = f[1].parse().ok()?;
        let pid: usize = f[2].parse().ok()?;
        let nodes: Vec<usize> = f[3..7]
            .iter()
            .filter_map(|s| s.parse::<usize>().ok())
            .collect();
        if nodes.len() < 4 {
            return None;
        }
        Some(NastranElement {
            eid,
            pid,
            element_type: FeElementType::Tet4,
            nodes,
        })
    }
    fn parse_chexa(line: &str) -> Option<NastranElement> {
        let f = Self::split_bulk(line);
        if f.len() < 10 {
            return None;
        }
        let eid: usize = f[1].parse().ok()?;
        let pid: usize = f[2].parse().ok()?;
        let nodes: Vec<usize> = f[3..11]
            .iter()
            .filter_map(|s| s.parse::<usize>().ok())
            .collect();
        Some(NastranElement {
            eid,
            pid,
            element_type: FeElementType::Hex8,
            nodes,
        })
    }
    fn parse_mat1(line: &str) -> Option<NastranMat1> {
        let f = Self::split_bulk(line);
        if f.len() < 4 {
            return None;
        }
        let mid: usize = f[1].parse().ok()?;
        let e: f64 = f[2].parse().unwrap_or(0.0);
        let g: f64 = f[3].parse().unwrap_or(0.0);
        let nu: f64 = f.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.3);
        let rho: f64 = f.get(5).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        Some(NastranMat1 { mid, e, g, nu, rho })
    }
    fn parse_spc(line: &str) -> Option<NastranSpc> {
        let f = Self::split_bulk(line);
        if f.len() < 4 {
            return None;
        }
        let sid: usize = f[1].parse().ok()?;
        let gid: usize = f[2].parse().ok()?;
        let components = f[3].clone();
        let value: f64 = f.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        Some(NastranSpc {
            sid,
            gid,
            components,
            value,
        })
    }
    fn parse_force(line: &str) -> Option<NastranForce> {
        let f = Self::split_bulk(line);
        if f.len() < 8 {
            return None;
        }
        let sid: usize = f[1].parse().ok()?;
        let gid: usize = f[2].parse().ok()?;
        let _cid: usize = f[3].parse().unwrap_or(0);
        let scale: f64 = f[4].parse().ok()?;
        let nx: f64 = f[5].parse().ok()?;
        let ny: f64 = f[6].parse().ok()?;
        let nz: f64 = f[7].parse().ok()?;
        Some(NastranForce {
            sid,
            gid,
            scale,
            direction: [nx, ny, nz],
        })
    }
    /// Write bulk data in free-field comma-separated format.
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
    /// Supported keywords: `*NODE`, `*ELEMENT`, `*NSET`, `*ELSET`,
    /// `*MATERIAL`, `*ELASTIC`, `*BOUNDARY`, `*CLOAD`, `*STEP`, `*END STEP`.
    pub fn parse(source: &str) -> FeMesh {
        let mut mesh = FeMesh::new();
        let mut state = AbaqusParseState::default();
        for raw_line in source.lines() {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with("**") {
                continue;
            }
            if line.starts_with('*') {
                let kw_upper = line.to_uppercase();
                state.node_section = false;
                state.element_section = false;
                state.node_set_section = None;
                state.elem_set_section = None;
                state.in_elastic = false;
                if kw_upper.starts_with("*NODE") && !kw_upper.starts_with("*NSET") {
                    state.node_section = true;
                } else if kw_upper.starts_with("*ELEMENT") {
                    state.element_section = true;
                    state.current_elem_type = Self::parse_element_type(line);
                } else if kw_upper.starts_with("*NSET") {
                    let name = Self::parse_keyword_param(line, "NSET");
                    state.node_set_section = Some(name);
                } else if kw_upper.starts_with("*ELSET") {
                    let name = Self::parse_keyword_param(line, "ELSET");
                    state.elem_set_section = Some(name);
                } else if kw_upper.starts_with("*MATERIAL") {
                    let name = Self::parse_keyword_param(line, "NAME");
                    state.material_name = Some(name);
                } else if kw_upper.starts_with("*ELASTIC") {
                    state.in_elastic = true;
                } else if kw_upper.starts_with("*STEP") {
                    let step_name = Self::parse_keyword_param(line, "NAME");
                    mesh.steps.push(AnalysisStep::new(step_name));
                } else if kw_upper.starts_with("*BOUNDARY") {
                } else if kw_upper.starts_with("*CLOAD") {
                }
                continue;
            }
            if state.node_section {
                if let Some(node) = Self::parse_node_line(line) {
                    mesh.nodes.push(node);
                }
            } else if state.element_section {
                if let Some(elem) = Self::parse_element_line(line, state.current_elem_type.clone())
                {
                    mesh.elements.push(elem);
                }
            } else if let Some(ref set_name) = state.node_set_section.clone() {
                let ids: Vec<usize> = line
                    .split(',')
                    .filter_map(|s| s.trim().parse::<usize>().ok())
                    .collect();
                mesh.node_sets
                    .entry(set_name.clone())
                    .or_default()
                    .extend(ids);
            } else if let Some(ref set_name) = state.elem_set_section.clone() {
                let ids: Vec<usize> = line
                    .split(',')
                    .filter_map(|s| s.trim().parse::<usize>().ok())
                    .collect();
                mesh.element_sets
                    .entry(set_name.clone())
                    .or_default()
                    .extend(ids);
            } else if state.in_elastic
                && let Some(mat_name) = &state.material_name.clone()
            {
                if let Some(mat) = Self::parse_elastic_line(line, mat_name) {
                    mesh.materials.push(mat);
                }
                state.in_elastic = false;
            }
        }
        mesh
    }
    fn parse_node_line(line: &str) -> Option<FeNode> {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            return None;
        }
        let id: usize = parts[0].trim().parse().ok()?;
        let x: f64 = parts[1].trim().parse().ok()?;
        let y: f64 = parts[2].trim().parse().ok()?;
        let z: f64 = parts[3].trim().parse().ok()?;
        Some(FeNode::new(id, [x, y, z]))
    }
    fn parse_element_line(line: &str, elem_type: FeElementType) -> Option<FeElement> {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 2 {
            return None;
        }
        let id: usize = parts[0].trim().parse().ok()?;
        let connectivity: Vec<usize> = parts[1..]
            .iter()
            .filter_map(|s| s.trim().parse::<usize>().ok())
            .collect();
        if connectivity.is_empty() {
            return None;
        }
        Some(FeElement::new(id, elem_type, connectivity))
    }
    fn parse_elastic_line(line: &str, name: &str) -> Option<LinearElasticMaterial> {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 2 {
            return None;
        }
        let e: f64 = parts[0].trim().parse().ok()?;
        let nu: f64 = parts[1].trim().parse().ok()?;
        Some(LinearElasticMaterial::new(0, name, e, nu, 0.0))
    }
    fn parse_element_type(line: &str) -> FeElementType {
        let up = line.to_uppercase();
        if up.contains("C3D4") || up.contains("TET4") {
            FeElementType::Tet4
        } else if up.contains("C3D10") || up.contains("TET10") {
            FeElementType::Tet10
        } else if up.contains("C3D8") || up.contains("HEX8") {
            FeElementType::Hex8
        } else if up.contains("C3D20") || up.contains("HEX20") {
            FeElementType::Hex20
        } else if up.contains("S3") || up.contains("TRI3") {
            FeElementType::Tri3
        } else if up.contains("S4") || up.contains("QUAD4") {
            FeElementType::Quad4
        } else if up.contains("T3D2") {
            FeElementType::Line2
        } else {
            FeElementType::Unknown(line.to_string())
        }
    }
    fn parse_keyword_param(line: &str, param: &str) -> String {
        let up = line.to_uppercase();
        let search = format!("{}=", param.to_uppercase());
        if let Some(pos) = up.find(&search) {
            let rest = &line[pos + search.len()..];
            rest.split(',').next().unwrap_or("").trim().to_string()
        } else {
            String::new()
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
            if t.starts_with("TITLE:") {
                title = t["TITLE:".len()..].trim().to_string();
            } else if t.starts_with("STEP:") {
                step = t["STEP:".len()..].trim().parse().unwrap_or(0);
            } else if t.starts_with("TIME:") {
                time = t["TIME:".len()..].trim().parse().unwrap_or(0.0);
            } else if t.starts_with("RESIDUAL_NORM:") {
                residual_norm = t["RESIDUAL_NORM:".len()..].trim().parse().unwrap_or(0.0);
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
