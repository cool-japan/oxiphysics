//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::collections::HashMap;
use std::fmt::Write;

/// VMD molecular representation style.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VmdRepresentation {
    /// Licorice (stick) representation.
    Licorice,
    /// CPK (ball and stick) representation.
    Cpk,
    /// Van der Waals (sphere) representation.
    Vdw,
    /// New ribbon for proteins.
    NewRibbons,
    /// Surf representation (MSMS surface).
    Surf,
    /// Points representation.
    Points,
    /// Lines representation.
    Lines,
    /// Dynamic bonds.
    DynamicBonds,
}
impl VmdRepresentation {
    fn as_str(self) -> &'static str {
        match self {
            VmdRepresentation::Licorice => "Licorice 0.1 12 12",
            VmdRepresentation::Cpk => "CPK 0.7 0.3 12 12",
            VmdRepresentation::Vdw => "VDW 1.0 12",
            VmdRepresentation::NewRibbons => "NewRibbons 0.3 10 1.5 0",
            VmdRepresentation::Surf => "Surf 1.4",
            VmdRepresentation::Points => "Points 1.0",
            VmdRepresentation::Lines => "Lines 1.0",
            VmdRepresentation::DynamicBonds => "DynamicBonds 1.6 0.1 12",
        }
    }
}
/// A named group of atom indices (e.g., a residue, chain, or selection).
#[derive(Debug, Clone)]
pub struct AtomGroup {
    /// Name of this group (e.g., "backbone", "hydrophobic").
    pub name: String,
    /// Atom indices belonging to this group.
    pub indices: Vec<usize>,
    /// Optional display color for this group.
    pub color: Option<[f32; 3]>,
}
impl AtomGroup {
    /// Create an empty group with the given name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            indices: Vec::new(),
            color: None,
        }
    }
    /// Create a group from a name and list of atom indices.
    pub fn from_indices(name: &str, indices: Vec<usize>) -> Self {
        Self {
            name: name.to_string(),
            indices,
            color: None,
        }
    }
    /// Add an atom index to this group.
    pub fn add(&mut self, index: usize) {
        self.indices.push(index);
    }
    /// Return the number of atoms in this group.
    pub fn len(&self) -> usize {
        self.indices.len()
    }
    /// Return `true` if the group has no atoms.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}
/// A complete molecular scene with atoms, bonds, and named groups.
///
/// Acts as the central data model consumed by all writers in this module.
#[derive(Debug, Clone)]
pub struct MolecularScene {
    /// Scene / molecule name.
    pub name: String,
    /// All atoms in the scene.
    pub atoms: Vec<MolAtom>,
    /// All bonds in the scene.
    pub bonds: Vec<MolBond>,
    /// Named atom groups (e.g., chains, residues, selections).
    pub groups: HashMap<String, AtomGroup>,
    /// Scene-level scalar properties (e.g., total energy, temperature).
    pub properties: HashMap<String, f64>,
    /// Scene-level string metadata.
    pub metadata: HashMap<String, String>,
    /// Unit cell parameters \[a, b, c, alpha, beta, gamma\] if periodic.
    pub unit_cell: Option<[f64; 6]>,
}
impl MolecularScene {
    /// Create an empty scene with the given name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            groups: HashMap::new(),
            properties: HashMap::new(),
            metadata: HashMap::new(),
            unit_cell: None,
        }
    }
    /// Add an atom and return its index.
    pub fn add_atom(&mut self, atom: MolAtom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        idx
    }
    /// Add a bond.
    pub fn add_bond(&mut self, bond: MolBond) {
        self.bonds.push(bond);
    }
    /// Add a named atom group.
    pub fn add_group(&mut self, group: AtomGroup) {
        self.groups.insert(group.name.clone(), group);
    }
    /// Set a scalar property.
    pub fn set_property(&mut self, key: &str, value: f64) {
        self.properties.insert(key.to_string(), value);
    }
    /// Set a metadata string.
    pub fn set_metadata(&mut self, key: &str, value: &str) {
        self.metadata.insert(key.to_string(), value.to_string());
    }
    /// Set unit cell parameters \[a, b, c, α, β, γ\] (Å and degrees).
    pub fn set_unit_cell(&mut self, params: [f64; 6]) {
        self.unit_cell = Some(params);
    }
    /// Return the atom at the given index, if valid.
    pub fn atom(&self, idx: usize) -> Option<&MolAtom> {
        self.atoms.get(idx)
    }
    /// Return the number of atoms.
    pub fn natoms(&self) -> usize {
        self.atoms.len()
    }
    /// Return the number of bonds.
    pub fn nbonds(&self) -> usize {
        self.bonds.len()
    }
    /// Compute the centroid of all atoms.
    pub fn centroid(&self) -> [f64; 3] {
        if self.atoms.is_empty() {
            return [0.0; 3];
        }
        let n = self.atoms.len() as f64;
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        let mut cz = 0.0_f64;
        for a in &self.atoms {
            cx += a.position[0];
            cy += a.position[1];
            cz += a.position[2];
        }
        [cx / n, cy / n, cz / n]
    }
    /// Compute the axis-aligned bounding box `(min, max)`.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for a in &self.atoms {
            for k in 0..3 {
                if a.position[k] < mn[k] {
                    mn[k] = a.position[k];
                }
                if a.position[k] > mx[k] {
                    mx[k] = a.position[k];
                }
            }
        }
        (mn, mx)
    }
    /// Determine bonds automatically from covalent radii (tolerance 0.4 Å).
    pub fn auto_bonds(&mut self) {
        self.bonds.clear();
        let n = self.atoms.len();
        let tolerance = 0.4_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let r_sum = self.atoms[i].covalent_radius() + self.atoms[j].covalent_radius();
                let dx = self.atoms[i].position[0] - self.atoms[j].position[0];
                let dy = self.atoms[i].position[1] - self.atoms[j].position[1];
                let dz = self.atoms[i].position[2] - self.atoms[j].position[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < r_sum + tolerance {
                    self.bonds.push(MolBond::single(i, j));
                }
            }
        }
    }
    /// Group atoms by chain ID.
    pub fn group_by_chain(&mut self) {
        let mut chains: HashMap<char, Vec<usize>> = HashMap::new();
        for (i, atom) in self.atoms.iter().enumerate() {
            chains.entry(atom.chain_id).or_default().push(i);
        }
        for (ch, indices) in chains {
            let name = format!("chain_{}", ch);
            self.groups
                .insert(name.clone(), AtomGroup::from_indices(&name, indices));
        }
    }
    /// Group atoms by residue sequence number.
    pub fn group_by_residue(&mut self) {
        let mut residues: HashMap<(char, i32, String), Vec<usize>> = HashMap::new();
        for (i, atom) in self.atoms.iter().enumerate() {
            residues
                .entry((atom.chain_id, atom.residue_seq, atom.residue_name.clone()))
                .or_default()
                .push(i);
        }
        for ((ch, seq, resname), indices) in residues {
            let name = format!("{}_{}_{}", ch, seq, resname);
            self.groups
                .insert(name.clone(), AtomGroup::from_indices(&name, indices));
        }
    }
}
/// Generates VMD Tcl scripts for molecular visualization.
///
/// Writes `mol load`, `mol representation`, `mol color`, and `mol material`
/// commands that can be executed in VMD's Tcl interpreter.
#[derive(Debug, Clone)]
pub struct VmdScriptWriter {
    /// Color scheme for the primary representation.
    pub color_scheme: ColorScheme,
    /// Primary representation style (e.g., "Licorice", "CPK", "VDW").
    pub representation: VmdRepresentation,
    /// Whether to add a second representation for the backbone.
    pub show_backbone: bool,
    /// Material name (e.g., "Opaque", "Transparent", "Glossy").
    pub material: String,
    /// Whether to write a `render` command at the end.
    pub add_render: bool,
    /// Output render filename.
    pub render_file: String,
}
impl VmdScriptWriter {
    /// Create a writer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Generate a complete VMD Tcl script for the given scene.
    pub fn write_scene(&self, scene: &MolecularScene) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "# VMD Tcl script generated by OxiPhysics");
        let _ = writeln!(out, "# Scene: {}", scene.name);
        writeln!(
            out,
            "# Atoms: {}  Bonds: {}",
            scene.natoms(),
            scene.nbonds()
        )
        .expect("operation should succeed");
        let _ = writeln!(out);
        let _ = writeln!(out, "# Write inline XYZ and load");
        writeln!(
            out,
            "set tmpfile \"/tmp/oxiphysics_{}.xyz\"",
            sanitize_name(&scene.name)
        )
        .expect("operation should succeed");
        let _ = writeln!(out, "set fh [open $tmpfile w]");
        let _ = writeln!(out, "puts $fh \"{}\"", scene.natoms());
        let _ = writeln!(out, "puts $fh \"{}\"", scene.name);
        for atom in &scene.atoms {
            writeln!(
                out,
                "puts $fh \"{} {:.6} {:.6} {:.6}\"",
                atom.element, atom.position[0], atom.position[1], atom.position[2],
            )
            .expect("operation should succeed");
        }
        let _ = writeln!(out, "close $fh");
        let _ = writeln!(out, "mol load xyz $tmpfile");
        let _ = writeln!(out, "set molid [molinfo top get id]");
        let _ = writeln!(out);
        let _ = writeln!(out, "mol delrep 0 $molid");
        let _ = writeln!(out, "# Primary representation");
        let _ = writeln!(out, "mol representation {}", self.representation.as_str());
        let _ = writeln!(out, "mol color {}", self.vmd_color_method());
        let _ = writeln!(out, "mol material {}", self.material);
        let _ = writeln!(out, "mol selection all");
        let _ = writeln!(out, "mol addrep $molid");
        if self.show_backbone {
            let _ = writeln!(out);
            let _ = writeln!(out, "# Backbone ribbon representation");
            let _ = writeln!(out, "mol representation NewRibbons 0.3 10 1.5 0");
            let _ = writeln!(out, "mol color Chain");
            let _ = writeln!(out, "mol material Opaque");
            let _ = writeln!(out, "mol selection backbone");
            let _ = writeln!(out, "mol addrep $molid");
        }
        for (grp_name, grp) in &scene.groups {
            if !grp.indices.is_empty() {
                let idx_list: Vec<String> = grp.indices.iter().map(|i| i.to_string()).collect();
                let _ = writeln!(out, "# Group: {grp_name}");
                writeln!(
                    out,
                    "set sel_{} [atomselect $molid \"index {}\"]",
                    sanitize_name(grp_name),
                    idx_list.join(" ")
                )
                .expect("operation should succeed");
            }
        }
        let _ = writeln!(out);
        let _ = writeln!(out, "display resetview");
        let _ = writeln!(out, "rotate x by -90");
        if self.add_render {
            let _ = writeln!(out);
            let _ = writeln!(out, "render TachyonInternal {}", self.render_file);
        }
        out
    }
    /// Generate a Tcl script to animate multiple XYZ frames.
    pub fn write_trajectory_script(&self, frames: &[Vec<MolAtom>], scene_name: &str) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "# VMD trajectory script");
        let _ = writeln!(out, "# Frames: {}", frames.len());
        let _ = writeln!(out, "set tmpbase \"/tmp/oxiphysics_traj\"");
        for (fi, frame) in frames.iter().enumerate() {
            let _ = writeln!(out, "# Frame {fi}");
            let _ = writeln!(out, "set ff [open ${{tmpbase}}_{fi}.xyz w]");
            let _ = writeln!(out, "puts $ff \"{}\"", frame.len());
            let _ = writeln!(out, "puts $ff \"frame {fi} {scene_name}\"");
            for atom in frame {
                writeln!(
                    out,
                    "puts $ff \"{} {:.6} {:.6} {:.6}\"",
                    atom.element, atom.position[0], atom.position[1], atom.position[2]
                )
                .expect("operation should succeed");
            }
            let _ = writeln!(out, "close $ff");
        }
        let _ = writeln!(out, "mol load xyz ${{tmpbase}}_0.xyz");
        for fi in 1..frames.len() {
            let _ = writeln!(out, "mol addfile ${{tmpbase}}_{fi}.xyz");
        }
        out
    }
    fn vmd_color_method(&self) -> &'static str {
        match self.color_scheme {
            ColorScheme::ByElement => "Element",
            ColorScheme::ByChain => "Chain",
            ColorScheme::ByResidue => "ResName",
            ColorScheme::ByScalar => "Beta",
            ColorScheme::Uniform => "ColorID 8",
        }
    }
}
/// A chemical bond between two atoms.
#[derive(Debug, Clone)]
pub struct MolBond {
    /// Index of the first atom.
    pub atom1: usize,
    /// Index of the second atom.
    pub atom2: usize,
    /// Bond order.
    pub order: BondOrder,
    /// Whether this bond is in an aromatic ring.
    pub is_aromatic: bool,
    /// Whether this bond is in a ring.
    pub in_ring: bool,
    /// Stereo designation (0 = not stereo, 1 = Up, 6 = Down).
    pub stereo: u8,
}
impl MolBond {
    /// Create a single bond.
    pub fn single(atom1: usize, atom2: usize) -> Self {
        Self {
            atom1,
            atom2,
            order: BondOrder::Single,
            is_aromatic: false,
            in_ring: false,
            stereo: 0,
        }
    }
    /// Create a double bond.
    pub fn double(atom1: usize, atom2: usize) -> Self {
        Self {
            atom1,
            atom2,
            order: BondOrder::Double,
            is_aromatic: false,
            in_ring: false,
            stereo: 0,
        }
    }
    /// Create a triple bond.
    pub fn triple(atom1: usize, atom2: usize) -> Self {
        Self {
            atom1,
            atom2,
            order: BondOrder::Triple,
            is_aromatic: false,
            in_ring: false,
            stereo: 0,
        }
    }
    /// Create an aromatic bond.
    pub fn aromatic(atom1: usize, atom2: usize) -> Self {
        Self {
            atom1,
            atom2,
            order: BondOrder::Aromatic,
            is_aromatic: true,
            in_ring: true,
            stereo: 0,
        }
    }
}
/// A single frame in an XYZ trajectory.
#[derive(Debug, Clone)]
pub struct XyzFrame {
    /// Atom positions \[x, y, z\] in ångströms.
    pub positions: Vec<[f64; 3]>,
    /// Atom velocities \[vx, vy, vz\] in Å/ps (optional).
    pub velocities: Option<Vec<[f64; 3]>>,
    /// Per-atom forces \[fx, fy, fz\] in kcal/mol/Å (optional).
    pub forces: Option<Vec<[f64; 3]>>,
    /// Element symbols for each atom.
    pub elements: Vec<String>,
    /// Simulation time in picoseconds.
    pub time_ps: f64,
    /// Step number.
    pub step: usize,
    /// Total energy (optional).
    pub total_energy: Option<f64>,
    /// Temperature in K (optional).
    pub temperature: Option<f64>,
}
impl XyzFrame {
    /// Create a frame from positions and elements.
    pub fn new(positions: Vec<[f64; 3]>, elements: Vec<String>) -> Self {
        assert_eq!(
            positions.len(),
            elements.len(),
            "positions and elements must have the same length"
        );
        Self {
            positions,
            velocities: None,
            forces: None,
            elements,
            time_ps: 0.0,
            step: 0,
            total_energy: None,
            temperature: None,
        }
    }
    /// Build an XyzFrame from a MolecularScene.
    pub fn from_scene(scene: &MolecularScene, step: usize, time_ps: f64) -> Self {
        let positions: Vec<[f64; 3]> = scene.atoms.iter().map(|a| a.position).collect();
        let elements: Vec<String> = scene.atoms.iter().map(|a| a.element.clone()).collect();
        Self {
            positions,
            velocities: None,
            forces: None,
            elements,
            time_ps,
            step,
            total_energy: None,
            temperature: None,
        }
    }
    /// Return the number of atoms in this frame.
    pub fn natoms(&self) -> usize {
        self.positions.len()
    }
}
/// Generates PyMOL .pml command scripts for molecular visualization.
///
/// Writes `load`, `show`, `color`, `set` and selection commands that can be
/// executed directly in PyMOL.
#[derive(Debug, Clone)]
pub struct PymolScriptWriter {
    /// Color scheme to use when generating color commands.
    pub color_scheme: ColorScheme,
    /// Whether to write sphere representation in addition to sticks.
    pub show_spheres: bool,
    /// Whether to include a `ray` command at the end for ray-tracing.
    pub add_ray: bool,
    /// Background color name (PyMOL color string, e.g., "white", "black").
    pub background: String,
    /// Stick radius for `show sticks`.
    pub stick_radius: f64,
    /// Sphere scale for `show spheres`.
    pub sphere_scale: f64,
    /// Whether to write `set orthoscopic` for orthographic projection.
    pub orthoscopic: bool,
}
impl PymolScriptWriter {
    /// Create a writer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Generate a complete PyMOL script for the given scene.
    pub fn write_scene(&self, scene: &MolecularScene) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "# PyMOL script generated by OxiPhysics");
        let _ = writeln!(out, "# Scene: {}", scene.name);
        writeln!(
            out,
            "# Atoms: {}  Bonds: {}",
            scene.natoms(),
            scene.nbonds()
        )
        .expect("operation should succeed");
        let _ = writeln!(out);
        let _ = writeln!(out, "bg_color {}", self.background);
        if self.orthoscopic {
            let _ = writeln!(out, "set orthoscopic, 1");
        }
        let obj_name = sanitize_name(&scene.name);
        let _ = writeln!(out, "# Create molecule object");
        let _ = writeln!(out, "create {obj_name}, none");
        for atom in &scene.atoms {
            writeln!(
                out,
                "pseudoatom {obj_name}, name={}, elem={}, resi={}, chain={}, pos=[{:.4},{:.4},{:.4}]",
                atom.atom_name, atom.element, atom.residue_seq, atom.chain_id, atom
                .position[0], atom.position[1], atom.position[2],
            )
                .expect("value should be present");
        }
        let _ = writeln!(out);
        let _ = writeln!(out, "hide everything, {obj_name}");
        let _ = writeln!(out, "show sticks, {obj_name}");
        writeln!(
            out,
            "set stick_radius, {:.4}, {obj_name}",
            self.stick_radius
        )
        .expect("operation should succeed");
        if self.show_spheres {
            let _ = writeln!(out, "show spheres, {obj_name}");
            writeln!(
                out,
                "set sphere_scale, {:.4}, {obj_name}",
                self.sphere_scale
            )
            .expect("operation should succeed");
        }
        self.write_colors(&mut out, scene, &obj_name);
        for (grp_name, grp) in &scene.groups {
            if !grp.indices.is_empty() {
                let sel_name = sanitize_name(grp_name);
                let indices_str: Vec<String> =
                    grp.indices.iter().map(|i| (i + 1).to_string()).collect();
                writeln!(
                    out,
                    "select {sel_name}, {obj_name} and index {}",
                    indices_str.join("+")
                )
                .expect("operation should succeed");
                if let Some(c) = grp.color {
                    writeln!(
                        out,
                        "set_color color_{sel_name}, [{:.2},{:.2},{:.2}]",
                        c[0], c[1], c[2]
                    )
                    .expect("operation should succeed");
                    let _ = writeln!(out, "color color_{sel_name}, {sel_name}");
                }
            }
        }
        let _ = writeln!(out);
        let _ = writeln!(out, "zoom {obj_name}");
        if self.add_ray {
            let _ = writeln!(out, "ray");
        }
        out
    }
    /// Generate a script that colors by B-factor / scalar_property using a spectrum.
    pub fn write_bfactor_spectrum(&self, scene: &MolecularScene) -> String {
        let mut out = String::new();
        let obj_name = sanitize_name(&scene.name);
        let _ = writeln!(out, "# B-factor spectrum script");
        let values: Vec<f64> = scene.atoms.iter().map(|a| a.scalar_property).collect();
        let min_v = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_v = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let _ = writeln!(out, "# scalar property range: {min_v:.4} .. {max_v:.4}");
        writeln!(
            out,
            "spectrum b, rainbow, {obj_name}, minimum={min_v:.4}, maximum={max_v:.4}"
        )
        .expect("operation should succeed");
        out
    }
    /// Write a PyMOL script that generates individual atom color commands.
    pub fn write_per_atom_colors(&self, scene: &MolecularScene) -> String {
        let mut out = String::new();
        let obj_name = sanitize_name(&scene.name);
        for (i, atom) in scene.atoms.iter().enumerate() {
            let c = atom.color;
            writeln!(
                out,
                "set_color atom_color_{i}, [{:.3},{:.3},{:.3}]",
                c[0], c[1], c[2]
            )
            .expect("operation should succeed");
            let _ = writeln!(out, "color atom_color_{i}, {obj_name} and index {}", i + 1);
        }
        out
    }
    fn write_colors(&self, out: &mut String, _scene: &MolecularScene, obj_name: &str) {
        match self.color_scheme {
            ColorScheme::ByElement => {
                let _ = writeln!(out, "util.cbaw {obj_name}");
            }
            ColorScheme::ByChain => {
                let _ = writeln!(out, "util.cbc {obj_name}");
            }
            ColorScheme::ByResidue => {
                let _ = writeln!(out, "util.cbss {obj_name}");
            }
            ColorScheme::ByScalar => {
                let _ = writeln!(out, "spectrum b, rainbow, {obj_name}");
            }
            ColorScheme::Uniform => {
                let _ = writeln!(out, "color gray, {obj_name}");
            }
        }
    }
}
/// Writes molecular data in the Chemical JSON (CJ) format.
///
/// Chemical JSON is an open format used by Avogadro and related tools to store
/// atoms, bonds, properties, and molecular geometry in JSON.
#[derive(Debug, Clone)]
pub struct ChemicalJsonWriter {
    /// Whether to include partial charges in output.
    pub include_charges: bool,
    /// Whether to include the unit cell in output.
    pub include_unit_cell: bool,
    /// Whether to pretty-print the JSON output.
    pub pretty_print: bool,
    /// Indent string used when pretty-printing.
    pub indent: String,
}
impl ChemicalJsonWriter {
    /// Create a writer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Write the scene as a Chemical JSON string.
    pub fn write_scene(&self, scene: &MolecularScene) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "{{");
        let i1 = &self.indent;
        let _ = writeln!(out, "{i1}\"chemicalJson\": 1,");
        let _ = writeln!(out, "{i1}\"name\": \"{}\",", escape_json(&scene.name));
        let _ = writeln!(out, "{i1}\"atoms\": {{");
        let i2 = format!("{i1}{i1}");
        let _ = writeln!(out, "{i2}\"elements\": {{");
        writeln!(
            out,
            "{i2}{i1}\"number\": [{}]",
            self.elements_numbers(scene)
        )
        .expect("operation should succeed");
        let _ = writeln!(out, "{i2}}},");
        let _ = writeln!(out, "{i2}\"coords\": {{");
        let _ = writeln!(out, "{i2}{i1}\"3d\": [{}]", self.coords_3d(scene));
        let _ = writeln!(out, "{i2}}},");
        writeln!(
            out,
            "{i2}\"formalCharges\": [{}],",
            self.formal_charges(scene)
        )
        .expect("operation should succeed");
        let _ = writeln!(out, "{i2}\"labels\": [{}],", self.atom_labels(scene));
        let _ = writeln!(out, "{i2}\"layer\": {{");
        let _ = writeln!(out, "{i2}{i1}\"scalars\": [{}]", self.scalars(scene));
        let _ = writeln!(out, "{i2}}}");
        let _ = writeln!(out, "{i1}}},");
        let _ = writeln!(out, "{i1}\"bonds\": {{");
        let _ = writeln!(out, "{i2}\"connections\": {{");
        let _ = writeln!(out, "{i2}{i1}\"index\": [{}]", self.bond_indices(scene));
        let _ = writeln!(out, "{i2}}},");
        let _ = writeln!(out, "{i2}\"order\": [{}]", self.bond_orders(scene));
        let _ = writeln!(out, "{i1}}},");
        let _ = writeln!(out, "{i1}\"properties\": {{");
        let mut props: Vec<String> = scene
            .properties
            .iter()
            .map(|(k, v)| format!("{i2}\"{}\": {:.6}", escape_json(k), v))
            .collect();
        props.sort();
        let _ = writeln!(out, "{}", props.join(",\n"));
        let _ = writeln!(out, "{i1}}}");
        let _ = write!(out, "}}");
        out
    }
    fn elements_numbers(&self, scene: &MolecularScene) -> String {
        scene
            .atoms
            .iter()
            .map(|a| atomic_number(&a.element).to_string())
            .collect::<Vec<_>>()
            .join(", ")
    }
    fn coords_3d(&self, scene: &MolecularScene) -> String {
        scene
            .atoms
            .iter()
            .flat_map(|a| a.position.iter().cloned())
            .map(|v| format!("{v:.6}"))
            .collect::<Vec<_>>()
            .join(", ")
    }
    fn formal_charges(&self, scene: &MolecularScene) -> String {
        if self.include_charges {
            scene
                .atoms
                .iter()
                .map(|a| a.partial_charge.unwrap_or(0.0))
                .map(|c| format!("{c:.4}"))
                .collect::<Vec<_>>()
                .join(", ")
        } else {
            scene
                .atoms
                .iter()
                .map(|_| "0")
                .collect::<Vec<_>>()
                .join(", ")
        }
    }
    fn atom_labels(&self, scene: &MolecularScene) -> String {
        scene
            .atoms
            .iter()
            .map(|a| format!("\"{}\"", escape_json(&a.atom_name)))
            .collect::<Vec<_>>()
            .join(", ")
    }
    fn scalars(&self, scene: &MolecularScene) -> String {
        scene
            .atoms
            .iter()
            .map(|a| format!("{:.6}", a.scalar_property))
            .collect::<Vec<_>>()
            .join(", ")
    }
    fn bond_indices(&self, scene: &MolecularScene) -> String {
        scene
            .bonds
            .iter()
            .flat_map(|b| [b.atom1, b.atom2])
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    }
    fn bond_orders(&self, scene: &MolecularScene) -> String {
        scene
            .bonds
            .iter()
            .map(|b| b.order.as_int().to_string())
            .collect::<Vec<_>>()
            .join(", ")
    }
}
/// Writes molecules in MDL Molfile format (.mol) using the V2000 connection table.
///
/// The V2000 format is the most widely supported molfile format and is accepted
/// by RDKit, OpenBabel, MarvinSketch, and most cheminformatics toolkits.
#[derive(Debug, Clone)]
pub struct MolfileWriter {
    /// Author / program name written to the header.
    pub author: String,
    /// Whether to write SDF record separator (`$$$$`) after each molecule.
    pub sdf_mode: bool,
    /// Whether to include stereo bond information.
    pub include_stereo: bool,
}
impl MolfileWriter {
    /// Create a writer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Write a single molecule as MDL V2000 Molfile.
    pub fn write_scene(&self, scene: &MolecularScene) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "{}", scene.name);
        let _ = writeln!(out, "  {}  3D", self.author);
        let _ = writeln!(out);
        let natoms = scene.natoms();
        let nbonds = scene.nbonds();
        if natoms > 999 || nbonds > 999 {
            let _ = writeln!(out, "# WARNING: V2000 supports max 999 atoms/bonds");
        }
        let natoms_clamped = natoms.min(999);
        let nbonds_clamped = nbonds.min(999);
        writeln!(
            out,
            "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
            natoms_clamped, nbonds_clamped
        )
        .expect("operation should succeed");
        for atom in &scene.atoms {
            let symbol = format_element_symbol(&atom.element);
            writeln!(
                out,
                "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                atom.position[0], atom.position[1], atom.position[2], symbol,
            )
            .expect("operation should succeed");
        }
        for bond in &scene.bonds {
            let stereo = if self.include_stereo { bond.stereo } else { 0 };
            writeln!(
                out,
                "{:>3}{:>3}{:>3}{:>3}  0  0  0",
                bond.atom1 + 1,
                bond.atom2 + 1,
                bond.order.as_int(),
                stereo,
            )
            .expect("operation should succeed");
        }
        let _ = writeln!(out, "M  END");
        if self.sdf_mode {
            for (key, val) in &scene.properties {
                let _ = writeln!(out, ">  <{key}>");
                let _ = writeln!(out, "{val:.6}");
                let _ = writeln!(out);
            }
            let _ = writeln!(out, "$$$$");
        }
        out
    }
    /// Write multiple scenes as an SDF file.
    pub fn write_sdf(&self, scenes: &[MolecularScene]) -> String {
        let sdf_writer = Self {
            sdf_mode: true,
            ..self.clone()
        };
        scenes
            .iter()
            .map(|s| sdf_writer.write_scene(s))
            .collect::<Vec<_>>()
            .join("")
    }
}
/// A single atom in a molecular scene.
#[derive(Debug, Clone)]
pub struct MolAtom {
    /// Atom index (0-based, unique within the scene).
    pub index: usize,
    /// Element symbol (e.g., "C", "O", "N").
    pub element: String,
    /// 3-D position in ångströms \[x, y, z\].
    pub position: [f64; 3],
    /// Display radius in ångströms (defaults to van-der-Waals radius).
    pub radius: f64,
    /// Display color \[r, g, b\] in \[0, 1\] (defaults to CPK color).
    pub color: [f32; 3],
    /// Atom name as in PDB convention (e.g., "CA", "CB").
    pub atom_name: String,
    /// Residue name (e.g., "ALA", "HOH").
    pub residue_name: String,
    /// Residue sequence number.
    pub residue_seq: i32,
    /// Chain identifier.
    pub chain_id: char,
    /// Partial charge (optional, used in Chemical JSON).
    pub partial_charge: Option<f64>,
    /// Arbitrary scalar property (e.g., B-factor, RMSD).
    pub scalar_property: f64,
}
impl MolAtom {
    /// Create an atom with default display properties derived from its element.
    pub fn new(index: usize, element: &str, position: [f64; 3]) -> Self {
        let (_cov_radius, color) = element_data(element);
        let vdw = vdw_radius(element);
        Self {
            index,
            element: element.to_string(),
            position,
            radius: vdw,
            color,
            atom_name: format!("{}{}", element, index + 1),
            residue_name: "MOL".to_string(),
            residue_seq: 1,
            chain_id: 'A',
            partial_charge: None,
            scalar_property: 0.0,
        }
    }
    /// Create an atom with fully specified display properties.
    #[allow(clippy::too_many_arguments)]
    pub fn with_display(
        index: usize,
        element: &str,
        position: [f64; 3],
        radius: f64,
        color: [f32; 3],
        atom_name: &str,
        residue_name: &str,
        residue_seq: i32,
        chain_id: char,
    ) -> Self {
        Self {
            index,
            element: element.to_string(),
            position,
            radius,
            color,
            atom_name: atom_name.to_string(),
            residue_name: residue_name.to_string(),
            residue_seq,
            chain_id,
            partial_charge: None,
            scalar_property: 0.0,
        }
    }
    /// Return the covalent radius for this atom's element.
    pub fn covalent_radius(&self) -> f64 {
        element_data(&self.element).0
    }
    /// Return the default CPK color for this atom's element.
    pub fn cpk_color(&self) -> [f32; 3] {
        element_data(&self.element).1
    }
}
/// Color scheme used when writing molecular scripts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColorScheme {
    /// Color atoms by element (CPK coloring).
    ByElement,
    /// Color atoms by residue type.
    ByResidue,
    /// Color atoms by chain.
    ByChain,
    /// Color atoms by B-factor / scalar property.
    ByScalar,
    /// Use a single uniform color.
    Uniform,
}
/// Writes multi-frame XYZ trajectories including optional velocities and forces.
///
/// The extended XYZ comment line includes time, step, energy, and temperature
/// as key=value pairs following the OVITO / ASE convention.
#[derive(Debug, Clone, Default)]
pub struct XyzTrajectoryWriter {
    /// Whether to include velocities in the output.
    pub include_velocities: bool,
    /// Whether to include forces in the output.
    pub include_forces: bool,
    /// Number of decimal places for coordinates.
    pub coord_decimals: usize,
}
impl XyzTrajectoryWriter {
    /// Create a writer with default settings.
    pub fn new() -> Self {
        Self {
            include_velocities: false,
            include_forces: false,
            coord_decimals: 6,
        }
    }
    /// Write a single XYZ frame to a string.
    pub fn write_frame(&self, frame: &XyzFrame) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "{}", frame.natoms());
        let mut comment = format!("step={} time={:.4}", frame.step, frame.time_ps);
        if let Some(e) = frame.total_energy {
            let _ = write!(comment, " energy={e:.6}");
        }
        if let Some(t) = frame.temperature {
            let _ = write!(comment, " temperature={t:.3}");
        }
        if self.include_velocities && frame.velocities.is_some() {
            let _ = write!(comment, " Properties=species:S:1:pos:R:3:vel:R:3");
        } else if self.include_forces && frame.forces.is_some() {
            let _ = write!(comment, " Properties=species:S:1:pos:R:3:forces:R:3");
        } else {
            let _ = write!(comment, " Properties=species:S:1:pos:R:3");
        }
        let _ = writeln!(out, "{comment}");
        let d = self.coord_decimals;
        for i in 0..frame.natoms() {
            let p = frame.positions[i];
            let elem = &frame.elements[i];
            let mut line = format!(
                "{:<3} {:>12.*} {:>12.*} {:>12.*}",
                elem, d, p[0], d, p[1], d, p[2]
            );
            if self.include_velocities
                && let Some(vels) = &frame.velocities
            {
                let v = vels[i];
                write!(
                    line,
                    " {:>12.*} {:>12.*} {:>12.*}",
                    d, v[0], d, v[1], d, v[2]
                )
                .expect("operation should succeed");
            }
            if self.include_forces
                && let Some(frc) = &frame.forces
            {
                let f = frc[i];
                write!(
                    line,
                    " {:>12.*} {:>12.*} {:>12.*}",
                    d, f[0], d, f[1], d, f[2]
                )
                .expect("operation should succeed");
            }
            let _ = writeln!(out, "{line}");
        }
        out
    }
    /// Write a vector of frames as a concatenated multi-frame XYZ file.
    pub fn write_trajectory(&self, frames: &[XyzFrame]) -> String {
        frames
            .iter()
            .map(|f| self.write_frame(f))
            .collect::<Vec<_>>()
            .join("")
    }
    /// Append a frame to an existing string buffer (useful for streaming).
    pub fn append_frame(&self, buf: &mut String, frame: &XyzFrame) {
        buf.push_str(&self.write_frame(frame));
    }
    /// Create a writer that includes velocities and forces.
    pub fn with_dynamics() -> Self {
        Self {
            include_velocities: true,
            include_forces: true,
            coord_decimals: 6,
        }
    }
}
/// Bond order for a chemical bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondOrder {
    /// Single bond (order 1).
    Single,
    /// Double bond (order 2).
    Double,
    /// Triple bond (order 3).
    Triple,
    /// Aromatic bond (order 1.5).
    Aromatic,
}
impl BondOrder {
    /// Integer bond order (aromatic → 1).
    pub fn as_int(self) -> u8 {
        match self {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 1,
        }
    }
    /// Float bond order (aromatic → 1.5).
    pub fn as_float(self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
        }
    }
}
