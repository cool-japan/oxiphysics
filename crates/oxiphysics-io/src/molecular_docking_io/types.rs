//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use std::fmt::Write as _;
/// A virtual screening results collection.
#[derive(Debug, Clone)]
pub struct VirtualScreeningResults {
    /// List of individual compound results.
    pub results: Vec<VirtualScreeningResult>,
    /// Receptor name used in the screen.
    pub receptor_name: String,
    /// Total compounds screened.
    pub total_screened: u64,
}
impl VirtualScreeningResults {
    /// Create a new `VirtualScreeningResults`.
    pub fn new(receptor_name: &str) -> Self {
        Self {
            results: Vec::new(),
            receptor_name: receptor_name.to_string(),
            total_screened: 0,
        }
    }
    /// Add a result.
    pub fn add_result(&mut self, result: VirtualScreeningResult) {
        self.total_screened += 1;
        self.results.push(result);
    }
    /// Sort results by best docking score (ascending, best = most negative).
    pub fn sort_by_score(&mut self) {
        self.results.sort_by(|a, b| {
            a.best_score
                .partial_cmp(&b.best_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    /// Return the top-N results by docking score.
    pub fn top_n(&self, n: usize) -> &[VirtualScreeningResult] {
        let end = n.min(self.results.len());
        &self.results[..end]
    }
    /// Serialize all results to CSV format.
    pub fn to_csv(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(
            out,
            "compound_id,smiles,best_score_kcal_mol,num_poses,lipinski_ok,mol_weight_da"
        );
        for r in &self.results {
            let _ = writeln!(out, "{}", r.to_csv_line());
        }
        out
    }
    /// Count Lipinski-compliant hits.
    pub fn lipinski_hit_count(&self) -> usize {
        self.results.iter().filter(|r| r.lipinski_ok).count()
    }
}
/// A simplified SMILES bond record.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesBond {
    /// Atom index of the first atom (0-based).
    pub atom_a: usize,
    /// Atom index of the second atom (0-based).
    pub atom_b: usize,
    /// Bond order: 1 = single, 2 = double, 3 = triple, 4 = aromatic.
    pub bond_order: u8,
}
impl SmilesBond {
    /// Create a new `SmilesBond`.
    pub fn new(atom_a: usize, atom_b: usize, bond_order: u8) -> Self {
        Self {
            atom_a,
            atom_b,
            bond_order,
        }
    }
}
/// Cluster of docking poses by RMSD.
#[derive(Debug, Clone)]
pub struct PoseCluster {
    /// Cluster ID (0-based).
    pub cluster_id: usize,
    /// Indices of poses in this cluster (into an external pose list).
    pub pose_indices: Vec<usize>,
    /// RMSD threshold used to define membership.
    pub rmsd_cutoff: f64,
    /// Best (lowest energy) pose index within the cluster.
    pub best_pose_index: usize,
}
impl PoseCluster {
    /// Create a new `PoseCluster`.
    pub fn new(cluster_id: usize, rmsd_cutoff: f64) -> Self {
        Self {
            cluster_id,
            pose_indices: Vec::new(),
            rmsd_cutoff,
            best_pose_index: 0,
        }
    }
    /// Cluster size.
    pub fn size(&self) -> usize {
        self.pose_indices.len()
    }
}
/// A MOL2 molecule.
#[derive(Debug, Clone)]
pub struct Mol2Molecule {
    /// Molecule name.
    pub mol_name: String,
    /// Molecule type (e.g. "SMALL", "PROTEIN").
    pub mol_type: String,
    /// Charge type (e.g. "GASTEIGER").
    pub charge_type: String,
    /// Atoms.
    pub atoms: Vec<Mol2Atom>,
    /// Bonds.
    pub bonds: Vec<Mol2Bond>,
}
impl Mol2Molecule {
    /// Create an empty `Mol2Molecule`.
    pub fn new(mol_name: &str, mol_type: &str, charge_type: &str) -> Self {
        Self {
            mol_name: mol_name.to_string(),
            mol_type: mol_type.to_string(),
            charge_type: charge_type.to_string(),
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }
    /// Serialize to MOL2 format string.
    pub fn to_mol2_string(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "@<TRIPOS>MOLECULE");
        let _ = writeln!(out, "{}", self.mol_name);
        let _ = writeln!(out, "{:4} {:4} 0 0 0", self.atoms.len(), self.bonds.len());
        let _ = writeln!(out, "{}", self.mol_type);
        let _ = writeln!(out, "{}", self.charge_type);
        let _ = writeln!(out);
        let _ = writeln!(out, "@<TRIPOS>ATOM");
        for atom in &self.atoms {
            let _ = writeln!(out, "{}", atom.to_mol2_line());
        }
        let _ = writeln!(out, "@<TRIPOS>BOND");
        for bond in &self.bonds {
            let _ = writeln!(out, "{}", bond.to_mol2_line());
        }
        out
    }
}
/// A docking restraint definition.
#[derive(Debug, Clone)]
pub struct DockingRestraint {
    /// Restraint ID.
    pub id: u32,
    /// Restraint type and parameters.
    pub restraint_type: RestraintType,
    /// Force constant (kcal/mol/Ų).
    pub force_constant: f64,
    /// Whether this restraint is active.
    pub active: bool,
}
impl DockingRestraint {
    /// Create a new `DockingRestraint`.
    pub fn new(id: u32, restraint_type: RestraintType, force_constant: f64) -> Self {
        Self {
            id,
            restraint_type,
            force_constant,
            active: true,
        }
    }
    /// Evaluate the restraint penalty for the given receptor and ligand atoms.
    ///
    /// Returns the energy penalty (kcal/mol).
    pub fn evaluate_penalty(
        &self,
        receptor_atoms: &[PdbqtAtom],
        ligand_atoms: &[PdbqtAtom],
    ) -> f64 {
        if !self.active {
            return 0.0;
        }
        match &self.restraint_type {
            RestraintType::Distance {
                receptor_atom,
                ligand_atom,
                target_distance,
                tolerance,
            } => {
                if *receptor_atom >= receptor_atoms.len() || *ligand_atom >= ligand_atoms.len() {
                    return 0.0;
                }
                let ra = &receptor_atoms[*receptor_atom];
                let la = &ligand_atoms[*ligand_atom];
                let dx = ra.position[0] - la.position[0];
                let dy = ra.position[1] - la.position[1];
                let dz = ra.position[2] - la.position[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let dev = (dist - target_distance).abs() - tolerance;
                if dev > 0.0 {
                    0.5 * self.force_constant * dev * dev
                } else {
                    0.0
                }
            }
            RestraintType::Position {
                ligand_atom,
                reference,
                tolerance,
            } => {
                if *ligand_atom >= ligand_atoms.len() {
                    return 0.0;
                }
                let la = &ligand_atoms[*ligand_atom];
                let dx = la.position[0] - reference[0];
                let dy = la.position[1] - reference[1];
                let dz = la.position[2] - reference[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let dev = (dist - tolerance).max(0.0);
                0.5 * self.force_constant * dev * dev
            }
            RestraintType::ExcludedSphere {
                ligand_atom,
                centre,
                radius,
            } => {
                if *ligand_atom >= ligand_atoms.len() {
                    return 0.0;
                }
                let la = &ligand_atoms[*ligand_atom];
                let dx = la.position[0] - centre[0];
                let dy = la.position[1] - centre[1];
                let dz = la.position[2] - centre[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let overlap = radius - dist;
                if overlap > 0.0 {
                    0.5 * self.force_constant * overlap * overlap
                } else {
                    0.0
                }
            }
        }
    }
}
/// A virtual screening result entry for a single compound.
#[derive(Debug, Clone)]
pub struct VirtualScreeningResult {
    /// Compound identifier.
    pub compound_id: String,
    /// SMILES string.
    pub smiles: String,
    /// Best docking score (kcal/mol).
    pub best_score: f64,
    /// Number of poses generated.
    pub num_poses: u32,
    /// Drug-likeness flags (Lipinski rule of 5: MW, HBD, HBA, logP).
    pub lipinski_ok: bool,
    /// Molecular weight (Da).
    pub molecular_weight: f64,
}
impl VirtualScreeningResult {
    /// Create a new `VirtualScreeningResult`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        compound_id: &str,
        smiles: &str,
        best_score: f64,
        num_poses: u32,
        lipinski_ok: bool,
        molecular_weight: f64,
    ) -> Self {
        Self {
            compound_id: compound_id.to_string(),
            smiles: smiles.to_string(),
            best_score,
            num_poses,
            lipinski_ok,
            molecular_weight,
        }
    }
    /// Format as a CSV line.
    pub fn to_csv_line(&self) -> String {
        format!(
            "{},{},{:.3},{},{},{}",
            self.compound_id,
            self.smiles,
            self.best_score,
            self.num_poses,
            self.lipinski_ok,
            self.molecular_weight,
        )
    }
}
/// A parsed docking run log.
#[derive(Debug, Clone)]
pub struct DockingLog {
    /// Software version string.
    pub version: String,
    /// Receptor file name.
    pub receptor_file: String,
    /// Ligand file name.
    pub ligand_file: String,
    /// Grid box used.
    pub grid_box: Option<GridBox>,
    /// Parsed Vina results.
    pub vina_results: Vec<VinaResult>,
    /// Exhaustiveness setting.
    pub exhaustiveness: u32,
    /// Seed used.
    pub seed: i64,
}
impl DockingLog {
    /// Create an empty `DockingLog`.
    pub fn new() -> Self {
        Self {
            version: String::new(),
            receptor_file: String::new(),
            ligand_file: String::new(),
            grid_box: None,
            vina_results: Vec::new(),
            exhaustiveness: 8,
            seed: 0,
        }
    }
    /// Parse a Vina log string into a `DockingLog`.
    pub fn from_vina_log(log: &str) -> Self {
        let mut dl = Self::new();
        for line in log.lines() {
            if line.starts_with("AutoDock Vina") || line.starts_with("Vina") {
                dl.version = line.trim().to_string();
            } else if line.contains("receptor") && dl.receptor_file.is_empty() {
                if let Some(part) = line.split_whitespace().last() {
                    dl.receptor_file = part.to_string();
                }
            } else if line.contains("ligand") && dl.ligand_file.is_empty() {
                if let Some(part) = line.split_whitespace().last() {
                    dl.ligand_file = part.to_string();
                }
            } else if line.contains("exhaustiveness") {
                if let Some(part) = line.split_whitespace().last()
                    && let Ok(e) = part.parse::<u32>()
                {
                    dl.exhaustiveness = e;
                }
            } else if line.contains("seed")
                && let Some(part) = line.split_whitespace().last()
                && let Ok(s) = part.parse::<i64>()
            {
                dl.seed = s;
            }
        }
        dl.vina_results = parse_vina_log(log);
        dl
    }
    /// Best affinity from all results (or 0.0 if none).
    pub fn best_affinity(&self) -> f64 {
        self.vina_results
            .iter()
            .map(|r| r.affinity)
            .fold(f64::INFINITY, f64::min)
    }
}
/// A connection table atom entry in SDF/MOL format.
#[derive(Debug, Clone)]
pub struct SdfAtom {
    /// Cartesian position (Angstrom).
    pub position: [f64; 3],
    /// Element symbol.
    pub element: String,
    /// Formal charge code (0 = no charge, 3 = +1, 5 = -1, ...).
    pub charge_code: u8,
}
impl SdfAtom {
    /// Create a new `SdfAtom`.
    pub fn new(position: [f64; 3], element: &str, charge_code: u8) -> Self {
        Self {
            position,
            element: element.to_string(),
            charge_code,
        }
    }
    /// Format as a V2000 atom line.
    pub fn to_v2000_line(&self) -> String {
        format!(
            "{:10.4}{:10.4}{:10.4} {:<3}  0  {}  0  0  0  0  0  0  0  0  0  0",
            self.position[0], self.position[1], self.position[2], self.element, self.charge_code,
        )
    }
}
/// A bond table entry in SDF/MOL format.
#[derive(Debug, Clone)]
pub struct SdfBond {
    /// First atom index (1-based).
    pub atom1: u32,
    /// Second atom index (1-based).
    pub atom2: u32,
    /// Bond type (1 = single, 2 = double, 3 = triple, 4 = aromatic).
    pub bond_type: u8,
    /// Bond stereo (0 = not stereo, 1 = up, 6 = down).
    pub stereo: u8,
}
impl SdfBond {
    /// Create a new `SdfBond`.
    pub fn new(atom1: u32, atom2: u32, bond_type: u8, stereo: u8) -> Self {
        Self {
            atom1,
            atom2,
            bond_type,
            stereo,
        }
    }
    /// Format as a V2000 bond line.
    pub fn to_v2000_line(&self) -> String {
        format!(
            "{:3}{:3}{:3}{:3}  0  0  0",
            self.atom1, self.atom2, self.bond_type, self.stereo
        )
    }
}
/// A bond in Tripos MOL2 format.
#[derive(Debug, Clone)]
pub struct Mol2Bond {
    /// Bond ID (1-based).
    pub bond_id: u32,
    /// Origin atom ID (1-based).
    pub origin_atom_id: u32,
    /// Target atom ID (1-based).
    pub target_atom_id: u32,
    /// Bond type: "1", "2", "3", "ar", "am", "nc".
    pub bond_type: String,
}
impl Mol2Bond {
    /// Create a new `Mol2Bond`.
    pub fn new(bond_id: u32, origin: u32, target: u32, bond_type: &str) -> Self {
        Self {
            bond_id,
            origin_atom_id: origin,
            target_atom_id: target,
            bond_type: bond_type.to_string(),
        }
    }
    /// Format as a MOL2 `@`TRIPOS`BOND` record.
    pub fn to_mol2_line(&self) -> String {
        format!(
            "{:6}  {:6}  {:6}  {}",
            self.bond_id, self.origin_atom_id, self.target_atom_id, self.bond_type
        )
    }
}
/// SDF (Structure Data File) molecule entry.
#[derive(Debug, Clone)]
pub struct SdfEntry {
    /// Molecule name.
    pub name: String,
    /// Program / timestamp line.
    pub info_line: String,
    /// Comment line.
    pub comment: String,
    /// Atoms.
    pub atoms: Vec<SdfAtom>,
    /// Bonds.
    pub bonds: Vec<SdfBond>,
    /// SD tags (key-value pairs).
    pub tags: Vec<(String, String)>,
}
impl SdfEntry {
    /// Create an empty `SdfEntry`.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            info_line: "  OxiPhysics        3D".to_string(),
            comment: String::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            tags: Vec::new(),
        }
    }
    /// Add an SD data tag.
    pub fn add_tag(&mut self, key: &str, value: &str) {
        self.tags.push((key.to_string(), value.to_string()));
    }
    /// Serialize to SDF V2000 format.
    pub fn to_sdf_string(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "{}", self.name);
        let _ = writeln!(out, "{}", self.info_line);
        let _ = writeln!(out, "{}", self.comment);
        let _ = writeln!(
            out,
            "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000",
            self.atoms.len(),
            self.bonds.len()
        );
        for atom in &self.atoms {
            let _ = writeln!(out, "{}", atom.to_v2000_line());
        }
        for bond in &self.bonds {
            let _ = writeln!(out, "{}", bond.to_v2000_line());
        }
        let _ = writeln!(out, "M  END");
        for (key, val) in &self.tags {
            let _ = writeln!(out, ">  <{key}>");
            let _ = writeln!(out, "{val}");
            let _ = writeln!(out);
        }
        let _ = writeln!(out, "$$$$");
        out
    }
}
/// SMILES molecule representation.
#[derive(Debug, Clone)]
pub struct SmilesMolecule {
    /// SMILES string.
    pub smiles: String,
    /// Parsed atoms.
    pub atoms: Vec<SmilesAtom>,
    /// Parsed bonds.
    pub bonds: Vec<SmilesBond>,
    /// Molecule name.
    pub name: String,
}
impl SmilesMolecule {
    /// Create a new `SmilesMolecule` from a SMILES string (lightweight tokenizer).
    ///
    /// This is a simplified tokenizer for demonstration; a complete SMILES parser
    /// requires a full grammar implementation.
    pub fn from_smiles(smiles: &str, name: &str) -> Self {
        let atoms = tokenize_smiles_atoms(smiles);
        Self {
            smiles: smiles.to_string(),
            atoms,
            bonds: Vec::new(),
            name: name.to_string(),
        }
    }
    /// Generate a canonical SMILES string from stored atoms.
    ///
    /// Returns the original SMILES string (round-trip is identity for simplified parser).
    pub fn to_smiles(&self) -> String {
        self.smiles.clone()
    }
    /// Molecular formula string (e.g. "C6H6").
    pub fn molecular_formula(&self) -> String {
        let mut counts = std::collections::HashMap::new();
        for atom in &self.atoms {
            *counts.entry(atom.element.clone()).or_insert(0u32) += 1;
        }
        let mut parts: Vec<String> = counts
            .iter()
            .map(|(el, &cnt)| {
                if cnt == 1 {
                    el.clone()
                } else {
                    format!("{el}{cnt}")
                }
            })
            .collect();
        parts.sort();
        parts.join("")
    }
}
/// Type of pharmacophoric feature.
#[derive(Debug, Clone, PartialEq)]
pub enum PharmacophoreFeatureType {
    /// Hydrogen bond donor.
    HbDonor,
    /// Hydrogen bond acceptor.
    HbAcceptor,
    /// Hydrophobic region.
    Hydrophobic,
    /// Aromatic ring.
    Aromatic,
    /// Positively ionizable group.
    PositiveIonizable,
    /// Negatively ionizable group.
    NegativeIonizable,
    /// Excluded volume sphere.
    ExcludedVolume,
}
/// PDBQT molecule (receptor or ligand).
#[derive(Debug, Clone)]
pub struct PdbqtMolecule {
    /// Molecule name / title.
    pub name: String,
    /// List of atoms.
    pub atoms: Vec<PdbqtAtom>,
    /// Remarks from the file header.
    pub remarks: Vec<String>,
    /// Number of rotatable bonds (from PDBQT TORSDOF record).
    pub torsional_degrees_of_freedom: u32,
}
impl PdbqtMolecule {
    /// Create an empty `PdbqtMolecule`.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            atoms: Vec::new(),
            remarks: Vec::new(),
            torsional_degrees_of_freedom: 0,
        }
    }
    /// Add an atom to the molecule.
    pub fn add_atom(&mut self, atom: PdbqtAtom) {
        self.atoms.push(atom);
    }
    /// Parse a PDBQT-formatted string into a `PdbqtMolecule`.
    ///
    /// Supports `ATOM`, `HETATM`, `REMARK`, `TORSDOF` records.
    pub fn from_pdbqt_str(name: &str, content: &str) -> Self {
        let mut mol = Self::new(name);
        for line in content.lines() {
            let record = &line[..line.len().min(6)];
            match record.trim() {
                "ATOM" | "HETATM" => {
                    if let Some(atom) = parse_pdbqt_atom_line(line) {
                        mol.atoms.push(atom);
                    }
                }
                "REMARK" => {
                    mol.remarks.push(line[6..].trim().to_string());
                }
                "TORSDOF" => {
                    if let Ok(n) = line[7..].trim().parse::<u32>() {
                        mol.torsional_degrees_of_freedom = n;
                    }
                }
                _ => {}
            }
        }
        mol
    }
    /// Serialize this molecule to a PDBQT-formatted string.
    pub fn to_pdbqt_string(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "REMARK  Name = {}", self.name);
        for remark in &self.remarks {
            let _ = writeln!(out, "REMARK  {remark}");
        }
        for atom in &self.atoms {
            let _ = writeln!(out, "{}", atom.to_pdbqt_line());
        }
        let _ = writeln!(out, "TORSDOF {}", self.torsional_degrees_of_freedom);
        let _ = writeln!(out, "END");
        out
    }
    /// Centre of mass of the molecule (uniform weighting).
    pub fn centre_of_mass(&self) -> [f64; 3] {
        if self.atoms.is_empty() {
            return [0.0; 3];
        }
        let n = self.atoms.len() as f64;
        let sum = self.atoms.iter().fold([0.0_f64; 3], |acc, a| {
            [
                acc[0] + a.position[0],
                acc[1] + a.position[1],
                acc[2] + a.position[2],
            ]
        });
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }
}
/// A docking pose, representing a ligand conformation with associated scores.
#[derive(Debug, Clone)]
pub struct DockingPose {
    /// Pose index (1-based).
    pub pose_index: u32,
    /// Binding affinity score (kcal/mol, lower = better).
    pub binding_affinity: f64,
    /// RMSD from best mode (Angstrom).
    pub rmsd_lb: f64,
    /// RMSD from best mode upper bound (Angstrom).
    pub rmsd_ub: f64,
    /// Ligand atoms in this pose.
    pub atoms: Vec<PdbqtAtom>,
    /// Cluster membership index (0 = unassigned).
    pub cluster_id: usize,
}
impl DockingPose {
    /// Create a new `DockingPose`.
    pub fn new(pose_index: u32, binding_affinity: f64, rmsd_lb: f64, rmsd_ub: f64) -> Self {
        Self {
            pose_index,
            binding_affinity,
            rmsd_lb,
            rmsd_ub,
            atoms: Vec::new(),
            cluster_id: 0,
        }
    }
    /// Add an atom to this pose.
    pub fn add_atom(&mut self, atom: PdbqtAtom) {
        self.atoms.push(atom);
    }
    /// Compute RMSD between this pose and another pose.
    ///
    /// Both poses must have the same number of atoms.
    /// Returns `None` if atom counts differ.
    pub fn rmsd_to(&self, other: &DockingPose) -> Option<f64> {
        if self.atoms.len() != other.atoms.len() {
            return None;
        }
        let n = self.atoms.len() as f64;
        if n < 1e-10 {
            return Some(0.0);
        }
        let sum: f64 = self
            .atoms
            .iter()
            .zip(other.atoms.iter())
            .map(|(a, b)| {
                let dx = a.position[0] - b.position[0];
                let dy = a.position[1] - b.position[1];
                let dz = a.position[2] - b.position[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        Some((sum / n).sqrt())
    }
}
/// A simplified SMILES atom entry.
#[derive(Debug, Clone, PartialEq)]
pub struct SmilesAtom {
    /// Elemental symbol.
    pub element: String,
    /// Formal charge.
    pub formal_charge: i32,
    /// Isotope mass number (0 = unspecified).
    pub isotope: u32,
    /// Whether the atom is aromatic.
    pub aromatic: bool,
    /// Number of implicit hydrogens.
    pub implicit_h: u32,
}
impl SmilesAtom {
    /// Create a new `SmilesAtom`.
    pub fn new(element: &str, formal_charge: i32, aromatic: bool) -> Self {
        Self {
            element: element.to_string(),
            formal_charge,
            isotope: 0,
            aromatic,
            implicit_h: 0,
        }
    }
}
/// A single atom in Tripos MOL2 format.
#[derive(Debug, Clone)]
pub struct Mol2Atom {
    /// Atom ID (1-based).
    pub atom_id: u32,
    /// Atom name.
    pub atom_name: String,
    /// Cartesian coordinates (Angstrom).
    pub position: [f64; 3],
    /// Tripos atom type (e.g. "C.3", "N.am").
    pub atom_type: String,
    /// Substructure ID.
    pub subst_id: u32,
    /// Substructure name.
    pub subst_name: String,
    /// Partial charge.
    pub charge: f64,
}
impl Mol2Atom {
    /// Create a new `Mol2Atom`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        atom_id: u32,
        atom_name: &str,
        position: [f64; 3],
        atom_type: &str,
        subst_id: u32,
        subst_name: &str,
        charge: f64,
    ) -> Self {
        Self {
            atom_id,
            atom_name: atom_name.to_string(),
            position,
            atom_type: atom_type.to_string(),
            subst_id,
            subst_name: subst_name.to_string(),
            charge,
        }
    }
    /// Format this atom as a MOL2 `@`TRIPOS`ATOM` record.
    pub fn to_mol2_line(&self) -> String {
        format!(
            "{:6}  {:<8}  {:9.4}  {:9.4}  {:9.4}  {:<8}  {:4}  {:<8}  {:9.4}",
            self.atom_id,
            self.atom_name,
            self.position[0],
            self.position[1],
            self.position[2],
            self.atom_type,
            self.subst_id,
            self.subst_name,
            self.charge,
        )
    }
}
/// Docking grid box specification for AutoDock Vina / AutoGrid.
#[derive(Debug, Clone)]
pub struct GridBox {
    /// Center of the grid box (Angstrom).
    pub center: [f64; 3],
    /// Size of the grid box in each dimension (Angstrom).
    pub size: [f64; 3],
    /// Grid spacing (Angstrom).
    pub spacing: f64,
}
impl GridBox {
    /// Create a new `GridBox`.
    pub fn new(center: [f64; 3], size: [f64; 3], spacing: f64) -> Self {
        Self {
            center,
            size,
            spacing: spacing.max(1e-6),
        }
    }
    /// Number of grid points in each dimension.
    pub fn grid_points(&self) -> [u32; 3] {
        [
            ((self.size[0] / self.spacing).ceil() as u32).max(1),
            ((self.size[1] / self.spacing).ceil() as u32).max(1),
            ((self.size[2] / self.spacing).ceil() as u32).max(1),
        ]
    }
    /// Total number of grid points.
    pub fn total_grid_points(&self) -> u64 {
        let gp = self.grid_points();
        gp[0] as u64 * gp[1] as u64 * gp[2] as u64
    }
    /// Check whether a position lies within the grid box.
    pub fn contains(&self, pos: [f64; 3]) -> bool {
        for i in 0..3 {
            let half = self.size[i] * 0.5;
            if (pos[i] - self.center[i]).abs() > half {
                return false;
            }
        }
        true
    }
    /// Write as an AutoGrid gpf-style configuration.
    pub fn to_gpf_string(&self, receptor_name: &str, ligand_types: &[&str]) -> String {
        let mut out = String::new();
        let _ = writeln!(
            out,
            "npts {} {} {}",
            self.grid_points()[0],
            self.grid_points()[1],
            self.grid_points()[2]
        );
        let _ = writeln!(out, "spacing {:.3}", self.spacing);
        let _ = writeln!(
            out,
            "gridcenter {:.3} {:.3} {:.3}",
            self.center[0], self.center[1], self.center[2]
        );
        let _ = writeln!(out, "receptor {receptor_name}");
        for lt in ligand_types {
            let _ = writeln!(out, "ligand_types {lt}");
        }
        out
    }
}
/// A pharmacophore hypothesis (collection of features).
#[derive(Debug, Clone)]
pub struct PharmacophoreHypothesis {
    /// Hypothesis name.
    pub name: String,
    /// Pharmacophore features.
    pub features: Vec<PharmacophoreFeature>,
}
impl PharmacophoreHypothesis {
    /// Create a new `PharmacophoreHypothesis`.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            features: Vec::new(),
        }
    }
    /// Add a feature.
    pub fn add_feature(&mut self, feature: PharmacophoreFeature) {
        self.features.push(feature);
    }
    /// Export to a simple pharmacophore file format (JSON-like text).
    pub fn to_export_string(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "PHARMACOPHORE {}", self.name);
        for (i, f) in self.features.iter().enumerate() {
            let _ = writeln!(
                out,
                "FEATURE {:3} {:4} {:8.3} {:8.3} {:8.3} {:6.2} {:6.3}",
                i + 1,
                f.feature_type_str(),
                f.position[0],
                f.position[1],
                f.position[2],
                f.tolerance,
                f.weight,
            );
        }
        out
    }
    /// Count matched features for a set of query positions.
    ///
    /// Returns number of features that are satisfied by at least one query point.
    pub fn count_matches(&self, query_positions: &[[f64; 3]]) -> usize {
        self.features
            .iter()
            .filter(|f| query_positions.iter().any(|&q| f.is_matched(q)))
            .count()
    }
}
/// A single atom record in PDBQT format.
#[derive(Debug, Clone, PartialEq)]
pub struct PdbqtAtom {
    /// Serial number.
    pub serial: u32,
    /// Atom name (e.g. "CA").
    pub name: String,
    /// Residue name (e.g. "ALA").
    pub res_name: String,
    /// Chain ID.
    pub chain_id: char,
    /// Residue sequence number.
    pub res_seq: i32,
    /// Cartesian position (Angstrom).
    pub position: [f64; 3],
    /// Occupancy factor.
    pub occupancy: f64,
    /// Temperature factor (B-factor).
    pub temp_factor: f64,
    /// Partial charge (AutoDock specific).
    pub charge: f64,
    /// AutoDock atom type (e.g. "C", "OA", "HD").
    pub atom_type: String,
}
impl PdbqtAtom {
    /// Create a new `PdbqtAtom` with all fields specified.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        serial: u32,
        name: &str,
        res_name: &str,
        chain_id: char,
        res_seq: i32,
        position: [f64; 3],
        occupancy: f64,
        temp_factor: f64,
        charge: f64,
        atom_type: &str,
    ) -> Self {
        Self {
            serial,
            name: name.to_string(),
            res_name: res_name.to_string(),
            chain_id,
            res_seq,
            position,
            occupancy,
            temp_factor,
            charge,
            atom_type: atom_type.to_string(),
        }
    }
    /// Format this atom as a PDBQT ATOM record line.
    pub fn to_pdbqt_line(&self) -> String {
        format!(
            "ATOM  {:5} {:<4} {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}    {:+7.3} {:<2}",
            self.serial,
            self.name,
            self.res_name,
            self.chain_id,
            self.res_seq,
            self.position[0],
            self.position[1],
            self.position[2],
            self.occupancy,
            self.temp_factor,
            self.charge,
            self.atom_type,
        )
    }
}
/// A single result entry parsed from Vina's log output.
#[derive(Debug, Clone)]
pub struct VinaResult {
    /// Mode number (1-based).
    pub mode: u32,
    /// Affinity (kcal/mol).
    pub affinity: f64,
    /// RMSD lower bound (Angstrom).
    pub rmsd_lb: f64,
    /// RMSD upper bound (Angstrom).
    pub rmsd_ub: f64,
}
impl VinaResult {
    /// Create a new `VinaResult`.
    pub fn new(mode: u32, affinity: f64, rmsd_lb: f64, rmsd_ub: f64) -> Self {
        Self {
            mode,
            affinity,
            rmsd_lb,
            rmsd_ub,
        }
    }
}
/// A receptor-ligand complex combining receptor and ligand PDBQT data.
#[derive(Debug, Clone)]
pub struct ReceptorLigandComplex {
    /// Receptor molecule.
    pub receptor: PdbqtMolecule,
    /// Ligand molecule (best docking pose).
    pub ligand: PdbqtMolecule,
    /// Best docking pose result.
    pub best_pose: Option<DockingPose>,
    /// Binding affinity score.
    pub score: Option<BindingAffinityScore>,
}
impl ReceptorLigandComplex {
    /// Create a new `ReceptorLigandComplex`.
    pub fn new(receptor: PdbqtMolecule, ligand: PdbqtMolecule) -> Self {
        Self {
            receptor,
            ligand,
            best_pose: None,
            score: None,
        }
    }
    /// Serialize the complex as a single PDBQT file string.
    pub fn to_complex_pdbqt(&self) -> String {
        let mut out = String::new();
        let _ = writeln!(out, "REMARK  RECEPTOR-LIGAND COMPLEX");
        let _ = writeln!(out, "REMARK  Receptor: {}", self.receptor.name);
        let _ = writeln!(out, "REMARK  Ligand:   {}", self.ligand.name);
        if let Some(ref s) = self.score {
            let _ = writeln!(out, "REMARK  Score: {:.2} kcal/mol", s.total);
        }
        out.push_str(&self.receptor.to_pdbqt_string());
        out.push_str("TER\n");
        out.push_str(&self.ligand.to_pdbqt_string());
        out
    }
    /// Number of intermolecular contacts (atom pairs within `cutoff` Angstrom).
    pub fn contact_count(&self, cutoff: f64) -> usize {
        let cutoff2 = cutoff * cutoff;
        let mut count = 0usize;
        for ra in &self.receptor.atoms {
            for la in &self.ligand.atoms {
                let dx = ra.position[0] - la.position[0];
                let dy = ra.position[1] - la.position[1];
                let dz = ra.position[2] - la.position[2];
                if dx * dx + dy * dy + dz * dz <= cutoff2 {
                    count += 1;
                }
            }
        }
        count
    }
}
/// Type of docking restraint.
#[derive(Debug, Clone, PartialEq)]
pub enum RestraintType {
    /// Distance restraint between two atom indices.
    Distance {
        /// Receptor atom index.
        receptor_atom: usize,
        /// Ligand atom index.
        ligand_atom: usize,
        /// Target distance (Angstrom).
        target_distance: f64,
        /// Tolerance (Angstrom).
        tolerance: f64,
    },
    /// Position restraint: keep ligand near a reference position.
    Position {
        /// Ligand atom index.
        ligand_atom: usize,
        /// Reference position (Angstrom).
        reference: [f64; 3],
        /// Tolerance radius (Angstrom).
        tolerance: f64,
    },
    /// Excluded volume: ligand atom must stay outside a sphere.
    ExcludedSphere {
        /// Ligand atom index.
        ligand_atom: usize,
        /// Sphere centre (Angstrom).
        centre: [f64; 3],
        /// Sphere radius (Angstrom).
        radius: f64,
    },
}
/// A pharmacophore feature point.
#[derive(Debug, Clone)]
pub struct PharmacophoreFeature {
    /// Feature type.
    pub feature_type: PharmacophoreFeatureType,
    /// Centroid position (Angstrom).
    pub position: [f64; 3],
    /// Tolerance radius (Angstrom).
    pub tolerance: f64,
    /// Feature weight for scoring.
    pub weight: f64,
    /// Optional direction vector for directional features.
    pub direction: Option<[f64; 3]>,
}
impl PharmacophoreFeature {
    /// Create a new `PharmacophoreFeature`.
    pub fn new(
        feature_type: PharmacophoreFeatureType,
        position: [f64; 3],
        tolerance: f64,
        weight: f64,
    ) -> Self {
        Self {
            feature_type,
            position,
            tolerance,
            weight,
            direction: None,
        }
    }
    /// Check whether a given position satisfies this pharmacophore constraint.
    pub fn is_matched(&self, query_pos: [f64; 3]) -> bool {
        let dx = query_pos[0] - self.position[0];
        let dy = query_pos[1] - self.position[1];
        let dz = query_pos[2] - self.position[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        dist <= self.tolerance
    }
    /// Feature type as a string identifier.
    pub fn feature_type_str(&self) -> &'static str {
        match self.feature_type {
            PharmacophoreFeatureType::HbDonor => "HBD",
            PharmacophoreFeatureType::HbAcceptor => "HBA",
            PharmacophoreFeatureType::Hydrophobic => "HYD",
            PharmacophoreFeatureType::Aromatic => "AR",
            PharmacophoreFeatureType::PositiveIonizable => "PI",
            PharmacophoreFeatureType::NegativeIonizable => "NI",
            PharmacophoreFeatureType::ExcludedVolume => "EV",
        }
    }
}
/// Binding affinity score record with breakdown of energy components.
#[derive(Debug, Clone)]
pub struct BindingAffinityScore {
    /// Total estimated free energy of binding (kcal/mol).
    pub total: f64,
    /// Intermolecular energy (kcal/mol).
    pub intermolecular: f64,
    /// Internal energy of the ligand (kcal/mol).
    pub internal_ligand: f64,
    /// Torsional entropy penalty (kcal/mol).
    pub torsional: f64,
    /// Unbound ligand internal energy (kcal/mol).
    pub unbound_ligand: f64,
}
impl BindingAffinityScore {
    /// Create a new `BindingAffinityScore`.
    pub fn new(
        total: f64,
        intermolecular: f64,
        internal_ligand: f64,
        torsional: f64,
        unbound_ligand: f64,
    ) -> Self {
        Self {
            total,
            intermolecular,
            internal_ligand,
            torsional,
            unbound_ligand,
        }
    }
    /// Estimated binding free energy: ΔG = intermolecular + torsional.
    pub fn estimated_delta_g(&self) -> f64 {
        self.intermolecular + self.torsional
    }
    /// Format the score as a PDBQT REMARK block.
    pub fn to_remark_block(&self) -> String {
        format!(
            "REMARK VINA RESULT:    {:.1}      {:.3}      {:.3}",
            self.total, 0.0_f64, 0.0_f64
        )
    }
}

/// Compute the pairwise RMSD lower and upper bounds across all pose pairs.
///
/// Iterates over all (i, j) pairs with i < j, computes RMSD for each pair
/// (skipping pairs with different atom counts), and returns `(min_rmsd, max_rmsd)`.
/// Returns `(0.0, 0.0)` when fewer than two poses are supplied.
pub fn compute_pose_cluster_rmsd_bounds(poses: &[DockingPose]) -> (f64, f64) {
    if poses.len() < 2 {
        return (0.0, 0.0);
    }
    let mut lower = f64::INFINITY;
    let mut upper = f64::NEG_INFINITY;
    for i in 0..poses.len() {
        for j in (i + 1)..poses.len() {
            if let Some(d) = poses[i].rmsd_to(&poses[j]) {
                if d < lower {
                    lower = d;
                }
                if d > upper {
                    upper = d;
                }
            }
        }
    }
    if lower.is_infinite() {
        return (0.0, 0.0);
    }
    (lower, upper)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_pose(positions: &[[f64; 3]]) -> DockingPose {
        let mut pose = DockingPose::new(1, -7.0, 0.0, 0.0);
        for (i, &pos) in positions.iter().enumerate() {
            pose.add_atom(PdbqtAtom::new(
                i as u32 + 1,
                "C",
                "LIG",
                'A',
                1,
                pos,
                1.0,
                0.0,
                0.0,
                "C",
            ));
        }
        pose
    }

    #[test]
    fn test_rmsd_bounds_empty() {
        let (lo, hi) = compute_pose_cluster_rmsd_bounds(&[]);
        assert_eq!(lo, 0.0);
        assert_eq!(hi, 0.0);
    }

    #[test]
    fn test_rmsd_bounds_single_pose() {
        let poses = vec![make_pose(&[[0.0, 0.0, 0.0]])];
        let (lo, hi) = compute_pose_cluster_rmsd_bounds(&poses);
        assert_eq!(lo, 0.0);
        assert_eq!(hi, 0.0);
    }

    #[test]
    fn test_rmsd_bounds_two_poses_known_distance() {
        // 3 atoms shifted by (1,0,0) → RMSD = sqrt(3/3) = 1.0
        let poses = vec![
            make_pose(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            make_pose(&[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        ];
        let (lo, hi) = compute_pose_cluster_rmsd_bounds(&poses);
        assert!(
            (lo - 1.0).abs() < 1e-10,
            "lower bound should be 1.0, got {lo}"
        );
        assert!(
            (hi - 1.0).abs() < 1e-10,
            "upper bound should be 1.0, got {hi}"
        );
    }

    #[test]
    fn test_rmsd_bounds_three_poses() {
        // d(0,1)=1, d(0,2)=2, d(1,2)=1 → lb=1, ub=2
        let poses = vec![
            make_pose(&[[0.0, 0.0, 0.0]]),
            make_pose(&[[1.0, 0.0, 0.0]]),
            make_pose(&[[2.0, 0.0, 0.0]]),
        ];
        let (lo, hi) = compute_pose_cluster_rmsd_bounds(&poses);
        assert!(
            (lo - 1.0).abs() < 1e-10,
            "lower bound should be 1.0, got {lo}"
        );
        assert!(
            (hi - 2.0).abs() < 1e-10,
            "upper bound should be 2.0, got {hi}"
        );
    }
}
