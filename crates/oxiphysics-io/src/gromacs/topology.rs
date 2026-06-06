// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GROMACS topology, force field, index, MDP, and run-input types.

// ─── Topology types ─────────────────────────────────────────────────────────

/// Topology section types in a GROMACS `.top` file.
#[derive(Debug, Clone, PartialEq)]
pub enum TopologySection {
    /// `[ atoms ]` section.
    Atoms,
    /// `[ bonds ]` section.
    Bonds,
    /// `[ angles ]` section.
    Angles,
    /// `[ dihedrals ]` section.
    Dihedrals,
    /// `[ pairs ]` section.
    Pairs,
    /// Any other section.
    Other(String),
}

impl TopologySection {
    /// Parse a section name string into a `TopologySection`.
    pub(crate) fn from_name(name: &str) -> Self {
        match name.trim().to_lowercase().as_str() {
            "atoms" => Self::Atoms,
            "bonds" => Self::Bonds,
            "angles" => Self::Angles,
            "dihedrals" => Self::Dihedrals,
            "pairs" => Self::Pairs,
            other => Self::Other(other.to_string()),
        }
    }
}

/// Parsed representation of a GROMACS topology (`.top`) file.
#[derive(Debug, Clone)]
pub struct TopFile {
    /// Sections with their content lines.
    pub sections: Vec<(TopologySection, Vec<String>)>,
}

impl TopFile {
    /// Parse a topology file from a string.
    pub fn parse(data: &str) -> Result<Self, String> {
        let mut sections: Vec<(TopologySection, Vec<String>)> = Vec::new();
        let mut current_section: Option<TopologySection> = None;
        let mut current_lines: Vec<String> = Vec::new();

        for line in data.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with(';') {
                continue;
            }
            if trimmed.starts_with('[') && trimmed.ends_with(']') {
                if let Some(sec) = current_section.take() {
                    sections.push((sec, std::mem::take(&mut current_lines)));
                }
                let name = trimmed.trim_start_matches('[').trim_end_matches(']').trim();
                current_section = Some(TopologySection::from_name(name));
            } else if current_section.is_some() {
                current_lines.push(trimmed.to_string());
            }
        }

        if let Some(sec) = current_section {
            sections.push((sec, current_lines));
        }

        Ok(TopFile { sections })
    }

    /// Count atoms in the `[ atoms ]` section.
    pub fn atom_count(&self) -> usize {
        self.sections
            .iter()
            .filter(|(sec, _)| *sec == TopologySection::Atoms)
            .map(|(_, lines)| lines.len())
            .sum()
    }

    /// Serialize the topology back to a string.
    pub fn to_top_string(&self) -> String {
        let mut s = String::new();
        for (sec, lines) in &self.sections {
            let name = match sec {
                TopologySection::Atoms => "atoms",
                TopologySection::Bonds => "bonds",
                TopologySection::Angles => "angles",
                TopologySection::Dihedrals => "dihedrals",
                TopologySection::Pairs => "pairs",
                TopologySection::Other(n) => n.as_str(),
            };
            s.push_str(&format!("[ {} ]\n", name));
            for line in lines {
                s.push_str(line);
                s.push('\n');
            }
            s.push('\n');
        }
        s
    }
}

// ─── Force Field Parameters ─────────────────────────────────────────────────

impl std::str::FromStr for TopFile {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

/// Lennard-Jones parameters for an atom type.
#[derive(Debug, Clone, PartialEq)]
pub struct LJParams {
    /// Atom type name.
    pub atom_type: String,
    /// sigma parameter (nm).
    pub sigma: f64,
    /// epsilon parameter (kJ/mol).
    pub epsilon: f64,
}

/// Bond parameters (harmonic).
#[derive(Debug, Clone, PartialEq)]
pub struct BondParams {
    /// Atom type i.
    pub type_i: String,
    /// Atom type j.
    pub type_j: String,
    /// Equilibrium bond length (nm).
    pub r0: f64,
    /// Force constant (kJ/mol/nm^2).
    pub k: f64,
}

/// Angle parameters (harmonic).
#[derive(Debug, Clone, PartialEq)]
pub struct AngleParams {
    /// Atom type i.
    pub type_i: String,
    /// Atom type j (center).
    pub type_j: String,
    /// Atom type k.
    pub type_k: String,
    /// Equilibrium angle (degrees).
    pub theta0: f64,
    /// Force constant (kJ/mol/rad^2).
    pub k: f64,
}

/// Collection of force field parameters.
#[derive(Debug, Clone, Default)]
pub struct ForceFieldParams {
    /// LJ parameters per atom type.
    pub lj_params: Vec<LJParams>,
    /// Bond parameters.
    pub bond_params: Vec<BondParams>,
    /// Angle parameters.
    pub angle_params: Vec<AngleParams>,
}

impl ForceFieldParams {
    /// Create empty force field parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add LJ parameters for an atom type.
    pub fn add_lj(&mut self, atom_type: &str, sigma: f64, epsilon: f64) {
        self.lj_params.push(LJParams {
            atom_type: atom_type.to_string(),
            sigma,
            epsilon,
        });
    }

    /// Add bond parameters.
    pub fn add_bond(&mut self, type_i: &str, type_j: &str, r0: f64, k: f64) {
        self.bond_params.push(BondParams {
            type_i: type_i.to_string(),
            type_j: type_j.to_string(),
            r0,
            k,
        });
    }

    /// Add angle parameters.
    pub fn add_angle(&mut self, type_i: &str, type_j: &str, type_k: &str, theta0: f64, k: f64) {
        self.angle_params.push(AngleParams {
            type_i: type_i.to_string(),
            type_j: type_j.to_string(),
            type_k: type_k.to_string(),
            theta0,
            k,
        });
    }

    /// Look up LJ parameters by atom type.
    pub fn get_lj(&self, atom_type: &str) -> Option<&LJParams> {
        self.lj_params.iter().find(|p| p.atom_type == atom_type)
    }

    /// Look up bond parameters by pair of atom types.
    pub fn get_bond(&self, type_i: &str, type_j: &str) -> Option<&BondParams> {
        self.bond_params.iter().find(|p| {
            (p.type_i == type_i && p.type_j == type_j) || (p.type_i == type_j && p.type_j == type_i)
        })
    }
}

// ─── Index Group Handling ───────────────────────────────────────────────────

/// A named group of atom indices (like GROMACS index groups).
#[derive(Debug, Clone, PartialEq)]
pub struct IndexGroup {
    /// Group name.
    pub name: String,
    /// Atom indices (0-based).
    pub indices: Vec<usize>,
}

/// Parsed index file (`.ndx`).
#[derive(Debug, Clone, Default)]
pub struct IndexFile {
    /// Named groups.
    pub groups: Vec<IndexGroup>,
}

impl IndexFile {
    /// Create an empty index file.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse an index file from a string.
    pub fn parse(data: &str) -> Result<Self, String> {
        let mut groups: Vec<IndexGroup> = Vec::new();
        let mut current_name: Option<String> = None;
        let mut current_indices: Vec<usize> = Vec::new();

        for line in data.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('[') && trimmed.ends_with(']') {
                if let Some(name) = current_name.take() {
                    groups.push(IndexGroup {
                        name,
                        indices: std::mem::take(&mut current_indices),
                    });
                }
                let name = trimmed
                    .trim_start_matches('[')
                    .trim_end_matches(']')
                    .trim()
                    .to_string();
                current_name = Some(name);
            } else if current_name.is_some() {
                for token in trimmed.split_whitespace() {
                    if let Ok(idx) = token.parse::<usize>() {
                        current_indices.push(idx.saturating_sub(1));
                    }
                }
            }
        }

        if let Some(name) = current_name {
            groups.push(IndexGroup {
                name,
                indices: current_indices,
            });
        }

        Ok(IndexFile { groups })
    }

    /// Serialize the index file back to a string.
    pub fn to_ndx_string(&self) -> String {
        let mut s = String::new();
        for group in &self.groups {
            s.push_str(&format!("[ {} ]\n", group.name));
            for (i, &idx) in group.indices.iter().enumerate() {
                if i > 0 && i % 15 == 0 {
                    s.push('\n');
                }
                s.push_str(&format!("{} ", idx + 1));
            }
            s.push('\n');
        }
        s
    }

    /// Get a group by name.
    pub fn get_group(&self, name: &str) -> Option<&IndexGroup> {
        self.groups.iter().find(|g| g.name == name)
    }

    /// Number of groups.
    pub fn group_count(&self) -> usize {
        self.groups.len()
    }
}

impl std::str::FromStr for IndexFile {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ─── MDP File Parsing ───────────────────────────────────────────────────────

/// Parsed MDP (molecular dynamics parameters) file.
#[derive(Debug, Clone, Default)]
pub struct MdpFile {
    /// Key-value parameters.
    pub params: Vec<(String, String)>,
}

impl MdpFile {
    /// Create an empty MDP file.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse an MDP file from a string.
    pub fn parse(data: &str) -> Result<Self, String> {
        let mut params = Vec::new();
        for line in data.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with(';') {
                continue;
            }
            if let Some(eq_pos) = trimmed.find('=') {
                let key = trimmed[..eq_pos].trim().to_string();
                let value = trimmed[eq_pos + 1..].trim().to_string();
                let value = if let Some(sc) = value.find(';') {
                    value[..sc].trim().to_string()
                } else {
                    value
                };
                params.push((key, value));
            }
        }
        Ok(MdpFile { params })
    }

    /// Get a parameter value by key.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.params
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// Get a parameter as f64.
    pub fn get_f64(&self, key: &str) -> Option<f64> {
        self.get(key).and_then(|v| v.parse().ok())
    }

    /// Get a parameter as i64.
    pub fn get_i64(&self, key: &str) -> Option<i64> {
        self.get(key).and_then(|v| v.parse().ok())
    }

    /// Serialize the MDP file to a string.
    pub fn to_mdp_string(&self) -> String {
        let mut s = String::new();
        for (k, v) in &self.params {
            s.push_str(&format!("{:<24} = {}\n", k, v));
        }
        s
    }

    /// Set a parameter (overwrite if exists, append if not).
    pub fn set(&mut self, key: &str, value: &str) {
        if let Some(entry) = self.params.iter_mut().find(|(k, _)| k == key) {
            entry.1 = value.to_string();
        } else {
            self.params.push((key.to_string(), value.to_string()));
        }
    }

    /// Number of parameters.
    pub fn param_count(&self) -> usize {
        self.params.len()
    }
}

impl std::str::FromStr for MdpFile {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ─── GROMACS Run Input Generation ───────────────────────────────────────────

/// Generates GROMACS-compatible run input files.
pub struct GromacsRunInputGenerator;

impl GromacsRunInputGenerator {
    /// Generate a minimal MDP file for energy minimization.
    pub fn energy_minimization_mdp() -> MdpFile {
        let mut mdp = MdpFile::new();
        mdp.set("integrator", "steep");
        mdp.set("nsteps", "5000");
        mdp.set("emtol", "100.0");
        mdp.set("emstep", "0.01");
        mdp.set("nstxout", "100");
        mdp.set("nstenergy", "100");
        mdp.set("cutoff-scheme", "Verlet");
        mdp.set("coulombtype", "PME");
        mdp.set("rcoulomb", "1.0");
        mdp.set("rvdw", "1.0");
        mdp
    }

    /// Generate a minimal MDP file for NVT equilibration.
    pub fn nvt_equilibration_mdp(nsteps: i64, dt: f64, temperature: f64) -> MdpFile {
        let mut mdp = MdpFile::new();
        mdp.set("integrator", "md");
        mdp.set("nsteps", &nsteps.to_string());
        mdp.set("dt", &dt.to_string());
        mdp.set("nstxout", "500");
        mdp.set("nstvout", "500");
        mdp.set("nstenergy", "500");
        mdp.set("cutoff-scheme", "Verlet");
        mdp.set("coulombtype", "PME");
        mdp.set("rcoulomb", "1.0");
        mdp.set("rvdw", "1.0");
        mdp.set("tcoupl", "V-rescale");
        mdp.set("ref-t", &temperature.to_string());
        mdp.set("tau-t", "0.1");
        mdp.set("tc-grps", "System");
        mdp.set("gen-vel", "yes");
        mdp.set("gen-temp", &temperature.to_string());
        mdp
    }

    /// Generate a minimal MDP file for NPT production.
    pub fn npt_production_mdp(nsteps: i64, dt: f64, temperature: f64, pressure: f64) -> MdpFile {
        let mut mdp = Self::nvt_equilibration_mdp(nsteps, dt, temperature);
        mdp.set("pcoupl", "Parrinello-Rahman");
        mdp.set("ref-p", &pressure.to_string());
        mdp.set("tau-p", "2.0");
        mdp.set("compressibility", "4.5e-5");
        mdp.set("gen-vel", "no");
        mdp
    }
}

// ============================================================================
// Resolved Topology
// ============================================================================

/// Result of resolving a topology with `#include` directives.
#[derive(Debug, Clone, Default)]
pub struct ResolvedTopology {
    /// Merged section list after #include expansion.
    pub sections: Vec<(TopologySection, Vec<String>)>,
    /// List of file names that were #included.
    pub included_files: Vec<String>,
}

impl ResolvedTopology {
    /// Parse a topology string, recording `#include` directives.
    pub fn parse(s: &str) -> Self {
        let mut sections: Vec<(TopologySection, Vec<String>)> = Vec::new();
        let mut included_files: Vec<String> = Vec::new();
        let mut current_section: Option<(TopologySection, Vec<String>)> = None;

        for line in s.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with(';') || trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with("#include") {
                if let Some(fname) = trimmed
                    .split_whitespace()
                    .nth(1)
                    .map(|s| s.trim_matches('"').to_string())
                {
                    included_files.push(fname);
                }
                continue;
            }
            if trimmed.starts_with('[') && trimmed.ends_with(']') {
                if let Some(sec) = current_section.take() {
                    sections.push(sec);
                }
                let name = trimmed[1..trimmed.len() - 1].trim().to_string();
                let sec_type = match name.as_str() {
                    "atoms" => TopologySection::Atoms,
                    "bonds" => TopologySection::Bonds,
                    "angles" => TopologySection::Angles,
                    _ => TopologySection::Other(name),
                };
                current_section = Some((sec_type, Vec::new()));
            } else if let Some(ref mut sec) = current_section {
                sec.1.push(trimmed.to_string());
            }
        }
        if let Some(sec) = current_section {
            sections.push(sec);
        }
        ResolvedTopology {
            sections,
            included_files,
        }
    }

    /// Count the number of sections of a given type.
    pub fn count_section(&self, t: &TopologySection) -> usize {
        self.sections.iter().filter(|(s, _)| s == t).count()
    }
}

impl From<&str> for ResolvedTopology {
    fn from(s: &str) -> Self {
        Self::parse(s)
    }
}

// ============================================================================
// Dihedral / Extended Force Field
// ============================================================================

/// A dihedral (torsion) parameter entry.
#[derive(Debug, Clone, PartialEq)]
pub struct DihedralParams {
    /// Atom type i.
    pub type_i: String,
    /// Atom type j.
    pub type_j: String,
    /// Atom type k.
    pub type_k: String,
    /// Atom type l.
    pub type_l: String,
    /// Phase angle (degrees).
    pub phi0: f64,
    /// Force constant (kJ/mol).
    pub k: f64,
    /// Multiplicity n.
    pub n: u32,
}

/// Extended force field parameter set with dihedrals.
#[derive(Debug, Clone, Default)]
pub struct ExtForceFieldParams {
    /// Base parameters (LJ + bonds + angles).
    pub base: ForceFieldParams,
    /// Dihedral parameters.
    pub dihedrals: Vec<DihedralParams>,
}

impl ExtForceFieldParams {
    /// Add a dihedral parameter.
    pub fn add_dihedral(
        &mut self,
        type_i: &str,
        type_j: &str,
        type_k: &str,
        type_l: &str,
        phi0: f64,
        k: f64,
        n: u32,
    ) {
        self.dihedrals.push(DihedralParams {
            type_i: type_i.to_string(),
            type_j: type_j.to_string(),
            type_k: type_k.to_string(),
            type_l: type_l.to_string(),
            phi0,
            k,
            n,
        });
    }

    /// Compute OPLS-AA style torsional energy for a dihedral angle `phi` (rad).
    pub fn torsion_energy(
        &self,
        type_i: &str,
        type_j: &str,
        type_k: &str,
        type_l: &str,
        phi: f64,
    ) -> f64 {
        self.dihedrals
            .iter()
            .filter(|d| {
                d.type_i == type_i && d.type_j == type_j && d.type_k == type_k && d.type_l == type_l
            })
            .map(|d| {
                let phi0_rad = d.phi0 * std::f64::consts::PI / 180.0;
                d.k * (1.0 + (d.n as f64 * phi - phi0_rad).cos())
            })
            .sum()
    }
}

// ============================================================================
// TopologyData
// ============================================================================

/// A parsed atom entry from a GROMACS `[ atoms ]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct TopAtom {
    /// Atom index (1-based in GRO convention).
    pub index: usize,
    /// Atom type string.
    pub atom_type: String,
    /// Residue index.
    pub residue_index: usize,
    /// Residue name.
    pub residue_name: String,
    /// Atom name within the residue.
    pub atom_name: String,
    /// Charge group number.
    pub charge_group: i32,
    /// Partial charge (e).
    pub charge: f64,
    /// Mass (amu).
    pub mass: f64,
}

/// A bond entry from a `[ bonds ]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct TopBond {
    /// First atom index (1-based).
    pub atom_i: usize,
    /// Second atom index (1-based).
    pub atom_j: usize,
    /// Function type (1 = harmonic).
    pub func_type: i32,
}

/// An angle entry from a `[ angles ]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct TopAngle {
    /// First atom index (1-based).
    pub atom_i: usize,
    /// Second (central) atom index (1-based).
    pub atom_j: usize,
    /// Third atom index (1-based).
    pub atom_k: usize,
    /// Function type (1 = harmonic).
    pub func_type: i32,
}

/// A dihedral entry from a `[ dihedrals ]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct TopDihedral {
    /// First atom index (1-based).
    pub atom_i: usize,
    /// Second atom index (1-based).
    pub atom_j: usize,
    /// Third atom index (1-based).
    pub atom_k: usize,
    /// Fourth atom index (1-based).
    pub atom_l: usize,
    /// Function type.
    pub func_type: i32,
}

/// A pair entry from a `[ pairs ]` section.
#[derive(Debug, Clone, PartialEq)]
pub struct TopPair {
    /// First atom index (1-based).
    pub atom_i: usize,
    /// Second atom index (1-based).
    pub atom_j: usize,
    /// Function type (1 = LJ pairs).
    pub func_type: i32,
}

/// A `[ moleculetype ]` entry.
#[derive(Debug, Clone, PartialEq)]
pub struct TopMoleculeType {
    /// Molecule name.
    pub name: String,
    /// Number of excluded neighbours.
    pub nrexcl: i32,
}

/// Rich parsed representation of a GROMACS topology.
#[derive(Debug, Clone, Default)]
pub struct TopologyData {
    /// Molecule type records.
    pub molecule_types: Vec<TopMoleculeType>,
    /// Atom records.
    pub atoms: Vec<TopAtom>,
    /// Bond records.
    pub bonds: Vec<TopBond>,
    /// Pair records.
    pub pairs: Vec<TopPair>,
    /// Angle records.
    pub angles: Vec<TopAngle>,
    /// Dihedral records.
    pub dihedrals: Vec<TopDihedral>,
}

impl TopologyData {
    /// Create empty topology data.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse topology data from a GROMACS `.top` string.
    pub fn parse(s: &str) -> Result<Self, String> {
        let mut data = TopologyData::new();
        let mut current_section: Option<TopologySection> = None;

        for (line_no, raw_line) in s.lines().enumerate() {
            let trimmed = raw_line.trim();
            if trimmed.is_empty() || trimmed.starts_with(';') {
                continue;
            }
            let trimmed = if let Some(pos) = trimmed.find(';') {
                trimmed[..pos].trim()
            } else {
                trimmed
            };
            if trimmed.is_empty() {
                continue;
            }

            if trimmed.starts_with('[') && trimmed.ends_with(']') {
                let name = trimmed[1..trimmed.len() - 1].trim();
                current_section = Some(TopologySection::from_name(name));
                continue;
            }

            match &current_section {
                Some(TopologySection::Other(n)) if n == "moleculetype" => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 2 {
                        let nrexcl: i32 = parts[1]
                            .parse()
                            .map_err(|e| format!("line {}: bad nrexcl: {e}", line_no + 1))?;
                        data.molecule_types.push(TopMoleculeType {
                            name: parts[0].to_string(),
                            nrexcl,
                        });
                    }
                }
                Some(TopologySection::Atoms) => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 7 {
                        let index: usize = parts[0]
                            .parse()
                            .map_err(|e| format!("line {}: bad atom index: {e}", line_no + 1))?;
                        let atom_type = parts[1].to_string();
                        let residue_index: usize = parts[2]
                            .parse()
                            .map_err(|e| format!("line {}: bad resnr: {e}", line_no + 1))?;
                        let residue_name = parts[3].to_string();
                        let atom_name = parts[4].to_string();
                        let charge_group: i32 = parts[5]
                            .parse()
                            .map_err(|e| format!("line {}: bad cgnr: {e}", line_no + 1))?;
                        let charge: f64 = parts[6]
                            .parse()
                            .map_err(|e| format!("line {}: bad charge: {e}", line_no + 1))?;
                        let mass: f64 = if parts.len() >= 8 {
                            parts[7].parse().unwrap_or(0.0)
                        } else {
                            0.0
                        };
                        data.atoms.push(TopAtom {
                            index,
                            atom_type,
                            residue_index,
                            residue_name,
                            atom_name,
                            charge_group,
                            charge,
                            mass,
                        });
                    }
                }
                Some(TopologySection::Bonds) => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 2 {
                        let i: usize = parts[0]
                            .parse()
                            .map_err(|e| format!("line {}: bad bond atom_i: {e}", line_no + 1))?;
                        let j: usize = parts[1]
                            .parse()
                            .map_err(|e| format!("line {}: bad bond atom_j: {e}", line_no + 1))?;
                        let func: i32 = parts.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
                        data.bonds.push(TopBond {
                            atom_i: i,
                            atom_j: j,
                            func_type: func,
                        });
                    }
                }
                Some(TopologySection::Pairs) => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 2 {
                        let i: usize = parts[0].parse().unwrap_or(0);
                        let j: usize = parts[1].parse().unwrap_or(0);
                        let func: i32 = parts.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
                        data.pairs.push(TopPair {
                            atom_i: i,
                            atom_j: j,
                            func_type: func,
                        });
                    }
                }
                Some(TopologySection::Angles) => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 3 {
                        let i: usize = parts[0].parse().unwrap_or(0);
                        let j: usize = parts[1].parse().unwrap_or(0);
                        let k: usize = parts[2].parse().unwrap_or(0);
                        let func: i32 = parts.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
                        data.angles.push(TopAngle {
                            atom_i: i,
                            atom_j: j,
                            atom_k: k,
                            func_type: func,
                        });
                    }
                }
                Some(TopologySection::Dihedrals) => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 4 {
                        let i: usize = parts[0].parse().unwrap_or(0);
                        let j: usize = parts[1].parse().unwrap_or(0);
                        let k: usize = parts[2].parse().unwrap_or(0);
                        let l: usize = parts[3].parse().unwrap_or(0);
                        let func: i32 = parts.get(4).and_then(|s| s.parse().ok()).unwrap_or(1);
                        data.dihedrals.push(TopDihedral {
                            atom_i: i,
                            atom_j: j,
                            atom_k: k,
                            atom_l: l,
                            func_type: func,
                        });
                    }
                }
                _ => {}
            }
        }
        Ok(data)
    }

    /// Serialize back to a GROMACS-style topology string.
    pub fn to_top_string(&self) -> String {
        let mut s = String::new();

        if !self.molecule_types.is_empty() {
            s.push_str("[ moleculetype ]\n; name  nrexcl\n");
            for mol in &self.molecule_types {
                s.push_str(&format!("{:<8} {}\n", mol.name, mol.nrexcl));
            }
            s.push('\n');
        }
        if !self.atoms.is_empty() {
            s.push_str("[ atoms ]\n; nr type resnr residue atom cgnr charge mass\n");
            for a in &self.atoms {
                s.push_str(&format!(
                    "{:5} {:6} {:5} {:6} {:6} {:5} {:9.5} {:9.4}\n",
                    a.index,
                    a.atom_type,
                    a.residue_index,
                    a.residue_name,
                    a.atom_name,
                    a.charge_group,
                    a.charge,
                    a.mass,
                ));
            }
            s.push('\n');
        }
        if !self.bonds.is_empty() {
            s.push_str("[ bonds ]\n; ai aj funct\n");
            for b in &self.bonds {
                s.push_str(&format!("{} {} {}\n", b.atom_i, b.atom_j, b.func_type));
            }
            s.push('\n');
        }
        if !self.pairs.is_empty() {
            s.push_str("[ pairs ]\n; ai aj funct\n");
            for p in &self.pairs {
                s.push_str(&format!("{} {} {}\n", p.atom_i, p.atom_j, p.func_type));
            }
            s.push('\n');
        }
        if !self.angles.is_empty() {
            s.push_str("[ angles ]\n; ai aj ak funct\n");
            for a in &self.angles {
                s.push_str(&format!(
                    "{} {} {} {}\n",
                    a.atom_i, a.atom_j, a.atom_k, a.func_type
                ));
            }
            s.push('\n');
        }
        if !self.dihedrals.is_empty() {
            s.push_str("[ dihedrals ]\n; ai aj ak al funct\n");
            for d in &self.dihedrals {
                s.push_str(&format!(
                    "{} {} {} {} {}\n",
                    d.atom_i, d.atom_j, d.atom_k, d.atom_l, d.func_type
                ));
            }
            s.push('\n');
        }
        s
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }
    /// Number of bonds.
    pub fn n_bonds(&self) -> usize {
        self.bonds.len()
    }
    /// Number of angles.
    pub fn n_angles(&self) -> usize {
        self.angles.len()
    }
    /// Number of dihedrals.
    pub fn n_dihedrals(&self) -> usize {
        self.dihedrals.len()
    }
}

impl std::str::FromStr for TopologyData {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ============================================================================
// TopWriter / DihedralEntry
// ============================================================================

/// A writer helper for GROMACS topology dihedral sections.
pub struct TopWriter;

/// A single dihedral entry for [`TopWriter::write_dihedrals`].
#[derive(Debug, Clone)]
pub struct DihedralEntry {
    /// First atom index (1-based).
    pub atom_i: usize,
    /// Second atom index (1-based).
    pub atom_j: usize,
    /// Third atom index (1-based).
    pub atom_k: usize,
    /// Fourth atom index (1-based).
    pub atom_l: usize,
    /// Function type.
    pub func_type: i32,
    /// Equilibrium dihedral angle (degrees).
    pub phi0_deg: f64,
    /// Force constant (kJ mol^-1 rad^-2).
    pub k_phi: f64,
    /// Multiplicity.
    pub multiplicity: i32,
}

impl TopWriter {
    /// Generate the `[ dihedrals ]` section text from a slice of entries.
    pub fn write_dihedrals(dihedrals: &[DihedralEntry]) -> String {
        let mut s = String::new();
        s.push_str("[ dihedrals ]\n");
        s.push_str("; ai  aj  ak  al  funct  phi0  kphi  mult\n");
        for d in dihedrals {
            s.push_str(&format!(
                "{:4} {:4} {:4} {:4} {:5}  {:8.3}  {:8.4}  {:4}\n",
                d.atom_i,
                d.atom_j,
                d.atom_k,
                d.atom_l,
                d.func_type,
                d.phi0_deg,
                d.k_phi,
                d.multiplicity
            ));
        }
        s
    }

    /// Parse the `[ dihedrals ]` section text back into `DihedralEntry` records.
    pub fn read_dihedrals(text: &str) -> std::result::Result<Vec<DihedralEntry>, String> {
        let mut entries = Vec::new();
        for line in text.lines() {
            let t = line.trim();
            if t.is_empty() || t.starts_with(';') || t.starts_with('[') {
                continue;
            }
            let parts: Vec<&str> = t.split_whitespace().collect();
            if parts.len() < 5 {
                continue;
            }
            let atom_i: usize = parts[0].parse().map_err(|e| format!("bad atom_i: {e}"))?;
            let atom_j: usize = parts[1].parse().map_err(|e| format!("bad atom_j: {e}"))?;
            let atom_k: usize = parts[2].parse().map_err(|e| format!("bad atom_k: {e}"))?;
            let atom_l: usize = parts[3].parse().map_err(|e| format!("bad atom_l: {e}"))?;
            let func_type: i32 = parts[4].parse().map_err(|e| format!("bad funct: {e}"))?;
            let phi0_deg: f64 = if parts.len() > 5 {
                parts[5].parse().unwrap_or(0.0)
            } else {
                0.0
            };
            let k_phi: f64 = if parts.len() > 6 {
                parts[6].parse().unwrap_or(0.0)
            } else {
                0.0
            };
            let multiplicity: i32 = if parts.len() > 7 {
                parts[7].parse().unwrap_or(0)
            } else {
                0
            };
            entries.push(DihedralEntry {
                atom_i,
                atom_j,
                atom_k,
                atom_l,
                func_type,
                phi0_deg,
                k_phi,
                multiplicity,
            });
        }
        Ok(entries)
    }
}

// ============================================================================
// GroTopBuilder / MoleculeEntry / GroAtomType / LJ helpers
// ============================================================================

/// Represents a molecule type entry in the `[ molecules ]` directive.
#[derive(Debug, Clone)]
pub struct MoleculeEntry {
    /// Molecule type name.
    pub mol_name: String,
    /// Number of instances.
    pub count: usize,
}

/// Simplified GROMACS topology builder.
#[derive(Debug, Clone, Default)]
pub struct GroTopBuilder {
    /// Title of the system.
    pub system_name: String,
    /// Molecule type entries.
    pub molecules: Vec<MoleculeEntry>,
    /// Atom types defined in `[ atomtypes ]`.
    pub atom_types: Vec<GroAtomType>,
}

/// One atom type record in `[ atomtypes ]`.
#[derive(Debug, Clone)]
pub struct GroAtomType {
    /// Atom type name.
    pub name: String,
    /// Atomic mass in u.
    pub mass: f64,
    /// Partial charge in e.
    pub charge: f64,
    /// LJ epsilon parameter (kJ/mol).
    pub epsilon: f64,
    /// LJ sigma parameter (nm).
    pub sigma: f64,
}

impl GroTopBuilder {
    /// Create an empty topology builder.
    pub fn new(system_name: &str) -> Self {
        GroTopBuilder {
            system_name: system_name.to_string(),
            molecules: Vec::new(),
            atom_types: Vec::new(),
        }
    }

    /// Add a molecule type entry.
    pub fn add_molecule(&mut self, name: &str, count: usize) {
        self.molecules.push(MoleculeEntry {
            mol_name: name.to_string(),
            count,
        });
    }

    /// Add an atom type definition.
    pub fn add_atom_type(&mut self, at: GroAtomType) {
        self.atom_types.push(at);
    }

    /// Total number of molecules.
    pub fn total_molecules(&self) -> usize {
        self.molecules.iter().map(|m| m.count).sum()
    }

    /// Write the minimal `[ system ]` and `[ molecules ]` directives.
    pub fn write_system<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writeln!(writer, "[ system ]")?;
        writeln!(writer, "; Name")?;
        writeln!(writer, "{}", self.system_name)?;
        writeln!(writer)?;
        writeln!(writer, "[ molecules ]")?;
        writeln!(writer, "; Compound       #mols")?;
        for m in &self.molecules {
            writeln!(writer, "{:<16}{}", m.mol_name, m.count)?;
        }
        Ok(())
    }

    /// Write the `[ atomtypes ]` directive.
    pub fn write_atom_types<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writeln!(writer, "[ atomtypes ]")?;
        writeln!(writer, "; name   mass    charge    sigma    epsilon")?;
        for at in &self.atom_types {
            writeln!(
                writer,
                "{:<8} {:>10.4} {:>8.4} {:>12.6} {:>12.6}",
                at.name, at.mass, at.charge, at.sigma, at.epsilon
            )?;
        }
        Ok(())
    }

    /// Find an atom type by name.
    pub fn find_atom_type(&self, name: &str) -> Option<&GroAtomType> {
        self.atom_types.iter().find(|at| at.name == name)
    }
}

// ─── LJ force-field parameter helpers ──────────────────────────────────────

/// Convert a GROMACS LJ (C6, C12) representation to (epsilon, sigma).
pub fn c6c12_to_epsilon_sigma(c6: f64, c12: f64) -> (f64, f64) {
    if c6 <= 0.0 || c12 <= 0.0 {
        return (0.0, 0.0);
    }
    let sigma = (c12 / c6).powf(1.0 / 6.0);
    let epsilon = c6 * c6 / (4.0 * c12);
    (epsilon, sigma)
}

/// Convert (epsilon, sigma) to GROMACS (C6, C12).
pub fn epsilon_sigma_to_c6c12(epsilon: f64, sigma: f64) -> (f64, f64) {
    let s6 = sigma.powi(6);
    let c6 = 4.0 * epsilon * s6;
    let c12 = 4.0 * epsilon * s6 * s6;
    (c6, c12)
}

/// Compute Lorentz-Berthelot combining rules for two atom types.
pub fn lorentz_berthelot(eps_i: f64, sig_i: f64, eps_j: f64, sig_j: f64) -> (f64, f64) {
    let sigma_ij = (sig_i + sig_j) / 2.0;
    let epsilon_ij = (eps_i * eps_j).sqrt();
    (epsilon_ij, sigma_ij)
}

/// Evaluate the Lennard-Jones 12-6 potential at distance `r`.
pub fn lj_potential(r: f64, epsilon: f64, sigma: f64) -> f64 {
    if r <= 0.0 {
        return f64::INFINITY;
    }
    let sr6 = (sigma / r).powi(6);
    4.0 * epsilon * (sr6 * sr6 - sr6)
}

/// Evaluate the LJ force magnitude (dU/dr) at distance `r`.
pub fn lj_force(r: f64, epsilon: f64, sigma: f64) -> f64 {
    if r <= 0.0 {
        return f64::INFINITY;
    }
    let sr = sigma / r;
    let sr6 = sr.powi(6);
    24.0 * epsilon / r * (2.0 * sr6 * sr6 - sr6)
}

/// Parse a single `[ atomtypes ]` line into a `GroAtomType`.
pub fn parse_atomtype_line(line: &str) -> Option<GroAtomType> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 5 {
        return None;
    }
    let name = parts[0].to_string();
    let mass = parts[1].parse::<f64>().ok()?;
    let charge = parts[2].parse::<f64>().ok()?;
    let sigma = parts[3].parse::<f64>().ok()?;
    let epsilon = parts[4].parse::<f64>().ok()?;
    Some(GroAtomType {
        name,
        mass,
        charge,
        epsilon,
        sigma,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests_topology {
    use super::*;

    const WATER_TOP: &str = r#"
[ moleculetype ]
; name nrexcl
SOL 3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1 OW 1 SOL OW 1 -0.834 15.9994
2 HW1 1 SOL HW1 1  0.417  1.0080
3 HW2 1 SOL HW2 1  0.417  1.0080

[ bonds ]
1 2 1
1 3 1

[ pairs ]
2 3 1

[ angles ]
2 1 3 1

[ dihedrals ]
"#;

    #[test]
    fn topology_data_atoms_parsed() {
        let td = TopologyData::parse(WATER_TOP).expect("parse");
        assert_eq!(td.n_atoms(), 3);
        assert_eq!(td.atoms[0].atom_name, "OW");
    }

    #[test]
    fn topology_data_bonds_parsed() {
        let td = TopologyData::parse(WATER_TOP).expect("parse");
        assert_eq!(td.n_bonds(), 2);
    }

    #[test]
    fn topology_data_roundtrip() {
        let td = TopologyData::parse(WATER_TOP).expect("parse");
        let s = td.to_top_string();
        let td2 = TopologyData::parse(&s).expect("reparse");
        assert_eq!(td2.n_atoms(), td.n_atoms());
        assert_eq!(td2.n_bonds(), td.n_bonds());
    }

    #[test]
    fn resolved_topology_include_detection() {
        let top_str = "#include \"amber99.itp\"\n#include \"spc.itp\"\n[ moleculetype ]\nWater 3\n";
        let rt = ResolvedTopology::parse(top_str);
        assert_eq!(rt.included_files.len(), 2);
    }

    #[test]
    fn dihedral_params_torsion_energy() {
        let mut ff = ExtForceFieldParams::default();
        ff.add_dihedral("CT", "CT", "CT", "CT", 0.0, 5.0, 1);
        let e = ff.torsion_energy("CT", "CT", "CT", "CT", 0.0);
        assert!((e - 10.0).abs() < 1e-8);
    }

    #[test]
    fn top_builder_total_molecules() {
        let mut tb = GroTopBuilder::new("Water box");
        tb.add_molecule("SOL", 1000);
        tb.add_molecule("NA", 10);
        assert_eq!(tb.total_molecules(), 1010);
    }

    #[test]
    fn c6c12_roundtrip() {
        let eps = 0.65;
        let sig = 0.315;
        let (c6, c12) = epsilon_sigma_to_c6c12(eps, sig);
        let (eps2, sig2) = c6c12_to_epsilon_sigma(c6, c12);
        assert!((eps2 - eps).abs() < 1e-10);
        assert!((sig2 - sig).abs() < 1e-10);
    }

    #[test]
    fn lj_potential_at_sigma_is_zero() {
        let u = lj_potential(1.0, 1.0, 1.0);
        assert!(u.abs() < 1e-12);
    }

    #[test]
    fn parse_atomtype_line_valid() {
        let line = "OW    15.999  -0.834  0.315890  0.636386";
        let at = parse_atomtype_line(line).expect("parse");
        assert_eq!(at.name, "OW");
    }

    #[test]
    fn parse_atomtype_line_too_short_returns_none() {
        assert!(parse_atomtype_line("OW 15.999 -0.834").is_none());
    }

    #[test]
    fn mdp_file_from_str_basic() {
        let data = "integrator = md\nnsteps = 50000\n";
        let mdp = MdpFile::parse(data).expect("parse");
        assert_eq!(mdp.get("integrator"), Some("md"));
    }

    #[test]
    fn index_file_from_str_basic() {
        let data = "[ System ]\n1 2 3 4\n[ Protein ]\n1 2\n";
        let ndx = IndexFile::parse(data).expect("parse");
        assert!(ndx.get_group("System").is_some());
        assert_eq!(ndx.group_count(), 2);
    }

    #[test]
    fn top_file_atom_count() {
        let data =
            "[ atoms ]\n1 opls_135 1 ETH C1 1 -0.18 12.011\n2 opls_140 1 ETH H1 1  0.06  1.008\n";
        let top = TopFile::parse(data).expect("parse");
        assert_eq!(top.atom_count(), 2);
    }

    #[test]
    fn top_writer_roundtrip() {
        let entries = vec![DihedralEntry {
            atom_i: 1,
            atom_j: 2,
            atom_k: 3,
            atom_l: 4,
            func_type: 1,
            phi0_deg: 180.0,
            k_phi: 4.184,
            multiplicity: 2,
        }];
        let text = TopWriter::write_dihedrals(&entries);
        let parsed = TopWriter::read_dihedrals(&text).expect("parse");
        assert_eq!(parsed.len(), 1);
    }

    #[test]
    fn ff_params_add_and_get_lj() {
        let mut ff = ForceFieldParams::new();
        ff.add_lj("CT", 0.35, 0.276);
        assert!(ff.get_lj("CT").is_some());
    }

    #[test]
    fn em_mdp_has_integrator() {
        let mdp = GromacsRunInputGenerator::energy_minimization_mdp();
        assert!(mdp.get("integrator").is_some());
    }
}
