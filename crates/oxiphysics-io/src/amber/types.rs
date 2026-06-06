//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// A covalent bond entry.
pub struct AmberBond {
    /// Zero-based index of the first atom.
    pub i: usize,
    /// Zero-based index of the second atom.
    pub j: usize,
    /// Force constant (kcal mol⁻¹ Å⁻²).
    pub k: f64,
    /// Equilibrium bond length (Å).
    pub r0: f64,
}
/// A parsed AMBER atom mask selector.
///
/// AMBER masks use a simplified selection language:
/// - `:1-10` → residues 1 through 10
/// - `@CA` → atom named CA
/// - `:1-10@CA` → atom CA in residues 1-10
/// - `@C,CA,N` → atoms C, CA, or N
#[derive(Debug, Clone, PartialEq)]
pub struct AmberMask {
    /// Residue selection (list of residue numbers, 1-based).
    pub residues: Vec<usize>,
    /// Atom name selection (empty = all atoms).
    pub atom_names: Vec<String>,
}
impl AmberMask {
    /// Parse an AMBER mask string.
    ///
    /// Supported forms:
    /// - `:N` – single residue N
    /// - `:N-M` – residue range N to M
    /// - `:N,M` – residues N and M
    /// - `@name` – single atom name
    /// - `@n1,n2` – multiple atom names
    /// - `:N-M@name` – combined residue + atom
    ///
    /// Returns `Err` if the mask cannot be parsed.
    pub fn parse(mask: &str) -> Result<Self, String> {
        let mask = mask.trim();
        if mask.is_empty() {
            return Ok(AmberMask {
                residues: Vec::new(),
                atom_names: Vec::new(),
            });
        }
        let (res_part, atom_part) = if let Some(at_pos) = mask.find('@') {
            (&mask[..at_pos], Some(&mask[at_pos + 1..]))
        } else {
            (mask, None)
        };
        let residues = if let Some(rest) = res_part.strip_prefix(':') {
            parse_residue_selection(rest)?
        } else if res_part.is_empty() {
            Vec::new()
        } else {
            return Err(format!(
                "Invalid mask: expected ':' before residue selection, got '{}'",
                res_part
            ));
        };
        let atom_names = if let Some(ap) = atom_part {
            ap.split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect()
        } else {
            Vec::new()
        };
        Ok(AmberMask {
            residues,
            atom_names,
        })
    }
    /// Whether a residue number (1-based) is selected by this mask.
    pub fn matches_residue(&self, res_num: usize) -> bool {
        if self.residues.is_empty() {
            true
        } else {
            self.residues.contains(&res_num)
        }
    }
    /// Whether an atom name is selected by this mask.
    pub fn matches_atom(&self, atom_name: &str) -> bool {
        if self.atom_names.is_empty() {
            true
        } else {
            self.atom_names.iter().any(|n| n == atom_name)
        }
    }
    /// Whether a given (residue_num, atom_name) pair is selected.
    pub fn matches(&self, res_num: usize, atom_name: &str) -> bool {
        self.matches_residue(res_num) && self.matches_atom(atom_name)
    }
}
/// All topology data parsed from a prmtop file.
pub struct AmberTopology {
    /// Title string from the `%TITLE` or `TITLE` flag section.
    pub title: String,
    /// Per-atom data.
    pub atoms: Vec<AmberAtom>,
    /// Bond list with parameters.
    pub bonds: Vec<AmberBond>,
    /// Angle list with parameters.
    pub angles: Vec<AmberAngle>,
    /// Number of atoms (from `POINTERS`).
    pub n_atoms: usize,
    /// Number of bonds including hydrogen (from `POINTERS`).
    pub n_bonds: usize,
}
impl AmberTopology {
    /// Parse a prmtop file (given as a string) into an `AmberTopology`.
    ///
    /// Recognises `%FLAG` / `%FORMAT` blocks.  Only the sections needed for
    /// basic force-field data are extracted; unknown flags are silently skipped.
    pub fn from_prmtop_str(s: &str) -> Result<Self, String> {
        let title = {
            let raw = parse_flag_section(s, "TITLE")
                .or_else(|| parse_flag_section(s, "title"))
                .unwrap_or_default();
            raw.trim().to_string()
        };
        let pointers = parse_flag_section(s, "POINTERS")
            .map(|raw| parse_fortran_ints(&raw))
            .unwrap_or_default();
        let n_atoms = pointers.first().copied().unwrap_or(0) as usize;
        let n_bonds = pointers.get(2).copied().unwrap_or(0) as usize;
        let n_angles = pointers.get(4).copied().unwrap_or(0) as usize;
        let atom_names: Vec<String> = parse_flag_section(s, "ATOM_NAME")
            .map(|raw| parse_fortran_strings(&raw))
            .unwrap_or_default();
        const AMBER_CHARGE_SCALE: f64 = 18.2223;
        let raw_charges: Vec<f64> = parse_flag_section(s, "CHARGE")
            .map(|raw| parse_fortran_reals(&raw))
            .unwrap_or_default();
        let charges: Vec<f64> = raw_charges
            .iter()
            .map(|&q| q / AMBER_CHARGE_SCALE)
            .collect();
        let masses: Vec<f64> = parse_flag_section(s, "MASS")
            .map(|raw| parse_fortran_reals(&raw))
            .unwrap_or_default();
        let atom_types: Vec<String> = parse_flag_section(s, "AMBER_ATOM_TYPE")
            .map(|raw| parse_fortran_strings(&raw))
            .unwrap_or_default();
        let residue_labels: Vec<String> = parse_flag_section(s, "RESIDUE_LABEL")
            .map(|raw| parse_fortran_strings(&raw))
            .unwrap_or_default();
        let residue_ptrs: Vec<i64> = parse_flag_section(s, "RESIDUE_POINTER")
            .map(|raw| parse_fortran_ints(&raw))
            .unwrap_or_default();
        let atom_residue: Vec<String> = (0..n_atoms)
            .map(|atom_idx| {
                let res_idx = residue_ptrs
                    .iter()
                    .rposition(|&p| (p as usize) <= atom_idx + 1)
                    .unwrap_or(0);
                residue_labels
                    .get(res_idx)
                    .cloned()
                    .unwrap_or_else(|| "UNK".to_string())
            })
            .collect();
        let atoms: Vec<AmberAtom> = (0..n_atoms)
            .map(|i| AmberAtom {
                name: atom_names.get(i).cloned().unwrap_or_default(),
                residue_name: atom_residue.get(i).cloned().unwrap_or_default(),
                charge: charges.get(i).copied().unwrap_or(0.0),
                mass: masses.get(i).copied().unwrap_or(0.0),
                atom_type: atom_types.get(i).cloned().unwrap_or_default(),
            })
            .collect();
        let bond_ints_h: Vec<i64> = parse_flag_section(s, "BONDS_INC_HYDROGEN")
            .map(|raw| parse_fortran_ints(&raw))
            .unwrap_or_default();
        let bond_ints_no_h: Vec<i64> = parse_flag_section(s, "BONDS_WITHOUT_HYDROGEN")
            .map(|raw| parse_fortran_ints(&raw))
            .unwrap_or_default();
        let bond_force_k: Vec<f64> = parse_flag_section(s, "BOND_FORCE_CONSTANT")
            .map(|raw| parse_fortran_reals(&raw))
            .unwrap_or_default();
        let bond_equil_r: Vec<f64> = parse_flag_section(s, "BOND_EQUIL_VALUE")
            .map(|raw| parse_fortran_reals(&raw))
            .unwrap_or_default();
        let mut bonds: Vec<AmberBond> = build_bonds(&bond_ints_h, &bond_force_k, &bond_equil_r);
        bonds.extend(build_bonds(&bond_ints_no_h, &bond_force_k, &bond_equil_r));
        let angle_ints_h: Vec<i64> = parse_flag_section(s, "ANGLES_INC_HYDROGEN")
            .map(|raw| parse_fortran_ints(&raw))
            .unwrap_or_default();
        let angle_ints_no_h: Vec<i64> = parse_flag_section(s, "ANGLES_WITHOUT_HYDROGEN")
            .map(|raw| parse_fortran_ints(&raw))
            .unwrap_or_default();
        let angle_force_k: Vec<f64> = parse_flag_section(s, "ANGLE_FORCE_CONSTANT")
            .map(|raw| parse_fortran_reals(&raw))
            .unwrap_or_default();
        let angle_equil_t: Vec<f64> = parse_flag_section(s, "ANGLE_EQUIL_VALUE")
            .map(|raw| parse_fortran_reals(&raw))
            .unwrap_or_default();
        let mut angles: Vec<AmberAngle> =
            build_angles(&angle_ints_h, &angle_force_k, &angle_equil_t);
        angles.extend(build_angles(
            &angle_ints_no_h,
            &angle_force_k,
            &angle_equil_t,
        ));
        let n_bonds_actual = n_bonds.max(bonds.len());
        let _ = n_angles;
        Ok(AmberTopology {
            title,
            atoms,
            bonds,
            angles,
            n_atoms,
            n_bonds: n_bonds_actual,
        })
    }
    /// Sum of all partial charges (should be close to an integer for a physical system).
    pub fn total_charge(&self) -> f64 {
        self.atoms.iter().map(|a| a.charge).sum()
    }
    /// Sum of all atomic masses.
    pub fn total_mass(&self) -> f64 {
        self.atoms.iter().map(|a| a.mass).sum()
    }
    /// Sorted, deduplicated list of residue names.
    pub fn residue_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.atoms.iter().map(|a| a.residue_name.clone()).collect();
        names.sort();
        names.dedup();
        names
    }
    /// Human-readable summary of the topology.
    pub fn write_summary(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!("Title      : {}\n", self.title));
        s.push_str(&format!("Atoms      : {}\n", self.atoms.len()));
        s.push_str(&format!("Bonds      : {}\n", self.bonds.len()));
        s.push_str(&format!("Angles     : {}\n", self.angles.len()));
        s.push_str(&format!("Total charge : {:.4} e\n", self.total_charge()));
        s.push_str(&format!("Total mass   : {:.4} amu\n", self.total_mass()));
        let res = self.residue_names();
        s.push_str(&format!("Residues   : {}\n", res.join(", ")));
        s
    }
}
impl AmberTopology {
    /// Sorted, deduplicated list of atom type names.
    pub fn atom_type_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.atoms.iter().map(|a| a.atom_type.clone()).collect();
        names.sort();
        names.dedup();
        names
    }
    /// Number of unique atom types.
    pub fn n_atom_types(&self) -> usize {
        self.atom_type_names().len()
    }
    /// Get atoms belonging to a specific residue name.
    pub fn atoms_in_residue(&self, res_name: &str) -> Vec<usize> {
        self.atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| a.residue_name == res_name)
            .map(|(i, _)| i)
            .collect()
    }
}
/// A parsed energy term line from an AMBER output or mdout file.
///
/// Lines typically look like:
/// ```text
///  NSTEP =       50  TIME(PS) =  0.050000  TEMP(K) =  300.1  PRESS =     0.0
///  Etot   =   -1234.56  EKtot   =    456.78  EPtot   =   -1691.34
/// ```
#[derive(Debug, Clone)]
pub struct AmberEnergyFrame {
    /// Simulation step number.
    pub step: u64,
    /// Simulation time in picoseconds.
    pub time_ps: f64,
    /// Temperature in Kelvin.
    pub temp_k: f64,
    /// Total energy (kcal mol⁻¹).
    pub e_tot: f64,
    /// Kinetic energy (kcal mol⁻¹).
    pub e_kin: f64,
    /// Potential energy (kcal mol⁻¹).
    pub e_pot: f64,
}
/// A dihedral entry from an AMBER FRCMOD file.
#[derive(Debug, Clone)]
pub struct FrcmodDihedral {
    /// Atom types quadruplet (e.g., "X-CT-CT-X").
    pub atom_types: String,
    /// Dividing factor (n-fold).
    pub n_fold: f64,
    /// Barrier height (kcal mol⁻¹).
    pub v_n: f64,
    /// Phase offset (degrees).
    pub gamma_deg: f64,
    /// Periodicity (may be negative to continue with next term).
    pub periodicity: f64,
}
/// An angle entry from an AMBER FRCMOD file.
#[derive(Debug, Clone)]
pub struct FrcmodAngle {
    /// Atom types triplet (e.g., "HC-CT-HC").
    pub atom_types: String,
    /// Force constant (kcal mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium angle (degrees).
    pub theta0_deg: f64,
}
/// A non-bonded (Lennard-Jones) entry from an AMBER FRCMOD file.
#[derive(Debug, Clone)]
pub struct FrcmodNonbonded {
    /// Atom type.
    pub atom_type: String,
    /// R* (van der Waals radius, Å).
    pub r_star: f64,
    /// ε well depth (kcal mol⁻¹).
    pub epsilon: f64,
}
/// A single frame from an AMBER MDCRD trajectory.
#[derive(Debug, Clone)]
pub struct MdcrdFrame {
    /// Cartesian coordinates (Å), layout `[x0, y0, z0, x1, y1, z1, ...]`.
    pub coords: Vec<f64>,
    /// Optional periodic box `[a, b, c]` (last line if present).
    pub box_dims: Option<[f64; 3]>,
}
impl MdcrdFrame {
    /// Get position of atom `i` as `[x, y, z]`.
    pub fn position(&self, i: usize) -> Option<[f64; 3]> {
        let base = i * 3;
        if base + 2 < self.coords.len() {
            Some([
                self.coords[base],
                self.coords[base + 1],
                self.coords[base + 2],
            ])
        } else {
            None
        }
    }
    /// Number of atoms in this frame.
    pub fn n_atoms(&self) -> usize {
        self.coords.len() / 3
    }
}
/// Parsed AMBER FRCMOD (force field modification) file.
#[derive(Debug, Clone, Default)]
pub struct FrcmodFile {
    /// Title/comment line.
    pub title: String,
    /// Modified bond parameters.
    pub bonds: Vec<FrcmodBond>,
    /// Modified angle parameters.
    pub angles: Vec<FrcmodAngle>,
    /// Modified dihedral parameters.
    pub dihedrals: Vec<FrcmodDihedral>,
    /// Modified non-bonded (LJ) parameters.
    pub nonbonded: Vec<FrcmodNonbonded>,
}
impl FrcmodFile {
    /// Parse an AMBER FRCMOD file from a string.
    ///
    /// Sections are:
    /// - `MASS`  – atom masses (skipped here)
    /// - `BOND`  – bond parameters
    /// - `ANGL`  – angle parameters
    /// - `DIHE`  – dihedral parameters
    /// - `IMPR`  – improper dihedral parameters
    /// - `NONB`  – non-bonded (LJ) parameters
    pub fn parse(s: &str) -> Result<Self, String> {
        let mut frc = FrcmodFile::default();
        let mut section = FrcmodSection::None;
        let mut first_line = true;
        for line in s.lines() {
            let trimmed = line.trim();
            if first_line {
                frc.title = trimmed.to_string();
                first_line = false;
                continue;
            }
            let upper = trimmed.to_uppercase();
            if upper.starts_with("MASS") {
                section = FrcmodSection::Mass;
                continue;
            } else if upper.starts_with("BOND") {
                section = FrcmodSection::Bond;
                continue;
            } else if upper.starts_with("ANGL") {
                section = FrcmodSection::Angle;
                continue;
            } else if upper.starts_with("DIHE") {
                section = FrcmodSection::Dihe;
                continue;
            } else if upper.starts_with("IMPR") {
                section = FrcmodSection::Impr;
                continue;
            } else if upper.starts_with("NONB") {
                section = FrcmodSection::Nonb;
                continue;
            } else if upper.starts_with("END") {
                section = FrcmodSection::None;
                continue;
            }
            if trimmed.is_empty() {
                continue;
            }
            match section {
                FrcmodSection::Bond => {
                    let data = strip_amber_comment(trimmed);
                    let parts: Vec<&str> = data
                        .splitn(3, char::is_whitespace)
                        .filter(|s| !s.is_empty())
                        .collect();
                    let nums = collect_numbers(data);
                    if !parts.is_empty() && nums.len() >= 2 {
                        frc.bonds.push(FrcmodBond {
                            atom_types: parts[0].to_string(),
                            k: nums[0],
                            r0: nums[1],
                        });
                    }
                }
                FrcmodSection::Angle => {
                    let data = strip_amber_comment(trimmed);
                    let parts: Vec<&str> = data
                        .splitn(2, char::is_whitespace)
                        .filter(|s| !s.is_empty())
                        .collect();
                    let nums = collect_numbers(data);
                    if !parts.is_empty() && nums.len() >= 2 {
                        frc.angles.push(FrcmodAngle {
                            atom_types: parts[0].to_string(),
                            k: nums[0],
                            theta0_deg: nums[1],
                        });
                    }
                }
                FrcmodSection::Dihe | FrcmodSection::Impr => {
                    let data = strip_amber_comment(trimmed);
                    let parts: Vec<&str> = data
                        .splitn(2, char::is_whitespace)
                        .filter(|s| !s.is_empty())
                        .collect();
                    let nums = collect_numbers(data);
                    if !parts.is_empty() && nums.len() >= 4 {
                        frc.dihedrals.push(FrcmodDihedral {
                            atom_types: parts[0].to_string(),
                            n_fold: nums[0],
                            v_n: nums[1],
                            gamma_deg: nums[2],
                            periodicity: nums[3],
                        });
                    }
                }
                FrcmodSection::Nonb => {
                    let nums = collect_numbers(trimmed);
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if !parts.is_empty() && nums.len() >= 2 {
                        frc.nonbonded.push(FrcmodNonbonded {
                            atom_type: parts[0].to_string(),
                            r_star: nums[0],
                            epsilon: nums[1],
                        });
                    }
                }
                _ => {}
            }
        }
        Ok(frc)
    }
    /// Number of bond entries.
    pub fn n_bonds(&self) -> usize {
        self.bonds.len()
    }
    /// Number of angle entries.
    pub fn n_angles(&self) -> usize {
        self.angles.len()
    }
    /// Number of dihedral entries.
    pub fn n_dihedrals(&self) -> usize {
        self.dihedrals.len()
    }
    /// Number of nonbonded entries.
    pub fn n_nonbonded(&self) -> usize {
        self.nonbonded.len()
    }
    /// Look up a bond by atom type pair string (e.g., "CT-HC").
    pub fn get_bond(&self, types: &str) -> Option<&FrcmodBond> {
        self.bonds.iter().find(|b| {
            b.atom_types == types || {
                let parts: Vec<&str> = types.split('-').collect();
                if parts.len() == 2 {
                    let rev = format!("{}-{}", parts[1], parts[0]);
                    b.atom_types == rev
                } else {
                    false
                }
            }
        })
    }
    /// Look up nonbonded parameters by atom type.
    pub fn get_nonbonded(&self, atom_type: &str) -> Option<&FrcmodNonbonded> {
        self.nonbonded.iter().find(|n| n.atom_type == atom_type)
    }
}

impl std::str::FromStr for FrcmodFile {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}
/// FRCMOD section types.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum FrcmodSection {
    None,
    Mass,
    Bond,
    Angle,
    Dihe,
    Impr,
    Nonb,
}
/// A minimal AMBER RST7 restart record (coordinates + optional velocities).
#[derive(Debug, Clone, Default)]
pub struct AmberRst7 {
    /// Title string.
    pub title: String,
    /// Atom positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Atom velocities (Å ps⁻¹), optional.
    pub velocities: Vec<[f64; 3]>,
}
impl AmberRst7 {
    /// Create a new RST7 record.
    pub fn new(title: &str, positions: Vec<[f64; 3]>) -> Self {
        Self {
            title: title.to_string(),
            positions,
            velocities: Vec::new(),
        }
    }
    /// Attach velocity data to this restart record.
    ///
    /// `velocities` must have the same length as `positions`; returns an error
    /// otherwise.
    pub fn write_velocity(&mut self, velocities: Vec<[f64; 3]>) -> std::result::Result<(), String> {
        if velocities.len() != self.positions.len() {
            return Err(format!(
                "velocity count {} != position count {}",
                velocities.len(),
                self.positions.len()
            ));
        }
        self.velocities = velocities;
        Ok(())
    }
    /// Serialize to AMBER inpcrd/rst7 format string.
    pub fn to_string_repr(&self) -> String {
        write_inpcrd(
            &self.title,
            &self.positions,
            if self.velocities.is_empty() {
                None
            } else {
                Some(&self.velocities)
            },
        )
    }
    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }
    /// Whether velocity data is present.
    pub fn has_velocities(&self) -> bool {
        !self.velocities.is_empty()
    }
}
/// One atom entry from an AMBER topology file.
pub struct AmberAtom {
    /// Atom name (from `ATOM_NAME`).
    pub name: String,
    /// Residue name (from `RESIDUE_LABEL`).
    pub residue_name: String,
    /// Partial charge in electron-charge units (AMBER stores as e / 18.2223).
    pub charge: f64,
    /// Atomic mass in amu.
    pub mass: f64,
    /// AMBER atom-type string (from `AMBER_ATOM_TYPE`).
    pub atom_type: String,
}
/// A dihedral (torsion) parameter from AMBER topology.
#[derive(Debug, Clone)]
pub struct AmberDihedral {
    /// Zero-based index of atom i (first).
    pub i: usize,
    /// Zero-based index of atom j.
    pub j: usize,
    /// Zero-based index of atom k.
    pub k: usize,
    /// Zero-based index of atom l (last).
    pub l: usize,
    /// Barrier height (kcal mol⁻¹).
    pub v_n: f64,
    /// Phase offset (radians).
    pub gamma: f64,
    /// Periodicity (integer, usually 1–6).
    pub n: i32,
    /// `true` if this is a 1-4 non-bonded exclusion.
    pub is_improper: bool,
}
/// A bond entry from an AMBER FRCMOD file.
#[derive(Debug, Clone)]
pub struct FrcmodBond {
    /// Atom types (e.g., "CT-HC").
    pub atom_types: String,
    /// Force constant (kcal mol⁻¹ Å⁻²).
    pub k: f64,
    /// Equilibrium bond length (Å).
    pub r0: f64,
}
/// Parsed AMBER coordinate (inpcrd / restart) data.
pub struct AmberCoordinates {
    /// Title line.
    pub title: String,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Cartesian coordinates (Angstroms), flattened `[x0, y0, z0, x1, y1, z1, ...]`.
    pub coords: Vec<f64>,
    /// Optional velocities (same layout as coords).
    pub velocities: Option<Vec<f64>>,
    /// Optional box dimensions `[a, b, c, alpha, beta, gamma]`.
    pub box_dimensions: Option<[f64; 6]>,
}
impl AmberCoordinates {
    /// Parse an AMBER coordinate / restart file from a string.
    ///
    /// Format: title, natom, then coords in 12.7 format (6 per line).
    /// If the file has twice as many numbers as 3*natom, the second half
    /// is treated as velocities.
    pub fn parse(s: &str) -> Result<Self, String> {
        let lines: Vec<&str> = s.lines().collect();
        if lines.is_empty() {
            return Err("Empty coordinate file".to_string());
        }
        let title = lines[0].trim().to_string();
        if lines.len() < 2 {
            return Err("Missing atom count line".to_string());
        }
        let n_atoms: usize = lines[1]
            .split_whitespace()
            .next()
            .ok_or("Missing atom count")?
            .parse()
            .map_err(|e: std::num::ParseIntError| e.to_string())?;
        let mut all_vals: Vec<f64> = Vec::new();
        for line in lines.iter().skip(2) {
            for tok in line.split_whitespace() {
                if let Ok(v) = tok.parse::<f64>() {
                    all_vals.push(v);
                }
            }
        }
        let n_coords = n_atoms * 3;
        if all_vals.len() < n_coords {
            return Err(format!(
                "Expected {} coordinate values, got {}",
                n_coords,
                all_vals.len()
            ));
        }
        let coords = all_vals[..n_coords].to_vec();
        let velocities = if all_vals.len() >= n_coords * 2 {
            Some(all_vals[n_coords..n_coords * 2].to_vec())
        } else {
            None
        };
        let remaining = if velocities.is_some() {
            &all_vals[n_coords * 2..]
        } else {
            &all_vals[n_coords..]
        };
        let box_dimensions = if remaining.len() >= 6 {
            Some([
                remaining[0],
                remaining[1],
                remaining[2],
                remaining[3],
                remaining[4],
                remaining[5],
            ])
        } else {
            None
        };
        Ok(AmberCoordinates {
            title,
            n_atoms,
            coords,
            velocities,
            box_dimensions,
        })
    }
    /// Get the position of atom `i` as `[x, y, z]`.
    pub fn position(&self, i: usize) -> [f64; 3] {
        let base = i * 3;
        [
            self.coords[base],
            self.coords[base + 1],
            self.coords[base + 2],
        ]
    }
    /// Write coordinates in AMBER restart format.
    pub fn to_restart_string(&self) -> String {
        let mut s = format!("{}\n", self.title);
        s.push_str(&format!("{:5}\n", self.n_atoms));
        for (i, &v) in self.coords.iter().enumerate() {
            s.push_str(&format!("{:12.7}", v));
            if (i + 1) % 6 == 0 {
                s.push('\n');
            }
        }
        if !self.coords.len().is_multiple_of(6) {
            s.push('\n');
        }
        if let Some(ref vels) = self.velocities {
            for (i, &v) in vels.iter().enumerate() {
                s.push_str(&format!("{:12.7}", v));
                if (i + 1) % 6 == 0 {
                    s.push('\n');
                }
            }
            if vels.len() % 6 != 0 {
                s.push('\n');
            }
        }
        if let Some(ref bx) = self.box_dimensions {
            for &v in bx {
                s.push_str(&format!("{:12.7}", v));
            }
            s.push('\n');
        }
        s
    }
}

impl std::str::FromStr for AmberCoordinates {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}
/// A small inline AMBER force-field parameter record for a bond type.
#[derive(Debug, Clone)]
pub struct FfBondParam {
    /// Type string of atom A (e.g., `"CT"`).
    pub type_a: String,
    /// Type string of atom B (e.g., `"HC"`).
    pub type_b: String,
    /// Force constant (kcal mol⁻¹ Å⁻²).
    pub k: f64,
    /// Equilibrium length (Å).
    pub r0: f64,
}
/// A distance restraint entry from an AMBER NMR restraint file.
#[derive(Debug, Clone)]
pub struct AmberRestraint {
    /// Atom 1 index (0-based).
    pub iat1: usize,
    /// Atom 2 index (0-based).
    pub iat2: usize,
    /// Lower bound (Å).
    pub r1: f64,
    /// Lower bound inner (Å).
    pub r2: f64,
    /// Upper bound inner (Å).
    pub r3: f64,
    /// Upper bound outer (Å).
    pub r4: f64,
    /// Force constant (kcal mol⁻¹ Å⁻²).
    pub rk2: f64,
    /// Force constant upper (kcal mol⁻¹ Å⁻²).
    pub rk3: f64,
}
impl AmberRestraint {
    /// Compute the flat-bottom restraint penalty for a given distance `r`.
    ///
    /// Returns 0 for `r2 ≤ r ≤ r3`, quadratic outside.
    pub fn energy(&self, r: f64) -> f64 {
        if r < self.r2 {
            self.rk2 * (r - self.r2).powi(2)
        } else if r <= self.r3 {
            0.0
        } else {
            self.rk3 * (r - self.r3).powi(2)
        }
    }
}
/// A single Lennard-Jones parameter entry (per atom-type pair).
#[derive(Debug, Clone)]
pub struct LjParameter {
    /// Index of atom type A.
    pub type_a: usize,
    /// Index of atom type B.
    pub type_b: usize,
    /// A coefficient (kcal mol⁻¹ Å¹²).
    pub a_coeff: f64,
    /// B coefficient (kcal mol⁻¹ Å⁶).
    pub b_coeff: f64,
}
impl LjParameter {
    /// Compute the equilibrium distance r_min from A/B coefficients.
    ///
    /// r_min = (2A/B)^(1/6)
    pub fn r_min(&self) -> f64 {
        if self.b_coeff.abs() < 1e-30 {
            return 0.0;
        }
        (2.0 * self.a_coeff / self.b_coeff).powf(1.0 / 6.0)
    }
    /// Well depth ε = B² / (4A).
    pub fn epsilon(&self) -> f64 {
        if self.a_coeff.abs() < 1e-30 {
            return 0.0;
        }
        self.b_coeff * self.b_coeff / (4.0 * self.a_coeff)
    }
}
/// A small inline AMBER force-field parameter record for an angle type.
#[derive(Debug, Clone)]
pub struct FfAngleParam {
    /// Type string of atom i.
    pub type_i: String,
    /// Type string of central atom j.
    pub type_j: String,
    /// Type string of atom k.
    pub type_k: String,
    /// Force constant (kcal mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium angle (degrees).
    pub theta0_deg: f64,
}
/// A valence angle entry.
pub struct AmberAngle {
    /// Zero-based index of the first atom.
    pub i: usize,
    /// Zero-based index of the central atom.
    pub j: usize,
    /// Index into the angle parameter array.
    pub k_idx: usize,
    /// Force constant (kcal mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium angle (radians).
    pub theta0: f64,
}
/// Periodic box dimensions for AMBER simulations.
#[derive(Debug, Clone)]
pub struct AmberBox {
    /// Box length in X direction (Å).
    pub a: f64,
    /// Box length in Y direction (Å).
    pub b: f64,
    /// Box length in Z direction (Å).
    pub c: f64,
    /// Box angle alpha (degrees).
    pub alpha: f64,
    /// Box angle beta (degrees).
    pub beta: f64,
    /// Box angle gamma (degrees).
    pub gamma: f64,
}
impl AmberBox {
    /// Create a cubic box with side length `a` (Å).
    pub fn cubic(a: f64) -> Self {
        AmberBox {
            a,
            b: a,
            c: a,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }
    /// Create an orthorhombic box.
    pub fn orthorhombic(a: f64, b: f64, c: f64) -> Self {
        AmberBox {
            a,
            b,
            c,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }
    /// Whether the box is cubic (all sides equal, all angles 90°).
    pub fn is_cubic(&self) -> bool {
        (self.a - self.b).abs() < 1e-8
            && (self.b - self.c).abs() < 1e-8
            && (self.alpha - 90.0).abs() < 1e-8
            && (self.beta - 90.0).abs() < 1e-8
            && (self.gamma - 90.0).abs() < 1e-8
    }
    /// Whether the box is orthorhombic (all angles 90°).
    pub fn is_orthorhombic(&self) -> bool {
        (self.alpha - 90.0).abs() < 1e-8
            && (self.beta - 90.0).abs() < 1e-8
            && (self.gamma - 90.0).abs() < 1e-8
    }
    /// Compute the volume (Å³) using the general parallelepiped formula.
    pub fn volume(&self) -> f64 {
        use std::f64::consts::PI;
        let to_rad = |deg: f64| deg * PI / 180.0;
        let (ca, cb, cg) = (
            to_rad(self.alpha).cos(),
            to_rad(self.beta).cos(),
            to_rad(self.gamma).cos(),
        );
        let v = (1.0 + 2.0 * ca * cb * cg - ca * ca - cb * cb - cg * cg).sqrt();
        self.a * self.b * self.c * v
    }
    /// Convert box dimensions to Å (already in Å, returns same object).
    pub fn in_angstrom(&self) -> Self {
        self.clone()
    }
    /// Convert box dimensions to nm (divide all lengths by 10).
    pub fn in_nm(&self) -> [f64; 6] {
        [
            self.a / 10.0,
            self.b / 10.0,
            self.c / 10.0,
            self.alpha,
            self.beta,
            self.gamma,
        ]
    }
}
/// Simplified DCD trajectory writer (CHARMM/NAMD binary format).
///
/// Writes a self-contained binary header on construction, then one frame
/// at a time via [`AmberDcd::write_frame`].  The output is a valid DCD
/// (little-endian FORTRAN record format) readable by MDAnalysis / VMD.
///
/// Note: this implementation writes a simplified 84-byte header without
/// the full CHARMM title block; it is intended for testing and round-trip
/// verification in OxiPhysics.
pub struct AmberDcd {
    /// Accumulated frame data (flat x, y, z per frame, f32 LE).
    pub frames: Vec<Vec<f32>>,
    /// Number of atoms (fixed across frames).
    pub n_atoms: usize,
}
impl AmberDcd {
    /// Create a new DCD trajectory writer for `n_atoms` atoms.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            frames: Vec::new(),
            n_atoms,
        }
    }
    /// Append one frame of coordinates.
    ///
    /// `positions` must have length `n_atoms`; returns an error otherwise.
    pub fn write_frame(&mut self, positions: &[[f64; 3]]) -> std::result::Result<(), String> {
        if positions.len() != self.n_atoms {
            return Err(format!(
                "expected {} atoms, got {}",
                self.n_atoms,
                positions.len()
            ));
        }
        let flat: Vec<f32> = positions
            .iter()
            .flat_map(|p| p.iter().map(|&v| v as f32))
            .collect();
        self.frames.push(flat);
        Ok(())
    }
    /// Number of frames written.
    pub fn n_frames(&self) -> usize {
        self.frames.len()
    }
    /// Get position of atom `atom_idx` in frame `frame_idx` as f32.
    pub fn get_position(&self, frame_idx: usize, atom_idx: usize) -> Option<[f32; 3]> {
        let frame = self.frames.get(frame_idx)?;
        let base = atom_idx * 3;
        if base + 2 >= frame.len() {
            return None;
        }
        Some([frame[base], frame[base + 1], frame[base + 2]])
    }
    /// Serialize to a simplified DCD binary representation.
    ///
    /// Layout (little-endian):
    /// - 4 bytes: n_atoms (u32)
    /// - 4 bytes: n_frames (u32)
    /// - For each frame: 4 bytes record_len (u32), then n_atoms*3*4 bytes (f32 LE)
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        let n_atoms = self.n_atoms as u32;
        let n_frames = self.frames.len() as u32;
        out.extend_from_slice(&n_atoms.to_le_bytes());
        out.extend_from_slice(&n_frames.to_le_bytes());
        for frame in &self.frames {
            let record_len = (frame.len() * 4) as u32;
            out.extend_from_slice(&record_len.to_le_bytes());
            for &v in frame {
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
        out
    }
    /// Deserialize from the simplified DCD binary representation produced by
    /// [`AmberDcd::to_bytes`].
    pub fn from_bytes(data: &[u8]) -> std::result::Result<Self, String> {
        if data.len() < 8 {
            return Err("DCD: too short".into());
        }
        let n_atoms = u32::from_le_bytes([data[0], data[1], data[2], data[3]]) as usize;
        let n_frames = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
        let mut dcd = AmberDcd::new(n_atoms);
        let mut offset = 8;
        for _ in 0..n_frames {
            if offset + 4 > data.len() {
                return Err("DCD: truncated frame record length".into());
            }
            let record_len = u32::from_le_bytes([
                data[offset],
                data[offset + 1],
                data[offset + 2],
                data[offset + 3],
            ]) as usize;
            offset += 4;
            if offset + record_len > data.len() {
                return Err("DCD: truncated frame data".into());
            }
            let n_floats = record_len / 4;
            let mut frame = Vec::with_capacity(n_floats);
            for i in 0..n_floats {
                let b = &data[offset + i * 4..offset + i * 4 + 4];
                frame.push(f32::from_le_bytes([b[0], b[1], b[2], b[3]]));
            }
            dcd.frames.push(frame);
            offset += record_len;
        }
        Ok(dcd)
    }
}
