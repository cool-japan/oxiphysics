// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Crystallography file format I/O.
//!
//! Provides parsers and writers for common crystallography formats:
//! - CIF (Crystallographic Information File): cell, atom sites and
//!   space-group symmetry expansion ([`CifParser`], [`write_cif`])
//! - VASP POSCAR/CONTCAR, including selective dynamics, Cartesian or direct
//!   coordinates and negative (volume) scale factors ([`VaspReader`],
//!   [`PoscarStructure`])
//! - powder XRD pattern simulation and analysis ([`simulate_xrd`],
//!   [`simulate_xrd_with_fwhm`], [`XrdPattern`])

use crate::{Error, Result};
use std::collections::HashMap;

mod cif_atoms;
mod cif_syntax;
mod lattice_math;
mod symmetry;
mod xrd;

pub use xrd::{
    DEFAULT_XRD_FWHM_DEG, lorentz_polarization_factor, simulate_xrd, simulate_xrd_with_fwhm,
};

#[cfg(test)]
mod cif_tests;
#[cfg(test)]
mod poscar_tests;

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// A single atom site in a crystal structure.
#[derive(Debug, Clone)]
pub struct CrystalAtom {
    /// Chemical element symbol (e.g. "Fe", "O").
    pub element: String,
    /// Fractional coordinates within the unit cell \[a, b, c\].
    pub fractional_pos: [f64; 3],
    /// Site occupancy (0.0 – 1.0).
    pub occupancy: f64,
    /// Isotropic displacement parameter (B-factor) in Å².
    pub b_factor: f64,
}

/// A periodic crystal structure described by its lattice and atom list.
#[derive(Debug, Clone)]
pub struct CrystalStructure {
    /// Lattice vectors as row matrix: `lattice[i]` is the i-th basis vector in Å.
    pub lattice: [[f64; 3]; 3],
    /// List of atoms in the asymmetric unit (or full unit cell).
    pub atoms: Vec<CrystalAtom>,
    /// International Tables space group number (1–230).
    pub space_group: u16,
}

/// Parsed X-ray diffraction pattern.
#[derive(Debug, Clone)]
pub struct XrdPattern {
    /// 2θ angles in degrees.
    pub two_theta: Vec<f64>,
    /// Intensity values (arbitrary units).
    pub intensity: Vec<f64>,
}

// ---------------------------------------------------------------------------
// CIF parser
// ---------------------------------------------------------------------------

/// Parser for CIF (Crystallographic Information File) format.
#[derive(Debug, Default)]
pub struct CifParser;

impl CifParser {
    /// Create a new `CifParser`.
    pub fn new() -> Self {
        Self
    }

    /// Parse a CIF file `content` and return a [`CrystalStructure`].
    ///
    /// Reads the first data block.  Handles `_cell_length_a/b/c`,
    /// `_cell_angle_alpha/beta/gamma`, the `_atom_site_*` table (fractional or
    /// Cartesian coordinates, `B_iso` or `U_iso` displacement) wherever it sits
    /// among the other tables, and items that follow a table.  When the file
    /// lists symmetry operations (`_space_group_symop_operation_xyz` or
    /// `_symmetry_equiv_pos_as_xyz`) the atom sites are expanded to the whole
    /// unit cell: every position is wrapped into `[0, 1)` and symmetry-equivalent
    /// duplicates (same element within 1e-4) are merged.
    ///
    /// # Errors
    ///
    /// Fails on malformed CIF syntax, an atom coordinate or occupancy that is
    /// not a number, a truncated atom table, or an unreadable symmetry
    /// operation.
    pub fn parse(&self, content: &str) -> Result<CrystalStructure> {
        let block = cif_syntax::parse_block(content)?;
        let kv = &block.items;

        // Extract cell parameters
        let a = Self::parse_float_kv(kv, "_cell_length_a").unwrap_or(1.0);
        let b = Self::parse_float_kv(kv, "_cell_length_b").unwrap_or(1.0);
        let c = Self::parse_float_kv(kv, "_cell_length_c").unwrap_or(1.0);
        let alpha = Self::parse_float_kv(kv, "_cell_angle_alpha")
            .unwrap_or(90.0)
            .to_radians();
        let beta = Self::parse_float_kv(kv, "_cell_angle_beta")
            .unwrap_or(90.0)
            .to_radians();
        let gamma = Self::parse_float_kv(kv, "_cell_angle_gamma")
            .unwrap_or(90.0)
            .to_radians();

        let lattice = Self::build_lattice(a, b, c, alpha, beta, gamma);

        // Space group
        let space_group = kv
            .get("_space_group_it_number")
            .or_else(|| kv.get("_symmetry_int_tables_number"))
            .and_then(|s| s.parse::<u16>().ok())
            .unwrap_or(1);

        // The atom-site table is the loop that carries the site coordinates;
        // other tables (e.g. anisotropic displacement parameters) repeat the
        // label column and must not be mistaken for it.
        let atoms = match block
            .loop_with_any(&cif_atoms::FRACT_X_TAGS)
            .or_else(|| block.loop_with_any(&cif_atoms::CARTN_X_TAGS))
        {
            Some(table) => cif_atoms::extract_atoms(table, &lattice)?,
            None => Vec::new(),
        };
        let operations = cif_atoms::collect_symmetry_ops(&block)?;
        let atoms = if operations.is_empty() {
            atoms
        } else {
            symmetry::expand_atoms(&atoms, &operations)
        };

        Ok(CrystalStructure {
            lattice,
            atoms,
            space_group,
        })
    }

    /// Build a lattice matrix from cell parameters (conventional setting).
    fn build_lattice(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> [[f64; 3]; 3] {
        let cos_a = alpha.cos();
        let cos_b = beta.cos();
        let cos_g = gamma.cos();
        let sin_g = gamma.sin();

        // Standard triclinic convention
        let cx = cos_b;
        let cy = (cos_a - cos_b * cos_g) / sin_g;
        let cz = (1.0 - cx * cx - cy * cy).max(0.0).sqrt();

        [
            [a, 0.0, 0.0],
            [b * cos_g, b * sin_g, 0.0],
            [c * cx, c * cy, c * cz],
        ]
    }

    /// Parse a float from a key-value map entry, ignoring standard uncertainties like "1.234(5)".
    fn parse_float_kv(kv: &HashMap<String, String>, key: &str) -> Option<f64> {
        let s = kv.get(key)?;
        Self::strip_su(s).parse().ok()
    }

    /// Strip standard uncertainty "(n)" from a numeric string.
    fn strip_su(s: &str) -> &str {
        if let Some(pos) = s.find('(') {
            &s[..pos]
        } else {
            s
        }
    }
}

// ---------------------------------------------------------------------------
// CrystalStructure methods
// ---------------------------------------------------------------------------

impl CrystalStructure {
    /// Create a simple cubic structure with given lattice parameter.
    pub fn new_cubic(a: f64) -> Self {
        Self {
            lattice: [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
            atoms: Vec::new(),
            space_group: 1,
        }
    }

    /// Convert all atom fractional coordinates to Cartesian coordinates (Å).
    pub fn to_cartesian(&self) -> Vec<[f64; 3]> {
        self.atoms
            .iter()
            .map(|atom| {
                let [u, v, w] = atom.fractional_pos;
                let x = u * self.lattice[0][0] + v * self.lattice[1][0] + w * self.lattice[2][0];
                let y = u * self.lattice[0][1] + v * self.lattice[1][1] + w * self.lattice[2][1];
                let z = u * self.lattice[0][2] + v * self.lattice[1][2] + w * self.lattice[2][2];
                [x, y, z]
            })
            .collect()
    }

    /// Volume of the unit cell in Å³ (scalar triple product of basis vectors).
    pub fn volume(&self) -> f64 {
        let a = self.lattice[0];
        let b = self.lattice[1];
        let c = self.lattice[2];
        // a · (b × c)
        let bxc = [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ];
        (a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2]).abs()
    }

    /// Crystal density in g/cm³.
    ///
    /// `molar_mass` is the molar mass of one formula unit in g/mol.
    /// The number of formula units Z is taken as the atom count.
    pub fn density(&self, molar_mass: f64) -> f64 {
        const NA: f64 = 6.022_140_76e23; // Avogadro
        const ANG3_TO_CM3: f64 = 1e-24;
        let z = self.atoms.len() as f64;
        let vol_cm3 = self.volume() * ANG3_TO_CM3;
        (z * molar_mass) / (NA * vol_cm3)
    }

    /// Compute the reciprocal lattice vectors (2π convention).
    ///
    /// Returns `[[b1x,b1y,b1z\],[b2x,b2y,b2z],[b3x,b3y,b3z]]` in Å⁻¹.
    pub fn reciprocal_lattice(&self) -> [[f64; 3]; 3] {
        let a1 = self.lattice[0];
        let a2 = self.lattice[1];
        let a3 = self.lattice[2];
        let vol = self.volume();

        let cross = |u: [f64; 3], v: [f64; 3]| -> [f64; 3] {
            [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]
        };
        let scale = |s: f64, u: [f64; 3]| -> [f64; 3] { [s * u[0], s * u[1], s * u[2]] };

        let b1 = scale(2.0 * std::f64::consts::PI / vol, cross(a2, a3));
        let b2 = scale(2.0 * std::f64::consts::PI / vol, cross(a3, a1));
        let b3 = scale(2.0 * std::f64::consts::PI / vol, cross(a1, a2));

        [b1, b2, b3]
    }

    /// Interplanar d-spacing for Miller indices (h,k,l) in Å.
    ///
    /// Uses the general triclinic formula via the reciprocal lattice.
    pub fn d_spacing(&self, h: i32, k: i32, l: i32) -> f64 {
        let rec = self.reciprocal_lattice();
        // G = h*b1 + k*b2 + l*b3
        let hf = h as f64;
        let kf = k as f64;
        let lf = l as f64;
        let gx = hf * rec[0][0] + kf * rec[1][0] + lf * rec[2][0];
        let gy = hf * rec[0][1] + kf * rec[1][1] + lf * rec[2][1];
        let gz = hf * rec[0][2] + kf * rec[1][2] + lf * rec[2][2];
        let g_mag = (gx * gx + gy * gy + gz * gz).sqrt();
        if g_mag < 1e-12 {
            f64::INFINITY
        } else {
            2.0 * std::f64::consts::PI / g_mag
        }
    }

    /// Compute the kinematic structure factor F(h,k,l) = (Re, Im).
    ///
    /// Uses atomic scattering factors approximated by atomic number (simplified).
    pub fn structure_factor(&self, h: i32, k: i32, l: i32) -> (f64, f64) {
        let mut re = 0.0_f64;
        let mut im = 0.0_f64;

        for atom in &self.atoms {
            let [u, v, w] = atom.fractional_pos;
            let phase = 2.0 * std::f64::consts::PI * (h as f64 * u + k as f64 * v + l as f64 * w);
            // Simplified scattering factor: f ≈ atomic_number(element) * exp(-B*sin²θ/λ²)
            let f = Self::approx_scattering_factor(&atom.element) * atom.occupancy;
            re += f * phase.cos();
            im += f * phase.sin();
        }
        (re, im)
    }

    /// Approximate atomic scattering factor by element symbol (f(0)).
    fn approx_scattering_factor(element: &str) -> f64 {
        // Very simplified: use atomic number as proxy for f(0)
        match element.trim_end_matches(|c: char| c.is_ascii_digit()) {
            "H" => 1.0,
            "He" => 2.0,
            "Li" => 3.0,
            "Be" => 4.0,
            "B" => 5.0,
            "C" => 6.0,
            "N" => 7.0,
            "O" => 8.0,
            "F" => 9.0,
            "Na" => 11.0,
            "Mg" => 12.0,
            "Al" => 13.0,
            "Si" => 14.0,
            "P" => 15.0,
            "S" => 16.0,
            "Cl" => 17.0,
            "Ar" => 18.0,
            "K" => 19.0,
            "Ca" => 20.0,
            "Ti" => 22.0,
            "V" => 23.0,
            "Cr" => 24.0,
            "Mn" => 25.0,
            "Fe" => 26.0,
            "Co" => 27.0,
            "Ni" => 28.0,
            "Cu" => 29.0,
            "Zn" => 30.0,
            "Ga" => 31.0,
            "Ge" => 32.0,
            "As" => 33.0,
            "Se" => 34.0,
            "Br" => 35.0,
            "Sr" => 38.0,
            "Mo" => 42.0,
            "Ag" => 47.0,
            "Sn" => 50.0,
            "Ba" => 56.0,
            "La" => 57.0,
            "W" => 74.0,
            "Pt" => 78.0,
            "Au" => 79.0,
            "Hg" => 80.0,
            "Pb" => 82.0,
            "U" => 92.0,
            _ => 10.0, // fallback
        }
    }
}

// ---------------------------------------------------------------------------
// CIF writer
// ---------------------------------------------------------------------------

/// Write a [`CrystalStructure`] to a CIF-formatted string.
pub fn write_cif(crystal: &CrystalStructure) -> String {
    let mut out = String::new();

    out.push_str("# CIF written by oxiphysics-io\n");
    out.push_str("data_oxiphysics\n\n");

    // Cell parameters
    let a = crystal.lattice[0];
    let b = crystal.lattice[1];
    let c = crystal.lattice[2];

    let a_len = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    let b_len = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
    let c_len = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();

    let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let alpha = if b_len > 0.0 && c_len > 0.0 {
        (dot(b, c) / (b_len * c_len))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    } else {
        90.0
    };
    let beta = if a_len > 0.0 && c_len > 0.0 {
        (dot(a, c) / (a_len * c_len))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    } else {
        90.0
    };
    let gamma = if a_len > 0.0 && b_len > 0.0 {
        (dot(a, b) / (a_len * b_len))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    } else {
        90.0
    };

    out.push_str(&format!("_cell_length_a    {a_len:.6}\n"));
    out.push_str(&format!("_cell_length_b    {b_len:.6}\n"));
    out.push_str(&format!("_cell_length_c    {c_len:.6}\n"));
    out.push_str(&format!("_cell_angle_alpha {alpha:.4}\n"));
    out.push_str(&format!("_cell_angle_beta  {beta:.4}\n"));
    out.push_str(&format!("_cell_angle_gamma {gamma:.4}\n"));
    out.push_str(&format!(
        "_symmetry_int_tables_number {}\n\n",
        crystal.space_group
    ));

    // Atom loop
    out.push_str("loop_\n");
    out.push_str("  _atom_site_type_symbol\n");
    out.push_str("  _atom_site_fract_x\n");
    out.push_str("  _atom_site_fract_y\n");
    out.push_str("  _atom_site_fract_z\n");
    out.push_str("  _atom_site_occupancy\n");
    out.push_str("  _atom_site_b_iso\n");

    for atom in &crystal.atoms {
        out.push_str(&format!(
            "  {:4} {:10.6} {:10.6} {:10.6} {:6.4} {:8.4}\n",
            atom.element,
            atom.fractional_pos[0],
            atom.fractional_pos[1],
            atom.fractional_pos[2],
            atom.occupancy,
            atom.b_factor,
        ));
    }

    out
}

// ---------------------------------------------------------------------------
// VASP POSCAR reader
// ---------------------------------------------------------------------------

/// A richer POSCAR structure that captures VASP-specific extensions.
///
/// Unlike [`CrystalStructure`] (which maps to CIF), this struct preserves
/// VASP-specific fields: selective dynamics flags, per-atom magnetic moments,
/// and ionic forces parsed from a companion OUTCAR file.
#[derive(Debug, Clone)]
pub struct PoscarStructure {
    /// Comment line (first line of the POSCAR file).
    pub comment: String,
    /// Effective (positive) scale factor that [`VaspReader`] applied to
    /// `lattice` and to any Cartesian input coordinates.  A negative factor in
    /// the file is a target cell volume and is reported here as the equivalent
    /// positive scale.
    pub scale: f64,
    /// Lattice vectors as row matrix (in Å after scale).
    pub lattice: [[f64; 3]; 3],
    /// Element symbols for each species (VASP5+ format; empty for VASP4).
    pub species: Vec<String>,
    /// Number of atoms per species.
    pub counts: Vec<usize>,
    /// Whether selective dynamics mode is active.
    pub selective_dynamics: bool,
    /// Coordinate mode of `positions`: `true` = Direct (fractional), `false` =
    /// Cartesian in Å.  [`VaspReader`] always returns direct coordinates; a
    /// Cartesian input block is converted with the full (possibly skewed)
    /// lattice.
    pub is_direct: bool,
    /// Atom positions in species order: fractional when `is_direct` is set,
    /// otherwise Cartesian in Å with the lattice scale already applied.
    pub positions: Vec<[f64; 3]>,
    /// Per-atom selective-dynamics flags `[x_free, y_free, z_free]`.
    /// `None` if selective dynamics were not specified.
    pub sd_flags: Option<Vec<[bool; 3]>>,
    /// Per-atom magnetic moments parsed from `# MAGMOM = ...` annotation.
    /// `None` if not present.
    pub magmom: Option<Vec<f64>>,
    /// Ionic forces in eV/Å parsed from a companion OUTCAR file.
    /// `None` if not available.
    pub forces: Option<Vec<[f64; 3]>>,
}

impl PoscarStructure {
    /// Atom positions as fractional coordinates.
    ///
    /// Returns the stored positions when [`is_direct`](Self::is_direct) is
    /// set, otherwise converts the Cartesian ones (Å) with the inverse of the
    /// full lattice.
    ///
    /// # Errors
    ///
    /// Fails when the structure is in Cartesian mode and the lattice is
    /// singular, so no conversion exists.
    pub fn fractional_positions(&self) -> Result<Vec<[f64; 3]>> {
        if self.is_direct {
            return Ok(self.positions.clone());
        }
        let inverse = lattice_math::inverse(&self.lattice).ok_or_else(|| {
            Error::Parse("POSCAR: lattice is singular, cannot convert Cartesian coordinates".into())
        })?;
        Ok(self
            .positions
            .iter()
            .map(|&r| lattice_math::cartesian_to_fractional(r, &inverse))
            .collect())
    }

    /// Convert to a [`CrystalStructure`] (drops VASP-specific fields).
    ///
    /// Positions of a Cartesian-mode structure are converted to fractional
    /// coordinates (see [`fractional_positions`](Self::fractional_positions));
    /// when the lattice is singular no such conversion exists and the stored
    /// positions are carried over unchanged.
    pub fn to_crystal_structure(&self) -> CrystalStructure {
        let fractional = self
            .fractional_positions()
            .unwrap_or_else(|_| self.positions.clone());
        let mut element_iter = self.species.iter();
        let mut count_iter = self.counts.iter();
        let mut atoms = Vec::with_capacity(fractional.len());
        let mut current_element = element_iter
            .next()
            .cloned()
            .unwrap_or_else(|| "X".to_string());
        let mut remaining = count_iter.next().copied().unwrap_or(0);

        for pos in &fractional {
            if remaining == 0 {
                current_element = element_iter
                    .next()
                    .cloned()
                    .unwrap_or_else(|| "X".to_string());
                remaining = count_iter.next().copied().unwrap_or(0);
            }
            atoms.push(CrystalAtom {
                element: current_element.clone(),
                fractional_pos: *pos,
                occupancy: 1.0,
                b_factor: 0.0,
            });
            remaining = remaining.saturating_sub(1);
        }
        CrystalStructure {
            lattice: self.lattice,
            atoms,
            space_group: 1,
        }
    }

    /// Write this structure to a POSCAR-formatted string.
    pub fn to_poscar_string(&self) -> String {
        let mut out = String::new();
        out.push_str(&self.comment);
        out.push('\n');
        out.push_str("1.0\n");
        for row in &self.lattice {
            out.push_str(&format!(
                "  {:.10}  {:.10}  {:.10}\n",
                row[0], row[1], row[2]
            ));
        }
        if !self.species.is_empty() {
            out.push_str(&self.species.join("  "));
            out.push('\n');
        }
        let count_strs: Vec<String> = self.counts.iter().map(|c| c.to_string()).collect();
        out.push_str(&count_strs.join("  "));
        out.push('\n');
        if self.selective_dynamics {
            out.push_str("Selective dynamics\n");
        }
        if self.is_direct {
            out.push_str("Direct\n");
        } else {
            out.push_str("Cartesian\n");
        }
        for (i, pos) in self.positions.iter().enumerate() {
            if self.selective_dynamics {
                let flags = self
                    .sd_flags
                    .as_ref()
                    .and_then(|f| f.get(i))
                    .copied()
                    .unwrap_or([true, true, true]);
                let fx = if flags[0] { "T" } else { "F" };
                let fy = if flags[1] { "T" } else { "F" };
                let fz = if flags[2] { "T" } else { "F" };
                out.push_str(&format!(
                    "  {:.10}  {:.10}  {:.10}  {fx}  {fy}  {fz}\n",
                    pos[0], pos[1], pos[2]
                ));
            } else {
                out.push_str(&format!(
                    "  {:.10}  {:.10}  {:.10}\n",
                    pos[0], pos[1], pos[2]
                ));
            }
        }
        if let Some(ref mm) = self.magmom {
            let mm_strs: Vec<String> = mm.iter().map(|v| format!("{v:.4}")).collect();
            out.push_str(&format!("# MAGMOM = {}\n", mm_strs.join(" ")));
        }
        out
    }
}

/// Parse OUTCAR ionic forces from a VASP OUTCAR file.
///
/// Looks for the last `TOTAL-FORCE (eV/Angst)` block and reads the
/// following `n_atoms` lines as `[fx, fy, fz]` force vectors.
pub fn parse_outcar_forces(content: &str, n_atoms: usize) -> Option<Vec<[f64; 3]>> {
    // Find the last occurrence of the TOTAL-FORCE block.
    let marker = "TOTAL-FORCE (eV/Angst)";
    let last_pos = content.rfind(marker)?;
    let after = &content[last_pos..];
    // Skip the header line and the dashed separator line.
    let mut lines = after.lines().skip(1);
    // Skip separator line (contains dashes)
    let sep = lines.next()?;
    if !sep.contains('-') {
        return None;
    }
    let mut forces = Vec::with_capacity(n_atoms);
    for line in lines.take(n_atoms) {
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() < 6 {
            return None;
        }
        let fx: f64 = toks.get(3).and_then(|s| s.parse().ok())?;
        let fy: f64 = toks.get(4).and_then(|s| s.parse().ok())?;
        let fz: f64 = toks.get(5).and_then(|s| s.parse().ok())?;
        forces.push([fx, fy, fz]);
    }
    if forces.len() == n_atoms {
        Some(forces)
    } else {
        None
    }
}

/// Reader for VASP POSCAR/CONTCAR format.
#[derive(Debug, Default)]
pub struct VaspReader;

impl VaspReader {
    /// Create a new `VaspReader`.
    pub fn new() -> Self {
        Self
    }

    /// Parse a POSCAR/CONTCAR file and return a [`CrystalStructure`].
    ///
    /// Handles the header, scale (a negative factor is the target cell volume
    /// in Å^3), lattice vectors, optional selective dynamics flags, and
    /// Cartesian/Direct coordinate blocks.  The returned atoms always carry
    /// fractional coordinates.
    pub fn parse_poscar(&self, content: &str) -> Result<CrystalStructure> {
        let ps = self.parse_poscar_full(content)?;
        Ok(ps.to_crystal_structure())
    }

    /// Parse a POSCAR/CONTCAR file and return a full [`PoscarStructure`],
    /// preserving selective dynamics flags, magnetic moments, and other
    /// VASP-specific information.
    ///
    /// The universal scale factor multiplies the lattice vectors; a negative
    /// factor is instead the cell volume in Å^3 the lattice is scaled to.
    /// Cartesian coordinates are multiplied by the same scale and converted to
    /// fractional coordinates with the full lattice, so the structure is
    /// returned in direct mode (`is_direct`) whichever mode the file used and
    /// [`PoscarStructure::scale`] holds the effective positive scale.
    ///
    /// # Errors
    ///
    /// Fails on a malformed header or lattice, a zero or non-finite scale, a
    /// volume scale (or Cartesian coordinates) combined with a singular
    /// lattice.
    pub fn parse_poscar_full(&self, content: &str) -> Result<PoscarStructure> {
        let mut lines = content.lines().peekable();

        // Line 1: comment
        let comment = lines.next().unwrap_or("").trim().to_string();
        // Line 2: universal scaling factor
        let scale_factor: f64 = lines
            .next()
            .unwrap_or("1.0")
            .trim()
            .parse()
            .map_err(|_| Error::Parse("POSCAR: invalid scale factor".into()))?;
        if !scale_factor.is_finite() || scale_factor == 0.0 {
            return Err(Error::Parse(
                "POSCAR: scale factor must be finite and non-zero".into(),
            ));
        }

        // Lines 3-5: lattice vectors (still unscaled)
        let mut lattice = [[0.0_f64; 3]; 3];
        for row in lattice.iter_mut() {
            let line = lines
                .next()
                .ok_or_else(|| Error::Parse("POSCAR: missing lattice vector".into()))?;
            let vals: Vec<f64> = line
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if vals.len() < 3 {
                return Err(Error::Parse("POSCAR: bad lattice vector".into()));
            }
            row[0] = vals[0];
            row[1] = vals[1];
            row[2] = vals[2];
        }

        // A positive factor multiplies the lattice.  A negative one is the
        // target cell volume in Å^3 (VASP convention): the unscaled lattice
        // is scaled uniformly until its volume matches.
        let scale = if scale_factor > 0.0 {
            scale_factor
        } else {
            let volume = lattice_math::determinant(&lattice).abs();
            if volume < lattice_math::MIN_LATTICE_DET {
                return Err(Error::Parse(
                    "POSCAR: a negative (volume) scale factor needs a non-singular lattice".into(),
                ));
            }
            (-scale_factor / volume).cbrt()
        };
        for row in lattice.iter_mut() {
            for value in row.iter_mut() {
                *value *= scale;
            }
        }

        // Line 6: element symbols (VASP5+) or direct species count
        let species_line = lines.next().unwrap_or("").trim().to_string();
        let (species, count_line) = if species_line
            .chars()
            .next()
            .map(|c| c.is_alphabetic())
            .unwrap_or(false)
        {
            let elems: Vec<String> = species_line
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            let cl = lines.next().unwrap_or("").trim().to_string();
            (elems, cl)
        } else {
            // VASP4: no element names — use placeholder
            (vec!["X".to_string()], species_line.clone())
        };

        let counts: Vec<usize> = count_line
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        let n_atoms: usize = counts.iter().sum();

        // Next line: optional "Selective dynamics", then coordinate mode
        let next_line = lines.next().unwrap_or("").trim().to_string();
        let (selective_dynamics, coord_type_line) = if next_line.to_lowercase().starts_with('s') {
            let cl = lines.next().unwrap_or("").trim().to_string();
            (true, cl)
        } else {
            (false, next_line)
        };

        let cartesian_input = coord_type_line.to_lowercase().starts_with('c')
            || coord_type_line.to_lowercase().starts_with('k');
        // Cartesian input is converted to fractional coordinates with the
        // inverse of the full scaled lattice (not just its diagonal).
        let lattice_inverse = if cartesian_input {
            Some(lattice_math::inverse(&lattice).ok_or_else(|| {
                Error::Parse(
                    "POSCAR: lattice is singular, cannot convert Cartesian coordinates".into(),
                )
            })?)
        } else {
            None
        };

        // Read atom positions (and optional selective dynamics flags)
        let mut positions: Vec<[f64; 3]> = Vec::with_capacity(n_atoms);
        let mut sd_flags: Vec<[bool; 3]> = Vec::with_capacity(n_atoms);
        // Also capture any trailing lines for MAGMOM
        let mut trailing: Vec<String> = Vec::new();

        for _ in 0..n_atoms {
            let line = lines.next().unwrap_or("0 0 0");
            let toks: Vec<&str> = line.split_whitespace().collect();
            let x: f64 = toks.first().and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let y: f64 = toks.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let z: f64 = toks.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let mut pos = [x, y, z];
            if let Some(inverse) = &lattice_inverse {
                // Cartesian coordinates carry the same universal scale factor
                // as the lattice vectors.
                let cartesian = [pos[0] * scale, pos[1] * scale, pos[2] * scale];
                pos = lattice_math::cartesian_to_fractional(cartesian, inverse);
            }
            positions.push(pos);

            if selective_dynamics && toks.len() >= 6 {
                let fx = toks[3] == "T";
                let fy = toks[4] == "T";
                let fz = toks[5] == "T";
                sd_flags.push([fx, fy, fz]);
            } else {
                sd_flags.push([true, true, true]);
            }
        }

        // Consume remaining lines to find MAGMOM annotations.
        for line in lines {
            trailing.push(line.trim().to_string());
        }

        // Parse MAGMOM from trailing comments: `# MAGMOM = 1.0 -1.0 ...`
        let magmom = trailing
            .iter()
            .find(|l| {
                let lower = l.to_lowercase();
                lower.contains("magmom") && lower.contains('=')
            })
            .and_then(|l| {
                // Strip leading `#` and everything up to and including `=`
                let after_eq = l.split('=').nth(1)?;
                let vals: Vec<f64> = after_eq
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.is_empty() { None } else { Some(vals) }
            });

        let sd_opt = if selective_dynamics {
            Some(sd_flags)
        } else {
            None
        };

        Ok(PoscarStructure {
            comment,
            scale,
            lattice,
            species,
            counts,
            selective_dynamics,
            is_direct: true,
            positions,
            sd_flags: sd_opt,
            magmom,
            forces: None,
        })
    }
}

/// Parse a string containing one or more POSCAR images (e.g., VASP NEB chain).
///
/// Images are separated by double-newlines or a blank line followed by
/// a new POSCAR header. Returns a `Vec<PoscarStructure>` with one entry per image.
pub fn read_poscar_images(content: &str) -> Result<Vec<PoscarStructure>> {
    // Split on blank lines (double newline) to find image boundaries.
    // Each image starts with a non-empty comment line followed by scale etc.
    let reader = VaspReader::new();

    // Split the content on sequences of two or more consecutive newlines.
    let blocks: Vec<&str> = content
        .split("\n\n")
        .map(|b| b.trim())
        .filter(|b| !b.is_empty())
        .collect();

    if blocks.is_empty() {
        return Err(Error::Parse("no POSCAR images found".into()));
    }

    let mut images = Vec::with_capacity(blocks.len());
    for block in &blocks {
        let ps = reader.parse_poscar_full(block)?;
        images.push(ps);
    }
    Ok(images)
}

// ---------------------------------------------------------------------------
// XRD pattern analysis
// ---------------------------------------------------------------------------

impl XrdPattern {
    /// Create an XRD pattern from parallel 2θ and intensity arrays.
    pub fn new(two_theta: Vec<f64>, intensity: Vec<f64>) -> Self {
        Self {
            two_theta,
            intensity,
        }
    }

    /// Find peak positions (local maxima) using a simple derivative method.
    ///
    /// Returns indices into `two_theta`/`intensity` where peaks occur.
    pub fn find_peaks(&self) -> Vec<usize> {
        let n = self.intensity.len();
        if n < 3 {
            return vec![];
        }
        let mut peaks = Vec::new();
        for i in 1..n - 1 {
            if self.intensity[i] > self.intensity[i - 1]
                && self.intensity[i] > self.intensity[i + 1]
            {
                peaks.push(i);
            }
        }
        peaks
    }

    /// Convert 2θ peak positions to d-spacings using Bragg's law: d = λ/(2 sin θ).
    ///
    /// Uses Cu Kα radiation (λ = 1.5406 Å) by default.
    pub fn d_spacings(&self) -> Vec<f64> {
        const LAMBDA: f64 = 1.5406; // Å, Cu Kα
        self.find_peaks()
            .iter()
            .map(|&i| {
                let theta = self.two_theta[i].to_radians() / 2.0;
                let sin_t = theta.sin();
                if sin_t > 1e-12 {
                    LAMBDA / (2.0 * sin_t)
                } else {
                    f64::INFINITY
                }
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- CifParser ---

    #[test]
    fn test_cif_parser_minimal_cubic() {
        let cif = r#"
data_test
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
_symmetry_int_tables_number 225
"#;
        let parser = CifParser::new();
        let crystal = parser.parse(cif).unwrap();
        assert!((crystal.volume() - 125.0).abs() < 0.01);
        assert_eq!(crystal.space_group, 225);
    }

    #[test]
    fn test_cif_parser_space_group_fallback() {
        let cif = "data_x\n_cell_length_a 3.0\n_cell_length_b 3.0\n_cell_length_c 3.0\n";
        let parser = CifParser::new();
        let crystal = parser.parse(cif).unwrap();
        assert_eq!(crystal.space_group, 1);
    }

    #[test]
    fn test_cif_parser_with_standard_uncertainty() {
        let cif =
            "data_x\n_cell_length_a 4.123(5)\n_cell_length_b 4.123(5)\n_cell_length_c 4.123(5)\n";
        let parser = CifParser::new();
        let crystal = parser.parse(cif).unwrap();
        assert!((crystal.lattice[0][0] - 4.123).abs() < 0.001);
    }

    #[test]
    fn test_volume_cubic() {
        let crystal = CrystalStructure::new_cubic(3.0);
        assert!((crystal.volume() - 27.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_orthorhombic() {
        let crystal = CrystalStructure {
            lattice: [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            atoms: vec![],
            space_group: 1,
        };
        assert!((crystal.volume() - 24.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_cartesian_fractional_origin() {
        let mut crystal = CrystalStructure::new_cubic(4.0);
        crystal.atoms.push(CrystalAtom {
            element: "Fe".into(),
            fractional_pos: [0.0, 0.0, 0.0],
            occupancy: 1.0,
            b_factor: 0.0,
        });
        let cart = crystal.to_cartesian();
        assert!((cart[0][0]).abs() < 1e-12);
        assert!((cart[0][1]).abs() < 1e-12);
        assert!((cart[0][2]).abs() < 1e-12);
    }

    #[test]
    fn test_to_cartesian_fractional_half() {
        let mut crystal = CrystalStructure::new_cubic(4.0);
        crystal.atoms.push(CrystalAtom {
            element: "O".into(),
            fractional_pos: [0.5, 0.5, 0.5],
            occupancy: 1.0,
            b_factor: 0.0,
        });
        let cart = crystal.to_cartesian();
        assert!((cart[0][0] - 2.0).abs() < 1e-10);
        assert!((cart[0][1] - 2.0).abs() < 1e-10);
        assert!((cart[0][2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_d_spacing_cubic_100() {
        // For cubic with a=3: d(1,0,0) = a = 3 Å
        let crystal = CrystalStructure::new_cubic(3.0);
        let d = crystal.d_spacing(1, 0, 0);
        assert!((d - 3.0).abs() < 0.001, "d_100 = {d}");
    }

    #[test]
    fn test_d_spacing_cubic_110() {
        // d(1,1,0) = a/√2
        let a = 4.0;
        let crystal = CrystalStructure::new_cubic(a);
        let d = crystal.d_spacing(1, 1, 0);
        let expected = a / 2.0_f64.sqrt();
        assert!(
            (d - expected).abs() < 0.001,
            "d_110 = {d}, expected {expected}"
        );
    }

    #[test]
    fn test_d_spacing_cubic_111() {
        // d(1,1,1) = a/√3
        let a = 3.0;
        let crystal = CrystalStructure::new_cubic(a);
        let d = crystal.d_spacing(1, 1, 1);
        let expected = a / 3.0_f64.sqrt();
        assert!(
            (d - expected).abs() < 0.001,
            "d_111 = {d}, expected {expected}"
        );
    }

    #[test]
    fn test_reciprocal_lattice_cubic() {
        let a = 2.0 * std::f64::consts::PI;
        let crystal = CrystalStructure::new_cubic(a);
        let rec = crystal.reciprocal_lattice();
        // For cubic: b1 = 2π/a * x̂ → magnitude = 1.0
        assert!((rec[0][0] - 1.0).abs() < 1e-10);
        assert!((rec[1][1] - 1.0).abs() < 1e-10);
        assert!((rec[2][2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_structure_factor_single_atom_at_origin() {
        let mut crystal = CrystalStructure::new_cubic(4.0);
        crystal.atoms.push(CrystalAtom {
            element: "Fe".into(),
            fractional_pos: [0.0, 0.0, 0.0],
            occupancy: 1.0,
            b_factor: 0.0,
        });
        // Any hkl: phase = 0, so F = f(Fe) + 0i
        let (re, im) = crystal.structure_factor(1, 0, 0);
        assert!(re > 0.0, "real part should be positive");
        assert!(im.abs() < 1e-10, "imaginary part should be 0");
    }

    #[test]
    fn test_density_positive() {
        let mut crystal = CrystalStructure::new_cubic(3.615); // copper
        for _ in 0..4 {
            crystal.atoms.push(CrystalAtom {
                element: "Cu".into(),
                fractional_pos: [0.0, 0.0, 0.0],
                occupancy: 1.0,
                b_factor: 0.0,
            });
        }
        let rho = crystal.density(63.546); // Cu molar mass
        assert!(rho > 0.0, "density must be positive");
        // copper density ~8.9 g/cm³
        assert!(
            rho > 5.0 && rho < 15.0,
            "density {rho} out of expected range"
        );
    }

    #[test]
    fn test_write_cif_roundtrip_cell() {
        let crystal = CrystalStructure::new_cubic(5.0);
        let cif = write_cif(&crystal);
        assert!(cif.contains("_cell_length_a"));
        assert!(cif.contains("5.000000"));
    }

    #[test]
    fn test_write_cif_contains_atom_loop() {
        let mut crystal = CrystalStructure::new_cubic(4.0);
        crystal.atoms.push(CrystalAtom {
            element: "Na".into(),
            fractional_pos: [0.0, 0.0, 0.0],
            occupancy: 1.0,
            b_factor: 0.5,
        });
        let cif = write_cif(&crystal);
        assert!(cif.contains("loop_"));
        assert!(cif.contains("Na"));
    }

    #[test]
    fn test_vasp_reader_simple_poscar() {
        let poscar = "Fe BCC\n1.0\n2.87 0.0 0.0\n0.0 2.87 0.0\n0.0 0.0 2.87\nFe\n2\nDirect\n0.0 0.0 0.0\n0.5 0.5 0.5\n";
        let reader = VaspReader::new();
        let crystal = reader.parse_poscar(poscar).unwrap();
        assert_eq!(crystal.atoms.len(), 2);
        assert_eq!(crystal.atoms[0].element, "Fe");
    }

    #[test]
    fn test_vasp_reader_lattice_scale() {
        let poscar =
            "Test\n2.0\n1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\nH\n1\nDirect\n0.0 0.0 0.0\n";
        let reader = VaspReader::new();
        let crystal = reader.parse_poscar(poscar).unwrap();
        // Scale=2 applied to lattice
        assert!((crystal.lattice[0][0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_xrd_pattern_find_peaks() {
        // Manually create a pattern with clear peaks
        let two_theta: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let mut intensity = vec![0.0; 20];
        intensity[5] = 100.0;
        intensity[10] = 200.0;
        let pattern = XrdPattern::new(two_theta, intensity);
        let peaks = pattern.find_peaks();
        assert!(peaks.contains(&5), "peak at index 5");
        assert!(peaks.contains(&10), "peak at index 10");
    }

    #[test]
    fn test_xrd_d_spacings_bragg() {
        // 2θ = 44.5° for Fe(110) with Cu Kα → d ~ 2.03 Å
        let two_theta = vec![40.0, 44.5, 50.0];
        let intensity = vec![5.0, 100.0, 5.0];
        let pattern = XrdPattern::new(two_theta, intensity);
        let ds = pattern.d_spacings();
        assert!(!ds.is_empty());
        // d for 2θ=44.5° ≈ 2.03 Å
        assert!(ds[0] > 1.5 && ds[0] < 3.0, "d = {}", ds[0]);
    }

    #[test]
    fn test_simulate_xrd_returns_correct_size() {
        let crystal = CrystalStructure::new_cubic(3.0);
        let pattern = simulate_xrd(&crystal, 1.5406, (20.0, 80.0), 100);
        assert_eq!(pattern.two_theta.len(), 100);
        assert_eq!(pattern.intensity.len(), 100);
    }

    #[test]
    fn test_simulate_xrd_non_negative_intensity() {
        let mut crystal = CrystalStructure::new_cubic(3.615);
        crystal.atoms.push(CrystalAtom {
            element: "Cu".into(),
            fractional_pos: [0.0, 0.0, 0.0],
            occupancy: 1.0,
            b_factor: 0.0,
        });
        let pattern = simulate_xrd(&crystal, 1.5406, (20.0, 80.0), 200);
        for &i in &pattern.intensity {
            assert!(i >= 0.0, "negative intensity: {i}");
        }
    }

    #[test]
    fn test_simulate_xrd_has_peaks() {
        let mut crystal = CrystalStructure::new_cubic(3.615);
        crystal.atoms.push(CrystalAtom {
            element: "Cu".into(),
            fractional_pos: [0.0, 0.0, 0.0],
            occupancy: 1.0,
            b_factor: 0.0,
        });
        let pattern = simulate_xrd(&crystal, 1.5406, (20.0, 100.0), 500);
        let peaks = pattern.find_peaks();
        assert!(!peaks.is_empty(), "simulated XRD should have peaks");
    }

    #[test]
    fn test_approx_scattering_factor_known() {
        assert!((CrystalStructure::approx_scattering_factor("Fe") - 26.0).abs() < 0.5);
        assert!((CrystalStructure::approx_scattering_factor("O") - 8.0).abs() < 0.5);
        assert!((CrystalStructure::approx_scattering_factor("Cu") - 29.0).abs() < 0.5);
    }

    #[test]
    fn test_build_lattice_cubic_is_diagonal() {
        let a = 4.0;
        let pi2 = std::f64::consts::PI / 2.0;
        let lat = CifParser::build_lattice(a, a, a, pi2, pi2, pi2);
        assert!((lat[0][0] - a).abs() < 1e-10);
        assert!(lat[0][1].abs() < 1e-10);
        assert!(lat[0][2].abs() < 1e-10);
        assert!((lat[1][1] - a).abs() < 1e-8, "b[1] = {}", lat[1][1]);
    }

    #[test]
    fn test_strip_su() {
        assert_eq!(CifParser::strip_su("1.234(5)"), "1.234");
        assert_eq!(CifParser::strip_su("5.000"), "5.000");
    }

    #[test]
    fn test_d_spacing_000_is_infinity() {
        let crystal = CrystalStructure::new_cubic(3.0);
        let d = crystal.d_spacing(0, 0, 0);
        assert!(d.is_infinite());
    }

    #[test]
    fn test_crystal_structure_no_atoms_cartesian() {
        let crystal = CrystalStructure::new_cubic(5.0);
        assert!(crystal.to_cartesian().is_empty());
    }

    #[test]
    fn test_cif_parser_empty_returns_default() {
        let parser = CifParser::new();
        let crystal = parser.parse("").unwrap();
        assert!((crystal.volume() - 1.0).abs() < 1e-10); // 1×1×1
    }

    // -----------------------------------------------------------------------
    // J2: Extended POSCAR reader tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_poscar_selective_dynamics_flags_preserved() {
        // POSCAR with selective dynamics: atom 0 is fully free, atom 1 has z fixed.
        let poscar = concat!(
            "BCC Fe\n",
            "1.0\n",
            "2.87 0.0 0.0\n",
            "0.0 2.87 0.0\n",
            "0.0 0.0 2.87\n",
            "Fe\n",
            "2\n",
            "Selective dynamics\n",
            "Direct\n",
            "0.0 0.0 0.0 T T T\n",
            "0.5 0.5 0.5 T T F\n",
        );
        let reader = VaspReader::new();
        let ps = reader.parse_poscar_full(poscar).expect("parse failed");

        assert!(ps.selective_dynamics, "selective_dynamics flag must be set");
        let flags = ps.sd_flags.as_ref().expect("sd_flags must be Some");
        assert_eq!(flags.len(), 2, "must have 2 flag entries");

        // Atom 0: all free
        assert!(flags[0][0], "atom0 x should be free (T)");
        assert!(flags[0][1], "atom0 y should be free (T)");
        assert!(flags[0][2], "atom0 z should be free (T)");

        // Atom 1: z fixed
        assert!(flags[1][0], "atom1 x should be free (T)");
        assert!(flags[1][1], "atom1 y should be free (T)");
        assert!(!flags[1][2], "atom1 z should be fixed (F)");
    }

    #[test]
    fn test_poscar_without_selective_dynamics() {
        let poscar = concat!(
            "Simple\n",
            "1.0\n",
            "3.0 0.0 0.0\n",
            "0.0 3.0 0.0\n",
            "0.0 0.0 3.0\n",
            "H\n",
            "1\n",
            "Direct\n",
            "0.0 0.0 0.0\n",
        );
        let reader = VaspReader::new();
        let ps = reader.parse_poscar_full(poscar).expect("parse failed");
        assert!(!ps.selective_dynamics, "selective_dynamics must be false");
        assert!(
            ps.sd_flags.is_none(),
            "sd_flags must be None when not specified"
        );
    }

    #[test]
    fn test_poscar_multi_image_vec_length() {
        // Two POSCAR images separated by double newline.
        let img1 = concat!(
            "Image1\n",
            "1.0\n",
            "3.0 0.0 0.0\n",
            "0.0 3.0 0.0\n",
            "0.0 0.0 3.0\n",
            "H\n",
            "1\n",
            "Direct\n",
            "0.0 0.0 0.0\n",
        );
        let img2 = concat!(
            "Image2\n",
            "1.0\n",
            "4.0 0.0 0.0\n",
            "0.0 4.0 0.0\n",
            "0.0 0.0 4.0\n",
            "O\n",
            "2\n",
            "Direct\n",
            "0.0 0.0 0.0\n",
            "0.5 0.5 0.5\n",
        );
        let combined = format!("{img1}\n{img2}");
        let images = read_poscar_images(&combined).expect("read_poscar_images failed");
        assert_eq!(images.len(), 2, "must parse 2 images");
        // First image has lattice 3.0, second 4.0
        assert!(
            (images[0].lattice[0][0] - 3.0).abs() < 1e-10,
            "image0 lattice a=3"
        );
        assert!(
            (images[1].lattice[0][0] - 4.0).abs() < 1e-10,
            "image1 lattice a=4"
        );
        // Species and counts
        assert_eq!(images[0].species[0], "H");
        assert_eq!(images[1].species[0], "O");
        assert_eq!(images[1].positions.len(), 2, "image1 should have 2 atoms");
    }

    #[test]
    fn test_poscar_multi_image_different_lattices() {
        // Verify fractional positions are preserved correctly across images.
        let img1 = concat!(
            "Img A\n",
            "1.0\n",
            "5.0 0.0 0.0\n",
            "0.0 5.0 0.0\n",
            "0.0 0.0 5.0\n",
            "Al\n",
            "1\n",
            "Direct\n",
            "0.25 0.25 0.25\n",
        );
        let img2 = concat!(
            "Img B\n",
            "1.0\n",
            "6.0 0.0 0.0\n",
            "0.0 6.0 0.0\n",
            "0.0 0.0 6.0\n",
            "Al\n",
            "1\n",
            "Direct\n",
            "0.75 0.75 0.75\n",
        );
        let combined = format!("{img1}\n{img2}");
        let images = read_poscar_images(&combined).expect("read_poscar_images failed");
        assert_eq!(images.len(), 2);
        assert!((images[0].positions[0][0] - 0.25).abs() < 1e-10);
        assert!((images[1].positions[0][0] - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_poscar_roundtrip_atom_positions() {
        // Build a PoscarStructure, write it, read it back, compare positions.
        let original = PoscarStructure {
            comment: "Test roundtrip".into(),
            scale: 1.0,
            lattice: [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]],
            species: vec!["Fe".into(), "O".into()],
            counts: vec![1, 2],
            selective_dynamics: false,
            is_direct: true,
            positions: vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5]],
            sd_flags: None,
            magmom: None,
            forces: None,
        };

        let poscar_str = original.to_poscar_string();
        let reader = VaspReader::new();
        let recovered = reader
            .parse_poscar_full(&poscar_str)
            .expect("roundtrip parse failed");

        assert_eq!(recovered.positions.len(), 3, "must recover 3 atoms");
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (recovered.positions[i][j] - original.positions[i][j]).abs() < 1e-8,
                    "position[{i}][{j}] mismatch: {} vs {}",
                    recovered.positions[i][j],
                    original.positions[i][j]
                );
            }
        }
    }

    #[test]
    fn test_poscar_magmom_parsing() {
        let poscar = concat!(
            "Fe NiO\n",
            "1.0\n",
            "4.2 0.0 0.0\n",
            "0.0 4.2 0.0\n",
            "0.0 0.0 4.2\n",
            "Fe Ni\n",
            "1 1\n",
            "Direct\n",
            "0.0 0.0 0.0\n",
            "0.5 0.5 0.5\n",
            "# MAGMOM = 4.0 -4.0\n",
        );
        let reader = VaspReader::new();
        let ps = reader.parse_poscar_full(poscar).expect("parse failed");
        let mm = ps.magmom.as_ref().expect("magmom must be Some");
        assert_eq!(mm.len(), 2, "should have 2 magnetic moments");
        assert!((mm[0] - 4.0).abs() < 1e-9, "mm[0] should be 4.0");
        assert!((mm[1] - (-4.0)).abs() < 1e-9, "mm[1] should be -4.0");
    }

    #[test]
    fn test_parse_outcar_forces_basic() {
        // Simulate a minimal OUTCAR excerpt with force data.
        let outcar = concat!(
            " some header\n",
            " TOTAL-FORCE (eV/Angst)\n",
            " ---------------\n",
            "  0.0  0.0  0.0  0.1  0.2  0.3\n",
            "  0.5  0.5  0.5 -0.1 -0.2 -0.3\n",
        );
        let forces = parse_outcar_forces(outcar, 2).expect("parse_outcar_forces failed");
        assert_eq!(forces.len(), 2, "must parse 2 force entries");
        assert!((forces[0][0] - 0.1).abs() < 1e-10, "fx[0] mismatch");
        assert!((forces[0][1] - 0.2).abs() < 1e-10, "fy[0] mismatch");
        assert!((forces[0][2] - 0.3).abs() < 1e-10, "fz[0] mismatch");
        assert!((forces[1][0] - (-0.1)).abs() < 1e-10, "fx[1] mismatch");
        assert!((forces[1][1] - (-0.2)).abs() < 1e-10, "fy[1] mismatch");
        assert!((forces[1][2] - (-0.3)).abs() < 1e-10, "fz[1] mismatch");
    }
}
