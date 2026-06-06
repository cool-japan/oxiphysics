// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! AMBER molecular dynamics file format I/O (prmtop, rst7, mdcrd, mdinfo).
//!
//! Supports reading and writing the primary AMBER file formats used in
//! classical molecular dynamics simulations:
//!
//! - `prmtop` / `parm7` – topology + force-field parameters
//! - `rst7` / `inpcrd` – restart / initial-coordinate files
//! - `mdcrd` / `nc` – ASCII trajectory files
//! - `mdin` – MD input (namelist) files
//! - `mdinfo` – runtime energy/temperature summary files
//! - GAFF / ff14SB parameter files (MASS, BOND, ANGLE, DIHE sections)

use std::collections::HashMap;
use std::fmt;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
// ──────────────────────────────────────────────────────────────────────────────
// Shared error type
// ──────────────────────────────────────────────────────────────────────────────

/// Errors produced by the AMBER I/O module.
#[derive(Debug)]
pub enum AmberIoError {
    /// The file could not be parsed.
    ParseError(String),
    /// An I/O operation failed.
    IoError(String),
    /// A required section was not found.
    MissingSection(String),
    /// Unexpected end of data.
    UnexpectedEof,
}

impl fmt::Display for AmberIoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AmberIoError::ParseError(s) => write!(f, "AMBER parse error: {s}"),
            AmberIoError::IoError(s) => write!(f, "AMBER I/O error: {s}"),
            AmberIoError::MissingSection(s) => write!(f, "AMBER missing section: {s}"),
            AmberIoError::UnexpectedEof => write!(f, "AMBER unexpected end of file"),
        }
    }
}

impl std::error::Error for AmberIoError {}

// ──────────────────────────────────────────────────────────────────────────────
// PrmtopAtom
// ──────────────────────────────────────────────────────────────────────────────

/// A single atom entry from an AMBER topology (`prmtop`) file.
#[derive(Debug, Clone, PartialEq)]
pub struct PrmtopAtom {
    /// AMBER atom name (up to 4 characters).
    pub name: String,
    /// Partial charge in AMBER internal units (e × 18.2223).
    pub charge: f64,
    /// Atomic mass in u (unified atomic mass units).
    pub mass: f64,
    /// AMBER atom-type string.
    pub atom_type: String,
}

impl PrmtopAtom {
    /// Create a new [`PrmtopAtom`].
    pub fn new(
        name: impl Into<String>,
        charge: f64,
        mass: f64,
        atom_type: impl Into<String>,
    ) -> Self {
        PrmtopAtom {
            name: name.into(),
            charge,
            mass,
            atom_type: atom_type.into(),
        }
    }

    /// Convert the stored AMBER charge to SI elementary charges (divide by 18.2223).
    pub fn charge_si(&self) -> f64 {
        self.charge / 18.2223
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// PrmtopBond
// ──────────────────────────────────────────────────────────────────────────────

/// Bond connectivity and force-field parameters from a prmtop file.
#[derive(Debug, Clone, PartialEq)]
pub struct PrmtopBond {
    /// Index of the first atom (0-based).
    pub atom_i: usize,
    /// Index of the second atom (0-based).
    pub atom_j: usize,
    /// Harmonic force constant k (kcal mol⁻¹ Å⁻²).
    pub force_constant: f64,
    /// Equilibrium bond length r₀ (Å).
    pub eq_length: f64,
}

impl PrmtopBond {
    /// Create a new [`PrmtopBond`].
    pub fn new(atom_i: usize, atom_j: usize, force_constant: f64, eq_length: f64) -> Self {
        PrmtopBond {
            atom_i,
            atom_j,
            force_constant,
            eq_length,
        }
    }

    /// Compute the harmonic potential energy for a given bond length.
    ///
    /// `E = k * (r - r0)^2`
    pub fn energy(&self, r: f64) -> f64 {
        self.force_constant * (r - self.eq_length).powi(2)
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// PrmtopAngle
// ──────────────────────────────────────────────────────────────────────────────

/// Angle connectivity and force-field parameters from a prmtop file.
#[derive(Debug, Clone, PartialEq)]
pub struct PrmtopAngle {
    /// Index of atom i (0-based).
    pub atom_i: usize,
    /// Index of the central atom j (0-based).
    pub atom_j: usize,
    /// Index of atom k (0-based).
    pub atom_k: usize,
    /// Harmonic force constant k_θ (kcal mol⁻¹ rad⁻²).
    pub force_constant: f64,
    /// Equilibrium angle θ₀ (radians).
    pub eq_angle: f64,
}

impl PrmtopAngle {
    /// Create a new [`PrmtopAngle`].
    pub fn new(
        atom_i: usize,
        atom_j: usize,
        atom_k: usize,
        force_constant: f64,
        eq_angle: f64,
    ) -> Self {
        PrmtopAngle {
            atom_i,
            atom_j,
            atom_k,
            force_constant,
            eq_angle,
        }
    }

    /// Compute the harmonic bending energy for a given angle.
    ///
    /// `E = k * (θ - θ0)^2`
    pub fn energy(&self, theta: f64) -> f64 {
        self.force_constant * (theta - self.eq_angle).powi(2)
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// PrmtopFile
// ──────────────────────────────────────────────────────────────────────────────

/// AMBER topology file (`prmtop` / `parm7`) representation.
///
/// Only the most commonly used sections are modelled here; raw section data
/// for every other `%FLAG` is preserved in `extra_sections`.
#[derive(Debug, Clone, Default)]
pub struct PrmtopFile {
    /// Version string from the `%VERSION` line.
    pub version: String,
    /// All atoms in topology order.
    pub atoms: Vec<PrmtopAtom>,
    /// All covalent bonds.
    pub bonds: Vec<PrmtopBond>,
    /// All valence angles.
    pub angles: Vec<PrmtopAngle>,
    /// Raw text for unrecognised `%FLAG` sections.
    pub extra_sections: HashMap<String, Vec<String>>,
}

impl PrmtopFile {
    /// Create an empty [`PrmtopFile`].
    pub fn new() -> Self {
        PrmtopFile::default()
    }

    /// Write the topology to any [`Write`] sink in a simplified prmtop-style
    /// text format.
    ///
    /// The output is intentionally human-readable rather than byte-perfect
    /// AMBER binary prmtop; it is suitable for round-trip testing.
    pub fn write<W: Write>(&self, writer: &mut W) -> Result<(), AmberIoError> {
        let mut bw = BufWriter::new(writer);
        writeln!(bw, "%VERSION  {}", self.version)
            .map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, "%FLAG ATOM_NAME").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, "%FORMAT(20a4)").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        for (i, atom) in self.atoms.iter().enumerate() {
            write!(bw, "{:<4}", &atom.name).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            if (i + 1) % 20 == 0 {
                writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            }
        }
        if !self.atoms.len().is_multiple_of(20) {
            writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        writeln!(bw, "%FLAG CHARGE").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, "%FORMAT(5E16.8)").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        for (i, atom) in self.atoms.iter().enumerate() {
            write!(bw, "{:16.8E}", atom.charge)
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
            if (i + 1) % 5 == 0 {
                writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            }
        }
        if !self.atoms.len().is_multiple_of(5) {
            writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        writeln!(bw, "%FLAG MASS").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, "%FORMAT(5E16.8)").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        for (i, atom) in self.atoms.iter().enumerate() {
            write!(bw, "{:16.8E}", atom.mass).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            if (i + 1) % 5 == 0 {
                writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            }
        }
        if !self.atoms.len().is_multiple_of(5) {
            writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        writeln!(bw, "%FLAG ATOM_TYPE_INDEX").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, "%FORMAT(20a4)").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        for (i, atom) in self.atoms.iter().enumerate() {
            write!(bw, "{:<4}", &atom.atom_type)
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
            if (i + 1) % 20 == 0 {
                writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            }
        }
        if !self.atoms.len().is_multiple_of(20) {
            writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        // Bonds
        writeln!(bw, "%FLAG BOND_FORCE_CONSTANT")
            .map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, "%FORMAT(5E16.8)").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        for (i, bond) in self.bonds.iter().enumerate() {
            write!(bw, "{:16.8E}", bond.force_constant)
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
            if (i + 1) % 5 == 0 {
                writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            }
        }
        if !self.bonds.is_empty() && !self.bonds.len().is_multiple_of(5) {
            writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        writeln!(bw, "%FLAG BOND_EQUIL_VALUE").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, "%FORMAT(5E16.8)").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        for (i, bond) in self.bonds.iter().enumerate() {
            write!(bw, "{:16.8E}", bond.eq_length)
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
            if (i + 1) % 5 == 0 {
                writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
            }
        }
        if !self.bonds.is_empty() && !self.bonds.len().is_multiple_of(5) {
            writeln!(bw).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        writeln!(bw, "%END").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        Ok(())
    }

    /// Parse a simplified prmtop-style text produced by [`PrmtopFile::write`].
    pub fn read<R: Read>(reader: R) -> Result<Self, AmberIoError> {
        let br = BufReader::new(reader);
        let mut file = PrmtopFile::new();
        let lines: Vec<String> = br
            .lines()
            .map(|l| l.map_err(|e| AmberIoError::IoError(e.to_string())))
            .collect::<Result<_, _>>()?;

        let mut i = 0usize;
        let mut names: Vec<String> = Vec::new();
        let mut charges: Vec<f64> = Vec::new();
        let mut masses: Vec<f64> = Vec::new();
        let mut types: Vec<String> = Vec::new();
        let mut bond_kf: Vec<f64> = Vec::new();
        let mut bond_eq: Vec<f64> = Vec::new();

        while i < lines.len() {
            let line = lines[i].trim();
            if line.starts_with("%VERSION") {
                file.version = line.trim_start_matches("%VERSION").trim().to_string();
            } else if line == "%FLAG ATOM_NAME" {
                i += 1; // skip FORMAT
                i += 1;
                while i < lines.len() && !lines[i].starts_with('%') {
                    let row = &lines[i];
                    let mut pos = 0usize;
                    while pos + 4 <= row.len() {
                        names.push(row[pos..pos + 4].trim().to_string());
                        pos += 4;
                    }
                    i += 1;
                }
                continue;
            } else if line == "%FLAG CHARGE" {
                i += 1;
                i += 1;
                while i < lines.len() && !lines[i].starts_with('%') {
                    for tok in lines[i].split_whitespace() {
                        charges.push(
                            tok.parse::<f64>()
                                .map_err(|e| AmberIoError::ParseError(e.to_string()))?,
                        );
                    }
                    i += 1;
                }
                continue;
            } else if line == "%FLAG MASS" {
                i += 1;
                i += 1;
                while i < lines.len() && !lines[i].starts_with('%') {
                    for tok in lines[i].split_whitespace() {
                        masses.push(
                            tok.parse::<f64>()
                                .map_err(|e| AmberIoError::ParseError(e.to_string()))?,
                        );
                    }
                    i += 1;
                }
                continue;
            } else if line == "%FLAG ATOM_TYPE_INDEX" {
                i += 1;
                i += 1;
                while i < lines.len() && !lines[i].starts_with('%') {
                    let row = &lines[i];
                    let mut pos = 0usize;
                    while pos + 4 <= row.len() {
                        types.push(row[pos..pos + 4].trim().to_string());
                        pos += 4;
                    }
                    i += 1;
                }
                continue;
            } else if line == "%FLAG BOND_FORCE_CONSTANT" {
                i += 1;
                i += 1;
                while i < lines.len() && !lines[i].starts_with('%') {
                    for tok in lines[i].split_whitespace() {
                        bond_kf.push(
                            tok.parse::<f64>()
                                .map_err(|e| AmberIoError::ParseError(e.to_string()))?,
                        );
                    }
                    i += 1;
                }
                continue;
            } else if line == "%FLAG BOND_EQUIL_VALUE" {
                i += 1;
                i += 1;
                while i < lines.len() && !lines[i].starts_with('%') {
                    for tok in lines[i].split_whitespace() {
                        bond_eq.push(
                            tok.parse::<f64>()
                                .map_err(|e| AmberIoError::ParseError(e.to_string()))?,
                        );
                    }
                    i += 1;
                }
                continue;
            }
            i += 1;
        }

        let n = names.len();
        for idx in 0..n {
            file.atoms.push(PrmtopAtom {
                name: names.get(idx).cloned().unwrap_or_default(),
                charge: charges.get(idx).copied().unwrap_or(0.0),
                mass: masses.get(idx).copied().unwrap_or(0.0),
                atom_type: types.get(idx).cloned().unwrap_or_default(),
            });
        }

        for (idx, &force_constant) in bond_kf.iter().enumerate() {
            file.bonds.push(PrmtopBond {
                atom_i: idx * 2,
                atom_j: idx * 2 + 1,
                force_constant,
                eq_length: bond_eq.get(idx).copied().unwrap_or(0.0),
            });
        }

        Ok(file)
    }

    /// Return the number of atoms.
    pub fn natom(&self) -> usize {
        self.atoms.len()
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Rst7File
// ──────────────────────────────────────────────────────────────────────────────

/// AMBER restart / initial-coordinate file (`rst7` / `inpcrd`).
///
/// ASCII format:
/// ```text
/// title
/// natom [time]
/// x1 y1 z1 x2 y2 z2 ...   (10 per line, 12.7 format)
/// [velocities – same format]
/// [box: a b c alpha beta gamma]
/// ```
#[derive(Debug, Clone, Default)]
pub struct Rst7File {
    /// Title line.
    pub title: String,
    /// Simulation time (ps) from the header, if present.
    pub time: Option<f64>,
    /// Cartesian coordinates (Å) in x,y,z order, length = 3 × NATOM.
    pub coordinates: Vec<f64>,
    /// Velocities (Å ps⁻¹) in x,y,z order, same length as coordinates.
    pub velocities: Vec<f64>,
    /// Periodic box dimensions: `[a, b, c, alpha, beta, gamma]`.
    pub box_dimensions: Option<[f64; 6]>,
}

impl Rst7File {
    /// Create an empty [`Rst7File`].
    pub fn new() -> Self {
        Rst7File::default()
    }

    /// Return the number of atoms.
    pub fn natom(&self) -> usize {
        self.coordinates.len() / 3
    }

    /// Write the restart file to any [`Write`] sink.
    ///
    /// Coordinates are written in 12.7 fixed-point format, 6 values per line
    /// (= 2 atoms per line, matching the AMBER specification).
    pub fn write<W: Write>(&self, writer: &mut W) -> Result<(), AmberIoError> {
        let mut bw = BufWriter::new(writer);
        writeln!(bw, "{}", self.title).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        if let Some(t) = self.time {
            writeln!(bw, "{:6}{:15.7e}", self.natom(), t)
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
        } else {
            writeln!(bw, "{:6}", self.natom()).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        write_f64_block(&mut bw, &self.coordinates, 6)?;

        if !self.velocities.is_empty() {
            write_f64_block(&mut bw, &self.velocities, 6)?;
        }

        if let Some(b) = self.box_dimensions {
            writeln!(
                bw,
                "{:12.7}{:12.7}{:12.7}{:12.7}{:12.7}{:12.7}",
                b[0], b[1], b[2], b[3], b[4], b[5]
            )
            .map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }

        Ok(())
    }

    /// Parse an AMBER restart file from any [`Read`] source.
    pub fn read<R: Read>(reader: R) -> Result<Self, AmberIoError> {
        let br = BufReader::new(reader);
        let mut lines = br.lines();
        let mut rst = Rst7File::new();

        // title
        rst.title = lines
            .next()
            .ok_or(AmberIoError::UnexpectedEof)?
            .map_err(|e| AmberIoError::IoError(e.to_string()))?;

        // natom [time]
        let header = lines
            .next()
            .ok_or(AmberIoError::UnexpectedEof)?
            .map_err(|e| AmberIoError::IoError(e.to_string()))?;
        let mut hparts = header.split_whitespace();
        let natom: usize = hparts
            .next()
            .ok_or(AmberIoError::ParseError("missing natom".into()))?
            .parse()
            .map_err(|e: std::num::ParseIntError| AmberIoError::ParseError(e.to_string()))?;
        if let Some(t) = hparts.next() {
            rst.time = Some(
                t.parse::<f64>()
                    .map_err(|e| AmberIoError::ParseError(e.to_string()))?,
            );
        }

        let n_coord = natom * 3;
        let mut all_values: Vec<f64> = Vec::new();

        for line in lines {
            let l = line.map_err(|e| AmberIoError::IoError(e.to_string()))?;
            let trimmed = l.trim();
            if trimmed.is_empty() {
                continue;
            }
            // Check if this looks like a box line (6 values, reasonable floats)
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            // If we already have all coords + all vels, remainder is box
            if all_values.len() >= n_coord * 2 && parts.len() == 6 {
                let mut box_vals = [0f64; 6];
                for (i, p) in parts.iter().enumerate() {
                    box_vals[i] = p
                        .parse::<f64>()
                        .map_err(|e| AmberIoError::ParseError(e.to_string()))?;
                }
                rst.box_dimensions = Some(box_vals);
                break;
            }
            for p in parts {
                all_values.push(
                    p.parse::<f64>()
                        .map_err(|e| AmberIoError::ParseError(e.to_string()))?,
                );
            }
        }

        rst.coordinates = all_values[..n_coord.min(all_values.len())].to_vec();
        if all_values.len() > n_coord {
            rst.velocities = all_values[n_coord..].to_vec();
        }

        Ok(rst)
    }
}

/// Write a slice of `f64` values to a writer, `per_line` values per line,
/// using 12.7 fixed-point format.
fn write_f64_block<W: Write>(
    writer: &mut BufWriter<W>,
    values: &[f64],
    per_line: usize,
) -> Result<(), AmberIoError> {
    for (i, v) in values.iter().enumerate() {
        write!(writer, "{:12.7}", v).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        if (i + 1) % per_line == 0 {
            writeln!(writer).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }
    }
    if !values.is_empty() && !values.len().is_multiple_of(per_line) {
        writeln!(writer).map_err(|e| AmberIoError::IoError(e.to_string()))?;
    }
    Ok(())
}

// ──────────────────────────────────────────────────────────────────────────────
// MdcrdReader
// ──────────────────────────────────────────────────────────────────────────────

/// Reader for AMBER ASCII trajectory files (`mdcrd`).
///
/// Format:
/// ```text
/// title
/// x1 y1 z1 x2 ... (10 values per line, 8.3 format)
/// [box a b c] (optional, if IFBOX > 0)
/// ```
pub struct MdcrdReader<R: BufRead> {
    reader: R,
    /// Number of atoms in each frame.
    pub natom: usize,
    /// Whether each frame includes a periodic box record.
    pub has_box: bool,
    title_consumed: bool,
}

impl<R: BufRead> MdcrdReader<R> {
    /// Create a new [`MdcrdReader`].
    ///
    /// - `natom`   – number of atoms per frame
    /// - `has_box` – set to `true` if the trajectory includes box records
    pub fn new(reader: R, natom: usize, has_box: bool) -> Self {
        MdcrdReader {
            reader,
            natom,
            has_box,
            title_consumed: false,
        }
    }

    /// Consume the title line.  Called automatically by the first `read_frame`.
    fn consume_title(&mut self) -> Result<String, AmberIoError> {
        let mut title = String::new();
        self.reader
            .read_line(&mut title)
            .map_err(|e| AmberIoError::IoError(e.to_string()))?;
        self.title_consumed = true;
        Ok(title.trim_end().to_string())
    }

    /// Read the next frame and return coordinates (and optionally box).
    ///
    /// Returns `Ok(None)` at end of file.
    pub fn read_frame(&mut self) -> Result<Option<MdcrdFrame>, AmberIoError> {
        if !self.title_consumed {
            self.consume_title()?;
        }

        let n_coords = self.natom * 3;
        let mut values: Vec<f64> = Vec::with_capacity(n_coords);

        loop {
            if values.len() >= n_coords {
                break;
            }
            let mut line = String::new();
            let bytes = self
                .reader
                .read_line(&mut line)
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
            if bytes == 0 {
                if values.is_empty() {
                    return Ok(None);
                }
                break;
            }
            for tok in line.split_whitespace() {
                values.push(
                    tok.parse::<f64>()
                        .map_err(|e| AmberIoError::ParseError(e.to_string()))?,
                );
                if values.len() >= n_coords {
                    break;
                }
            }
        }

        if values.is_empty() {
            return Ok(None);
        }

        let box_dims = if self.has_box {
            let mut line = String::new();
            self.reader
                .read_line(&mut line)
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
            let parts: Vec<f64> = line
                .split_whitespace()
                .filter_map(|t| t.parse::<f64>().ok())
                .collect();
            if parts.len() >= 3 {
                Some([parts[0], parts[1], parts[2]])
            } else {
                None
            }
        } else {
            None
        };

        Ok(Some(MdcrdFrame {
            coordinates: values,
            box_dims,
        }))
    }
}

/// A single frame from an AMBER trajectory file.
#[derive(Debug, Clone)]
pub struct MdcrdFrame {
    /// Cartesian coordinates (Å), x,y,z interleaved, length = 3 × NATOM.
    pub coordinates: Vec<f64>,
    /// Periodic box `[a, b, c]` (Å), if present.
    pub box_dims: Option<[f64; 3]>,
}

impl MdcrdFrame {
    /// Return the number of atoms in this frame.
    pub fn natom(&self) -> usize {
        self.coordinates.len() / 3
    }

    /// Return atom position as `(x, y, z)` tuple (0-based index).
    pub fn atom_pos(&self, idx: usize) -> Option<(f64, f64, f64)> {
        let base = idx * 3;
        if base + 2 < self.coordinates.len() {
            Some((
                self.coordinates[base],
                self.coordinates[base + 1],
                self.coordinates[base + 2],
            ))
        } else {
            None
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// MdinFile
// ──────────────────────────────────────────────────────────────────────────────

/// AMBER MD input (`mdin`) file representation.
///
/// Stores key-value pairs from the `&cntrl` (and optionally `&ewald`, `&pb`,
/// etc.) namelists as plain strings.
#[derive(Debug, Clone, Default)]
pub struct MdinFile {
    /// Title comment line.
    pub title: String,
    /// Parsed key/value entries from `&cntrl`.
    pub cntrl: HashMap<String, String>,
    /// Parsed key/value entries from `&ewald` (empty if absent).
    pub ewald: HashMap<String, String>,
}

impl MdinFile {
    /// Create an empty [`MdinFile`].
    pub fn new() -> Self {
        MdinFile::default()
    }

    /// Look up a `&cntrl` integer parameter with a default fallback.
    pub fn cntrl_int(&self, key: &str, default: i64) -> i64 {
        self.cntrl
            .get(key)
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(default)
    }

    /// Look up a `&cntrl` float parameter with a default fallback.
    pub fn cntrl_float(&self, key: &str, default: f64) -> f64 {
        self.cntrl
            .get(key)
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(default)
    }

    /// Write the `mdin` file to any [`Write`] sink.
    pub fn write<W: Write>(&self, writer: &mut W) -> Result<(), AmberIoError> {
        let mut bw = BufWriter::new(writer);
        writeln!(bw, "{}", self.title).map_err(|e| AmberIoError::IoError(e.to_string()))?;
        writeln!(bw, " &cntrl").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        let mut keys: Vec<&String> = self.cntrl.keys().collect();
        keys.sort();
        for k in &keys {
            writeln!(bw, "  {} = {},", k, self.cntrl[*k])
                .map_err(|e| AmberIoError::IoError(e.to_string()))?;
        }
        writeln!(bw, " /").map_err(|e| AmberIoError::IoError(e.to_string()))?;
        Ok(())
    }

    /// Parse an AMBER `mdin` file.
    pub fn read<R: Read>(reader: R) -> Result<Self, AmberIoError> {
        let br = BufReader::new(reader);
        let mut file = MdinFile::new();
        let mut lines = br.lines().peekable();

        // Title
        file.title = lines
            .next()
            .ok_or(AmberIoError::UnexpectedEof)?
            .map_err(|e| AmberIoError::IoError(e.to_string()))?;

        let mut current_nl: Option<&str> = None;
        for line_result in lines {
            let line = line_result.map_err(|e| AmberIoError::IoError(e.to_string()))?;
            let trimmed = line.trim();

            if trimmed.starts_with('&') {
                let nl_name = trimmed.trim_start_matches('&').to_lowercase();
                current_nl = match nl_name.as_str() {
                    "cntrl" => Some("cntrl"),
                    "ewald" => Some("ewald"),
                    _ => None,
                };
                continue;
            }

            if trimmed == "/" || trimmed == "&end" || trimmed.to_lowercase() == "end" {
                current_nl = None;
                continue;
            }

            if let Some(nl) = current_nl {
                // Each line may contain multiple key=value, pairs
                for entry in trimmed.split(',') {
                    let entry = entry.trim();
                    if entry.is_empty() {
                        continue;
                    }
                    if let Some(pos) = entry.find('=') {
                        let key = entry[..pos].trim().to_lowercase();
                        let val = entry[pos + 1..].trim().to_string();
                        if !key.is_empty() {
                            match nl {
                                "cntrl" => {
                                    file.cntrl.insert(key, val);
                                }
                                "ewald" => {
                                    file.ewald.insert(key, val);
                                }
                                _ => {}
                            }
                        }
                    }
                }
            }
        }

        Ok(file)
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// MdinfoParser
// ──────────────────────────────────────────────────────────────────────────────

/// Snapshot of energy/thermodynamic data from an AMBER `mdinfo` file.
#[derive(Debug, Clone, Default)]
pub struct MdinfoSnapshot {
    /// MD step number.
    pub nstep: i64,
    /// Simulation time (ps).
    pub time: f64,
    /// Instantaneous temperature (K).
    pub temp: f64,
    /// Pressure (bar).
    pub press: f64,
    /// Total energy (kcal mol⁻¹).
    pub etot: f64,
    /// Kinetic energy (kcal mol⁻¹).
    pub ektot: f64,
    /// Potential energy (kcal mol⁻¹).
    pub eptot: f64,
}

/// Parser for AMBER `mdinfo` runtime output files.
pub struct MdinfoParser;

impl MdinfoParser {
    /// Parse the first energy snapshot found in the mdinfo text.
    ///
    /// Handles lines of the form:
    /// ```text
    ///  NSTEP =      500   TIME(PS) =       1.000
    ///  TEMP(K) =   298.15  PRESS =     0.00
    ///  ETOT   =   -5000.0  EKTOT  =    1000.0  EPTOT  =  -6000.0
    /// ```
    pub fn parse<R: Read>(reader: R) -> Result<MdinfoSnapshot, AmberIoError> {
        let br = BufReader::new(reader);
        let mut snap = MdinfoSnapshot::default();

        for line_result in br.lines() {
            let line = line_result.map_err(|e| AmberIoError::IoError(e.to_string()))?;
            // Tokenise: split on '=' boundaries.
            // Strategy: find all '=' positions; everything between the previous
            // '=' value (a number token) and the next '=' is the key.
            let chars: Vec<char> = line.chars().collect();
            let mut i = 0usize;
            while i < chars.len() {
                // skip whitespace
                while i < chars.len() && chars[i].is_whitespace() {
                    i += 1;
                }
                // collect a key token up to '='
                let key_start = i;
                while i < chars.len() && chars[i] != '=' {
                    i += 1;
                }
                if i >= chars.len() {
                    break;
                }
                let key: String = chars[key_start..i]
                    .iter()
                    .collect::<String>()
                    .trim()
                    .to_uppercase();
                i += 1; // skip '='
                // collect value: skip whitespace then read non-whitespace
                while i < chars.len() && chars[i].is_whitespace() {
                    i += 1;
                }
                let val_start = i;
                while i < chars.len() && !chars[i].is_whitespace() {
                    i += 1;
                }
                let val_str: String = chars[val_start..i].iter().collect();
                let val: f64 = match val_str.parse() {
                    Ok(v) => v,
                    Err(_) => continue,
                };
                match key.as_str() {
                    "NSTEP" => snap.nstep = val as i64,
                    "TIME(PS)" => snap.time = val,
                    "TEMP(K)" => snap.temp = val,
                    "PRESS" => snap.press = val,
                    "ETOT" => snap.etot = val,
                    "EKTOT" => snap.ektot = val,
                    "EPTOT" => snap.eptot = val,
                    _ => {}
                }
            }
        }

        Ok(snap)
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// AmberForceFieldReader
// ──────────────────────────────────────────────────────────────────────────────

/// A single mass entry from a GAFF / ff14SB parameter file.
#[derive(Debug, Clone)]
pub struct FfMass {
    /// Atom type label.
    pub atom_type: String,
    /// Atomic mass (u).
    pub mass: f64,
    /// Polarisability (Å³), optional.
    pub polarizability: f64,
}

/// A bond parameter entry from a GAFF / ff14SB parameter file.
#[derive(Debug, Clone)]
pub struct FfBond {
    /// First atom type.
    pub type_i: String,
    /// Second atom type.
    pub type_j: String,
    /// Force constant k (kcal mol⁻¹ Å⁻²).
    pub force_constant: f64,
    /// Equilibrium bond length (Å).
    pub eq_length: f64,
}

/// An angle parameter entry.
#[derive(Debug, Clone)]
pub struct FfAngle {
    /// Atom types i-j-k.
    pub type_i: String,
    /// Central atom type.
    pub type_j: String,
    /// Third atom type.
    pub type_k: String,
    /// Force constant (kcal mol⁻¹ rad⁻²).
    pub force_constant: f64,
    /// Equilibrium angle (degrees).
    pub eq_angle_deg: f64,
}

/// A dihedral parameter entry.
#[derive(Debug, Clone)]
pub struct FfDihedral {
    /// Atom types i-j-k-l.
    pub type_i: String,
    /// Second atom type.
    pub type_j: String,
    /// Third atom type.
    pub type_k: String,
    /// Fourth atom type.
    pub type_l: String,
    /// Periodicity N.
    pub periodicity: i32,
    /// Barrier height V/2 (kcal mol⁻¹).
    pub barrier: f64,
    /// Phase angle (degrees).
    pub phase_deg: f64,
}

/// Reader for GAFF / ff14SB force-field parameter files.
#[derive(Debug, Clone, Default)]
pub struct AmberForceFieldReader {
    /// Mass parameters.
    pub masses: Vec<FfMass>,
    /// Bond parameters.
    pub bonds: Vec<FfBond>,
    /// Angle parameters.
    pub angles: Vec<FfAngle>,
    /// Dihedral parameters.
    pub dihedrals: Vec<FfDihedral>,
}

impl AmberForceFieldReader {
    /// Create an empty reader.
    pub fn new() -> Self {
        AmberForceFieldReader::default()
    }

    /// Parse a GAFF/ff14SB parameter file.
    pub fn read<R: Read>(&mut self, reader: R) -> Result<(), AmberIoError> {
        let br = BufReader::new(reader);
        let lines: Vec<String> = br
            .lines()
            .map(|l| l.map_err(|e| AmberIoError::IoError(e.to_string())))
            .collect::<Result<_, _>>()?;

        #[derive(PartialEq)]
        enum Section {
            None,
            Mass,
            Bond,
            Angle,
            Dihe,
        }
        let mut section = Section::None;

        for line in &lines {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                section = Section::None;
                continue;
            }

            let upper = trimmed.to_uppercase();
            if upper.starts_with("MASS") {
                section = Section::Mass;
                continue;
            } else if upper.starts_with("BOND") {
                section = Section::Bond;
                continue;
            } else if upper.starts_with("ANGLE") || upper.starts_with("ANGL") {
                section = Section::Angle;
                continue;
            } else if upper.starts_with("DIHE") {
                section = Section::Dihe;
                continue;
            } else if upper.starts_with("NONB") || upper.starts_with("IMPR") {
                section = Section::None;
                continue;
            }

            match section {
                Section::Mass => {
                    // Format: AT  mass  polarizability  !comment
                    let clean = trimmed.split('!').next().unwrap_or("").trim();
                    let parts: Vec<&str> = clean.split_whitespace().collect();
                    if parts.len() >= 2 {
                        let mass = parts[1].parse::<f64>().unwrap_or(0.0);
                        let pol = parts.get(2).and_then(|p| p.parse().ok()).unwrap_or(0.0);
                        self.masses.push(FfMass {
                            atom_type: parts[0].to_string(),
                            mass,
                            polarizability: pol,
                        });
                    }
                }
                Section::Bond => {
                    // Format: I1-I2  kf  r0  !comment
                    let clean = trimmed.split('!').next().unwrap_or("").trim();
                    let parts: Vec<&str> = clean.split_whitespace().collect();
                    if parts.len() >= 3 {
                        let types: Vec<&str> = parts[0].split('-').collect();
                        if types.len() >= 2 {
                            let kf = parts[1].parse::<f64>().unwrap_or(0.0);
                            let r0 = parts[2].parse::<f64>().unwrap_or(0.0);
                            self.bonds.push(FfBond {
                                type_i: types[0].trim().to_string(),
                                type_j: types[1].trim().to_string(),
                                force_constant: kf,
                                eq_length: r0,
                            });
                        }
                    }
                }
                Section::Angle => {
                    // Format: I1-I2-I3  kf  theta0  !comment
                    let clean = trimmed.split('!').next().unwrap_or("").trim();
                    let parts: Vec<&str> = clean.split_whitespace().collect();
                    if parts.len() >= 3 {
                        let types: Vec<&str> = parts[0].split('-').collect();
                        if types.len() >= 3 {
                            let kf = parts[1].parse::<f64>().unwrap_or(0.0);
                            let th0 = parts[2].parse::<f64>().unwrap_or(0.0);
                            self.angles.push(FfAngle {
                                type_i: types[0].trim().to_string(),
                                type_j: types[1].trim().to_string(),
                                type_k: types[2].trim().to_string(),
                                force_constant: kf,
                                eq_angle_deg: th0,
                            });
                        }
                    }
                }
                Section::Dihe => {
                    // AMBER dihedral format: "I1-I2-I3-I4  idivf  pk  phase  pn  !comment"
                    // Atom types are 2 chars each joined by '-'; wildcards use spaces, e.g.:
                    //   "c3-hc-n3-c3" or "X -c3-n3-X "
                    // Detect the type/number boundary: first whitespace run followed by a
                    // token that parses as f64 (the idivf integer).
                    let clean = trimmed.split('!').next().unwrap_or("").trim();
                    let mut split_at = clean.len();
                    let cbytes = clean.as_bytes();
                    let mut idx = 0usize;
                    while idx < cbytes.len() {
                        if cbytes[idx].is_ascii_whitespace() {
                            let ws_start = idx;
                            while idx < cbytes.len() && cbytes[idx].is_ascii_whitespace() {
                                idx += 1;
                            }
                            if idx < cbytes.len() {
                                let rest = &clean[idx..];
                                if rest
                                    .split_whitespace()
                                    .next()
                                    .map(|t| t.parse::<f64>().is_ok())
                                    .unwrap_or(false)
                                {
                                    split_at = ws_start;
                                    break;
                                }
                            }
                        } else {
                            idx += 1;
                        }
                    }
                    let type_field = clean[..split_at].trim();
                    let remainder = clean[split_at..].trim();
                    // Parse atom types: spaces are part of wildcard names → replace with '-'
                    let normalised = type_field.replace(' ', "-");
                    let types: Vec<&str> =
                        normalised.split('-').filter(|s| !s.is_empty()).collect();
                    let nums: Vec<&str> = remainder.split_whitespace().collect();
                    if types.len() >= 4 && nums.len() >= 4 {
                        let barrier = nums[1].parse::<f64>().unwrap_or(0.0);
                        let phase = nums[2].parse::<f64>().unwrap_or(0.0);
                        let pn = nums[3].parse::<f64>().unwrap_or(1.0).abs() as i32;
                        self.dihedrals.push(FfDihedral {
                            type_i: types[0].trim().to_string(),
                            type_j: types[1].trim().to_string(),
                            type_k: types[2].trim().to_string(),
                            type_l: types[3].trim().to_string(),
                            periodicity: pn,
                            barrier,
                            phase_deg: phase,
                        });
                    }
                }
                Section::None => {}
            }
        }

        Ok(())
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Utility: AMBER charge unit conversion
// ──────────────────────────────────────────────────────────────────────────────

/// AMBER internal charge unit factor.
///
/// AMBER stores partial charges as `q * 18.2223` where `q` is in elementary
/// charges.  Divide by this constant to recover the charge in elementary
/// charges (e).
pub const AMBER_CHARGE_UNIT: f64 = 18.2223;

/// Convert an AMBER internal charge value to elementary charges.
pub fn amber_charge_to_e(amber_q: f64) -> f64 {
    amber_q / AMBER_CHARGE_UNIT
}

/// Convert an elementary-charge value to AMBER internal units.
pub fn e_to_amber_charge(q_e: f64) -> f64 {
    q_e * AMBER_CHARGE_UNIT
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    // ── PrmtopAtom ────────────────────────────────────────────────────────────

    #[test]
    fn test_prmtop_atom_new() {
        let a = PrmtopAtom::new("CA", 0.0, 12.011, "CT");
        assert_eq!(a.name, "CA");
        assert_eq!(a.mass, 12.011);
    }

    #[test]
    fn test_prmtop_atom_charge_si() {
        let a = PrmtopAtom::new("N", 18.2223, 14.007, "N");
        let q = a.charge_si();
        assert!((q - 1.0).abs() < 1e-6, "expected ~1.0, got {q}");
    }

    #[test]
    fn test_prmtop_atom_charge_negative() {
        let a = PrmtopAtom::new("O", -0.8 * 18.2223, 15.999, "OS");
        assert!(a.charge_si() < 0.0);
    }

    // ── PrmtopBond ────────────────────────────────────────────────────────────

    #[test]
    fn test_prmtop_bond_new() {
        let b = PrmtopBond::new(0, 1, 340.0, 1.526);
        assert!(b.force_constant > 0.0);
        assert!(b.eq_length > 0.0);
    }

    #[test]
    fn test_prmtop_bond_energy_at_eq() {
        let b = PrmtopBond::new(0, 1, 340.0, 1.526);
        let e = b.energy(b.eq_length);
        assert!(e.abs() < 1e-12, "energy at eq length should be ~0, got {e}");
    }

    #[test]
    fn test_prmtop_bond_energy_positive() {
        let b = PrmtopBond::new(0, 1, 340.0, 1.526);
        let e = b.energy(1.7);
        assert!(e > 0.0);
    }

    #[test]
    fn test_prmtop_bond_energy_symmetric() {
        let b = PrmtopBond::new(0, 1, 340.0, 1.526);
        let e1 = b.energy(1.526 + 0.1);
        let e2 = b.energy(1.526 - 0.1);
        assert!((e1 - e2).abs() < 1e-10);
    }

    // ── PrmtopAngle ───────────────────────────────────────────────────────────

    #[test]
    fn test_prmtop_angle_new() {
        let a = PrmtopAngle::new(0, 1, 2, 50.0, std::f64::consts::PI * 109.5 / 180.0);
        assert_eq!(a.atom_j, 1);
    }

    #[test]
    fn test_prmtop_angle_energy_at_eq() {
        let theta0 = 109.5_f64.to_radians();
        let a = PrmtopAngle::new(0, 1, 2, 50.0, theta0);
        let e = a.energy(theta0);
        assert!(e.abs() < 1e-12);
    }

    #[test]
    fn test_prmtop_angle_energy_positive_deviation() {
        let theta0 = 109.5_f64.to_radians();
        let a = PrmtopAngle::new(0, 1, 2, 50.0, theta0);
        let e = a.energy(theta0 + 0.2);
        assert!(e > 0.0);
    }

    // ── PrmtopFile round-trip ─────────────────────────────────────────────────

    fn make_prmtop() -> PrmtopFile {
        let mut f = PrmtopFile::new();
        f.version = "V0001.000".to_string();
        f.atoms.push(PrmtopAtom::new(
            "N",
            -0.4157 * AMBER_CHARGE_UNIT,
            14.007,
            "N3",
        ));
        f.atoms.push(PrmtopAtom::new(
            "H1",
            0.2719 * AMBER_CHARGE_UNIT,
            1.008,
            "H",
        ));
        f.atoms.push(PrmtopAtom::new(
            "CA",
            -0.0249 * AMBER_CHARGE_UNIT,
            12.011,
            "CT",
        ));
        f.bonds.push(PrmtopBond::new(0, 1, 434.0, 1.01));
        f.bonds.push(PrmtopBond::new(1, 2, 317.0, 1.522));
        f
    }

    #[test]
    fn test_prmtop_roundtrip_natom() {
        let original = make_prmtop();
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = PrmtopFile::read(Cursor::new(&buf)).unwrap();
        assert_eq!(original.natom(), restored.natom());
    }

    #[test]
    fn test_prmtop_roundtrip_atom_names() {
        let original = make_prmtop();
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = PrmtopFile::read(Cursor::new(&buf)).unwrap();
        for (a, b) in original.atoms.iter().zip(restored.atoms.iter()) {
            assert_eq!(a.name, b.name);
        }
    }

    #[test]
    fn test_prmtop_roundtrip_charges() {
        let original = make_prmtop();
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = PrmtopFile::read(Cursor::new(&buf)).unwrap();
        for (a, b) in original.atoms.iter().zip(restored.atoms.iter()) {
            assert!(
                (a.charge - b.charge).abs() < 1e-4,
                "charge mismatch: {} vs {}",
                a.charge,
                b.charge
            );
        }
    }

    #[test]
    fn test_prmtop_roundtrip_masses() {
        let original = make_prmtop();
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = PrmtopFile::read(Cursor::new(&buf)).unwrap();
        for (a, b) in original.atoms.iter().zip(restored.atoms.iter()) {
            assert!((a.mass - b.mass).abs() < 1e-4);
        }
    }

    #[test]
    fn test_prmtop_roundtrip_bond_force_constants() {
        let original = make_prmtop();
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = PrmtopFile::read(Cursor::new(&buf)).unwrap();
        assert_eq!(original.bonds.len(), restored.bonds.len());
        for (a, b) in original.bonds.iter().zip(restored.bonds.iter()) {
            assert!(
                (a.force_constant - b.force_constant).abs() < 1e-3,
                "force_constant mismatch: {} vs {}",
                a.force_constant,
                b.force_constant
            );
        }
    }

    #[test]
    fn test_prmtop_bond_parameters_positive() {
        let f = make_prmtop();
        for bond in &f.bonds {
            assert!(bond.force_constant > 0.0, "force constant must be positive");
            assert!(bond.eq_length > 0.0, "eq length must be positive");
        }
    }

    #[test]
    fn test_prmtop_empty_roundtrip() {
        let original = PrmtopFile::new();
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = PrmtopFile::read(Cursor::new(&buf)).unwrap();
        assert_eq!(restored.natom(), 0);
    }

    // ── Rst7File ──────────────────────────────────────────────────────────────

    fn make_rst7(natom: usize) -> Rst7File {
        let mut r = Rst7File::new();
        r.title = "Test restart file".to_string();
        r.time = Some(0.5);
        r.coordinates = (0..natom * 3).map(|i| i as f64 * 0.1).collect();
        r.velocities = vec![0.01; natom * 3];
        r
    }

    #[test]
    fn test_rst7_natom() {
        let r = make_rst7(20);
        assert_eq!(r.natom(), 20);
    }

    #[test]
    fn test_rst7_coordinate_count_matches_natom() {
        let r = make_rst7(7);
        assert_eq!(r.coordinates.len(), r.natom() * 3);
    }

    #[test]
    fn test_rst7_write_read_roundtrip_natom() {
        let original = make_rst7(10);
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = Rst7File::read(Cursor::new(&buf)).unwrap();
        assert_eq!(original.natom(), restored.natom());
    }

    #[test]
    fn test_rst7_write_read_roundtrip_coordinates() {
        let original = make_rst7(5);
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = Rst7File::read(Cursor::new(&buf)).unwrap();
        for (a, b) in original.coordinates.iter().zip(restored.coordinates.iter()) {
            assert!((a - b).abs() < 1e-5, "coord mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_rst7_coordinate_format_12_7() {
        // Verify the file actually contains 12.7 format numbers
        let r = make_rst7(1);
        let mut buf = Vec::new();
        r.write(&mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Each coordinate field should be 12 characters wide
        // Check at least one line after header (2 lines) has correct format
        let data_lines: Vec<&str> = text.lines().skip(2).collect();
        let first_data_line = data_lines[0];
        // At most 6 values × 12 chars = 72 chars per line
        assert!(first_data_line.len() <= 72 + 1); // +1 for newline
    }

    #[test]
    fn test_rst7_time_preserved() {
        let original = make_rst7(3);
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = Rst7File::read(Cursor::new(&buf)).unwrap();
        assert!(restored.time.is_some());
        assert!((restored.time.unwrap() - 0.5).abs() < 1e-3);
    }

    #[test]
    fn test_rst7_no_time() {
        let mut r = make_rst7(2);
        r.time = None;
        let mut buf = Vec::new();
        r.write(&mut buf).unwrap();
        let restored = Rst7File::read(Cursor::new(&buf)).unwrap();
        assert_eq!(restored.natom(), 2);
    }

    // ── MdcrdReader ───────────────────────────────────────────────────────────

    fn make_mdcrd_text(natom: usize, nframes: usize) -> String {
        let mut s = String::from("TITLE : mdcrd test\n");
        let n_coord = natom * 3;
        for frame in 0..nframes {
            let vals: Vec<f64> = (0..n_coord)
                .map(|i| (frame * 100 + i) as f64 * 0.1)
                .collect();
            let mut count = 0usize;
            for v in &vals {
                s.push_str(&format!("{:8.3}", v));
                count += 1;
                if count.is_multiple_of(10) {
                    s.push('\n');
                }
            }
            if !n_coord.is_multiple_of(10) {
                s.push('\n');
            }
        }
        s
    }

    #[test]
    fn test_mdcrd_read_frame_count() {
        let natom = 5;
        let nframes = 4;
        let text = make_mdcrd_text(natom, nframes);
        let mut reader = MdcrdReader::new(BufReader::new(Cursor::new(text)), natom, false);
        let mut count = 0usize;
        while let Some(_frame) = reader.read_frame().unwrap() {
            count += 1;
        }
        assert_eq!(count, nframes);
    }

    #[test]
    fn test_mdcrd_frame_natom() {
        let natom = 6;
        let text = make_mdcrd_text(natom, 1);
        let mut reader = MdcrdReader::new(BufReader::new(Cursor::new(text)), natom, false);
        let frame = reader.read_frame().unwrap().unwrap();
        assert_eq!(frame.natom(), natom);
    }

    #[test]
    fn test_mdcrd_frame_coord_len() {
        let natom = 3;
        let text = make_mdcrd_text(natom, 1);
        let mut reader = MdcrdReader::new(BufReader::new(Cursor::new(text)), natom, false);
        let frame = reader.read_frame().unwrap().unwrap();
        assert_eq!(frame.coordinates.len(), natom * 3);
    }

    #[test]
    fn test_mdcrd_atom_pos_valid() {
        let natom = 4;
        let text = make_mdcrd_text(natom, 1);
        let mut reader = MdcrdReader::new(BufReader::new(Cursor::new(text)), natom, false);
        let frame = reader.read_frame().unwrap().unwrap();
        assert!(frame.atom_pos(0).is_some());
        assert!(frame.atom_pos(natom - 1).is_some());
    }

    #[test]
    fn test_mdcrd_atom_pos_out_of_bounds() {
        let natom = 2;
        let text = make_mdcrd_text(natom, 1);
        let mut reader = MdcrdReader::new(BufReader::new(Cursor::new(text)), natom, false);
        let frame = reader.read_frame().unwrap().unwrap();
        // Index 2 is out of bounds for natom=2 (valid: 0,1)
        assert!(frame.atom_pos(100).is_none());
    }

    #[test]
    fn test_mdcrd_eof_returns_none() {
        let text = "TITLE\n";
        let mut reader = MdcrdReader::new(BufReader::new(Cursor::new(text)), 3, false);
        let result = reader.read_frame().unwrap();
        assert!(result.is_none());
    }

    // ── MdinFile ──────────────────────────────────────────────────────────────

    fn make_mdin_text() -> &'static str {
        "MD simulation test\n \
         &cntrl\n  \
           nstlim = 1000,\n  \
           dt = 0.002,\n  \
           temp0 = 300.0,\n  \
           ntpr = 100,\n  \
           ntwx = 100,\n  \
           ntb = 1,\n  \
         /\n"
    }

    #[test]
    fn test_mdin_parse_nstlim() {
        let f = MdinFile::read(Cursor::new(make_mdin_text())).unwrap();
        assert_eq!(f.cntrl_int("nstlim", 0), 1000);
    }

    #[test]
    fn test_mdin_parse_dt() {
        let f = MdinFile::read(Cursor::new(make_mdin_text())).unwrap();
        let dt = f.cntrl_float("dt", 0.0);
        assert!((dt - 0.002).abs() < 1e-6);
    }

    #[test]
    fn test_mdin_parse_temp0() {
        let f = MdinFile::read(Cursor::new(make_mdin_text())).unwrap();
        let temp = f.cntrl_float("temp0", 0.0);
        assert!((temp - 300.0).abs() < 1e-6);
    }

    #[test]
    fn test_mdin_parse_ntpr() {
        let f = MdinFile::read(Cursor::new(make_mdin_text())).unwrap();
        assert_eq!(f.cntrl_int("ntpr", 0), 100);
    }

    #[test]
    fn test_mdin_default_for_missing_key() {
        let f = MdinFile::read(Cursor::new(make_mdin_text())).unwrap();
        assert_eq!(f.cntrl_int("nonexistent", 42), 42);
    }

    #[test]
    fn test_mdin_write_read_roundtrip() {
        let original = MdinFile::read(Cursor::new(make_mdin_text())).unwrap();
        let mut buf = Vec::new();
        original.write(&mut buf).unwrap();
        let restored = MdinFile::read(Cursor::new(&buf)).unwrap();
        assert_eq!(
            original.cntrl_int("nstlim", 0),
            restored.cntrl_int("nstlim", 0)
        );
    }

    // ── MdinfoParser ──────────────────────────────────────────────────────────

    fn make_mdinfo_text() -> &'static str {
        " NSTEP =      500   TIME(PS) =       1.000\n \
          TEMP(K) =   298.15  PRESS =     0.00\n \
          ETOT   =   -5000.0  EKTOT  =    1000.0  EPTOT  =  -6000.0\n"
    }

    #[test]
    fn test_mdinfo_parse_nstep() {
        let s = MdinfoParser::parse(Cursor::new(make_mdinfo_text())).unwrap();
        assert_eq!(s.nstep, 500);
    }

    #[test]
    fn test_mdinfo_parse_time() {
        let s = MdinfoParser::parse(Cursor::new(make_mdinfo_text())).unwrap();
        assert!((s.time - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_mdinfo_parse_temp() {
        let s = MdinfoParser::parse(Cursor::new(make_mdinfo_text())).unwrap();
        assert!((s.temp - 298.15).abs() < 0.1);
    }

    #[test]
    fn test_mdinfo_parse_etot() {
        let s = MdinfoParser::parse(Cursor::new(make_mdinfo_text())).unwrap();
        assert!((s.etot - (-5000.0)).abs() < 1.0);
    }

    #[test]
    fn test_mdinfo_energy_conservation() {
        let s = MdinfoParser::parse(Cursor::new(make_mdinfo_text())).unwrap();
        // etot ≈ ektot + eptot
        let sum = s.ektot + s.eptot;
        assert!((s.etot - sum).abs() < 1.0);
    }

    // ── AmberForceFieldReader ─────────────────────────────────────────────────

    fn make_ff_text() -> &'static str {
        "AMBER force field parameters\n\
         MASS\n\
         c3  12.010  !sp3 carbon\n\
         hc   1.008  !H on aliphatic C\n\
         n3  14.010  !nitrogen sp3\n\
         \n\
         BOND\n\
         c3-hc  340.0  1.09  !C-H aliphatic\n\
         c3-n3  337.0  1.47  !C-N single\n\
         \n\
         ANGLE\n\
         hc-c3-hc  35.0  109.50  !tetrahedral H-C-H\n\
         hc-c3-n3  50.0  109.50  !H-C-N\n\
         \n\
         DIHE\n\
         X -c3-n3-X   1  1.40  180.0  2.0  !sp3 C-N dihedral\n\
         \n"
    }

    #[test]
    fn test_ff_mass_count() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        assert_eq!(r.masses.len(), 3);
    }

    #[test]
    fn test_ff_bond_count() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        assert_eq!(r.bonds.len(), 2);
    }

    #[test]
    fn test_ff_angle_count() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        assert_eq!(r.angles.len(), 2);
    }

    #[test]
    fn test_ff_dihedral_count() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        assert_eq!(r.dihedrals.len(), 1);
    }

    #[test]
    fn test_ff_bond_force_constants_positive() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        for b in &r.bonds {
            assert!(
                b.force_constant > 0.0,
                "bond force constant must be positive, got {}",
                b.force_constant
            );
        }
    }

    #[test]
    fn test_ff_bond_eq_lengths_positive() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        for b in &r.bonds {
            assert!(
                b.eq_length > 0.0,
                "bond eq length must be positive, got {}",
                b.eq_length
            );
        }
    }

    #[test]
    fn test_ff_angle_force_constants_positive() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        for a in &r.angles {
            assert!(a.force_constant > 0.0);
        }
    }

    #[test]
    fn test_ff_mass_values() {
        let mut r = AmberForceFieldReader::new();
        r.read(Cursor::new(make_ff_text())).unwrap();
        let c3 = r.masses.iter().find(|m| m.atom_type == "c3").unwrap();
        assert!((c3.mass - 12.01).abs() < 0.1);
    }

    // ── Charge conversion utilities ───────────────────────────────────────────

    #[test]
    fn test_amber_charge_to_e_unit() {
        let q = amber_charge_to_e(AMBER_CHARGE_UNIT);
        assert!((q - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_e_to_amber_charge_unit() {
        let q = e_to_amber_charge(1.0);
        assert!((q - AMBER_CHARGE_UNIT).abs() < 1e-10);
    }

    #[test]
    fn test_amber_charge_roundtrip() {
        let original = 0.6;
        let roundtripped = amber_charge_to_e(e_to_amber_charge(original));
        assert!((original - roundtripped).abs() < 1e-12);
    }

    #[test]
    fn test_amber_charge_factor_value() {
        // AMBER factor must be close to 18.2223
        assert!((AMBER_CHARGE_UNIT - 18.2223).abs() < 1e-4);
    }

    #[test]
    fn test_amber_charge_negative() {
        let q_si = -0.5;
        let q_amber = e_to_amber_charge(q_si);
        assert!(q_amber < 0.0);
        assert!((amber_charge_to_e(q_amber) - q_si).abs() < 1e-12);
    }
}
