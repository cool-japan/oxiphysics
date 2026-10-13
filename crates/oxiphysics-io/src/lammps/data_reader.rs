// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reader for LAMMPS data files (the input of `read_data` and the output of
//! `write_data`).

use super::types::LammpsAtomStyle;
use crate::{Error, Result};
use std::collections::HashMap;

/// Reader for the LAMMPS data-file format.
///
/// The first line of the file is a title and is always skipped.  The header
/// that follows supplies the atom, bond and atom-type counts and the box
/// (`xlo xhi`, `ylo yhi`, `zlo zhi` and, for triclinic cells, `xy xz yz`).  The
/// body is a sequence of sections, each opened by a keyword line:
///
/// * `Masses`, `Atoms`, `Velocities` and `Bonds` are read.
/// * Every other section (`Angles`, `Dihedrals`, `Impropers`, `Pair Coeffs`,
///   `Bond Coeffs`, ... and any section a fix defines) is skipped.
///
/// The `Atoms` columns depend on the atom style.  The style is taken from the
/// comment on the keyword line (`Atoms # full`, as written by `write_data` and by
/// [`LammpsDataWriter`](super::LammpsDataWriter)), or from
/// [`LammpsDataReader::from_data_str_with_style`].  Supported styles are
/// `atomic` (`id type x y z`), `charge` (`id type q x y z`), `molecular`, `bond`
/// and `angle` (`id mol type x y z`) and `full` (`id mol type q x y z`); trailing
/// image flags are ignored.  Without a style comment only a five- or
/// eight-column (atomic) layout can be recognised; other layouts are ambiguous
/// and rejected rather than guessed.
#[derive(Debug, Clone)]
pub struct LammpsDataReader {
    pub(super) atom_ids: Vec<u64>,
    pub(super) positions: Vec<[f64; 3]>,
    pub(super) masses: Vec<f64>,
    pub(super) types: Vec<u32>,
    pub(super) charges: Vec<f64>,
    pub(super) molecule_ids: Vec<u64>,
    pub(super) velocities: Option<Vec<[f64; 3]>>,
    pub(super) bonds: Vec<(u32, u32, u32)>,
    pub(super) box_lo: [f64; 3],
    pub(super) box_hi: [f64; 3],
    pub(super) tilt: Option<[f64; 3]>,
    pub(super) style: Option<LammpsAtomStyle>,
}

/// Column layout of one `Atoms` row.
#[derive(Debug, Clone, Copy)]
struct AtomColumns {
    style: LammpsAtomStyle,
    molecule: Option<usize>,
    atom_type: usize,
    charge: Option<usize>,
    x: usize,
}

impl AtomColumns {
    fn for_style(style: LammpsAtomStyle) -> Self {
        let (molecule, atom_type, charge, x) = match style {
            LammpsAtomStyle::Atomic => (None, 1, None, 2),
            LammpsAtomStyle::Charge => (None, 1, Some(2), 3),
            LammpsAtomStyle::Molecular => (Some(1), 2, None, 3),
            LammpsAtomStyle::Full => (Some(1), 2, Some(3), 4),
        };
        Self {
            style,
            molecule,
            atom_type,
            charge,
            x,
        }
    }

    /// Number of columns before the optional image flags.
    fn width(&self) -> usize {
        self.x + 3
    }
}

/// The section a body line belongs to.
enum Section {
    Masses,
    /// `hint` is the style fixed by the keyword line or the caller; `columns`
    /// is filled in from the first row when there is no hint.
    Atoms {
        hint: Option<LammpsAtomStyle>,
        columns: Option<AtomColumns>,
    },
    Velocities,
    Bonds,
    Skipped,
}

/// Counts and box read from the header.
#[derive(Default)]
struct Header {
    n_atoms: Option<usize>,
    n_atom_types: Option<usize>,
    n_bonds: Option<usize>,
    lo: [Option<f64>; 3],
    hi: [Option<f64>; 3],
    tilt: Option<[f64; 3]>,
}

/// Whether a non-blank line opens a section.
///
/// Data rows and header lines start with a number; a section keyword (`Atoms`,
/// `Pair Coeffs`, ...) never does.  A line ending in a header marker
/// (`xlo xhi`, `atoms`, ...) is always a header line, so a damaged one is
/// reported as such instead of silently opening a section.
fn is_keyword_line(tokens: &[&str]) -> bool {
    let non_numeric_start = tokens
        .first()
        .is_some_and(|first| first.parse::<f64>().is_err());
    non_numeric_start && !is_header_line(tokens)
}

/// Whether `tokens` end in one of the header markers read by [`Header::accept`].
fn is_header_line(tokens: &[&str]) -> bool {
    matches!(
        tokens,
        [.., "xlo", "xhi"]
            | [.., "ylo", "yhi"]
            | [.., "zlo", "zhi"]
            | [.., "xy", "xz", "yz"]
            | [_, "atoms"]
            | [_, "bonds"]
            | [_, "atom", "types"]
    )
}

fn parse_error(line_no: usize, message: impl std::fmt::Display) -> Error {
    Error::Parse(format!("LAMMPS data file, line {line_no}: {message}"))
}

fn parse_number<T: std::str::FromStr>(token: &str, what: &str, line_no: usize) -> Result<T> {
    token
        .parse::<T>()
        .map_err(|_| parse_error(line_no, format!("cannot parse {what} from '{token}'")))
}

impl Header {
    /// Interpret one header line; unknown header lines are ignored.
    fn accept(&mut self, tokens: &[&str], line_no: usize) -> Result<()> {
        match tokens {
            [n, "atoms"] => self.n_atoms = Some(parse_number(n, "the atom count", line_no)?),
            [n, "atom", "types"] => {
                self.n_atom_types = Some(parse_number(n, "the atom type count", line_no)?);
            }
            [n, "bonds"] => self.n_bonds = Some(parse_number(n, "the bond count", line_no)?),
            [lo, hi, "xlo", "xhi"] => self.set_bounds(0, lo, hi, line_no)?,
            [lo, hi, "ylo", "yhi"] => self.set_bounds(1, lo, hi, line_no)?,
            [lo, hi, "zlo", "zhi"] => self.set_bounds(2, lo, hi, line_no)?,
            [xy, xz, yz, "xy", "xz", "yz"] => {
                self.tilt = Some([
                    parse_number(xy, "the xy tilt", line_no)?,
                    parse_number(xz, "the xz tilt", line_no)?,
                    parse_number(yz, "the yz tilt", line_no)?,
                ]);
            }
            _ => {}
        }
        Ok(())
    }

    fn set_bounds(&mut self, axis: usize, lo: &str, hi: &str, line_no: usize) -> Result<()> {
        let lo: f64 = parse_number(lo, "a box lower bound", line_no)?;
        let hi: f64 = parse_number(hi, "a box upper bound", line_no)?;
        if !(lo.is_finite() && hi.is_finite() && hi > lo) {
            return Err(parse_error(
                line_no,
                format!("invalid box bounds {lo} {hi}: the upper bound must exceed the lower one"),
            ));
        }
        self.lo[axis] = Some(lo);
        self.hi[axis] = Some(hi);
        Ok(())
    }

    /// The box, or an error naming the first missing bound line.
    fn bounds(&self) -> Result<([f64; 3], [f64; 3])> {
        let mut lo = [0.0; 3];
        let mut hi = [0.0; 3];
        for (axis, name) in ["x", "y", "z"].into_iter().enumerate() {
            match (self.lo[axis], self.hi[axis]) {
                (Some(l), Some(h)) => {
                    lo[axis] = l;
                    hi[axis] = h;
                }
                _ => {
                    return Err(Error::Parse(format!(
                        "LAMMPS data file: missing '{name}lo {name}hi' box line"
                    )));
                }
            }
        }
        Ok((lo, hi))
    }
}

/// Style named by the comment of an `Atoms` keyword line.
fn style_from_comment(comment: &str, line_no: usize) -> Result<Option<LammpsAtomStyle>> {
    let Some(name) = comment.split_whitespace().next() else {
        return Ok(None);
    };
    match name.to_ascii_lowercase().as_str() {
        "atomic" => Ok(Some(LammpsAtomStyle::Atomic)),
        "charge" => Ok(Some(LammpsAtomStyle::Charge)),
        "molecular" | "bond" | "angle" => Ok(Some(LammpsAtomStyle::Molecular)),
        "full" => Ok(Some(LammpsAtomStyle::Full)),
        other => Err(parse_error(
            line_no,
            format!(
                "atom style '{other}' is not supported; supported styles are atomic, charge, \
                 molecular (bond, angle) and full"
            ),
        )),
    }
}

/// Atom style implied by the column count of the first `Atoms` row.
fn infer_columns(n_tokens: usize, line_no: usize) -> Result<AtomColumns> {
    match n_tokens {
        5 | 8 => Ok(AtomColumns::for_style(LammpsAtomStyle::Atomic)),
        _ => Err(parse_error(
            line_no,
            format!(
                "cannot infer the atom style from {n_tokens} columns; name it on the keyword \
                 line (for example 'Atoms # full') or use from_data_str_with_style"
            ),
        )),
    }
}

/// Per-atom columns collected from the `Atoms` section.
#[derive(Default)]
struct AtomTable {
    ids: Vec<u64>,
    positions: Vec<[f64; 3]>,
    types: Vec<u32>,
    charges: Vec<f64>,
    molecule_ids: Vec<u64>,
}

impl AtomTable {
    fn push_row(&mut self, tokens: &[&str], cols: &AtomColumns, line_no: usize) -> Result<()> {
        if tokens.len() < cols.width() {
            return Err(parse_error(
                line_no,
                format!(
                    "{:?} atom row needs at least {} columns, found {}",
                    cols.style,
                    cols.width(),
                    tokens.len()
                ),
            ));
        }
        let id: u64 = parse_number(tokens[0], "the atom id", line_no)?;
        let atom_type: u32 = parse_number(tokens[cols.atom_type], "the atom type", line_no)?;
        if atom_type == 0 {
            return Err(parse_error(line_no, "atom types start at 1"));
        }
        let position = [
            parse_number(tokens[cols.x], "x", line_no)?,
            parse_number(tokens[cols.x + 1], "y", line_no)?,
            parse_number(tokens[cols.x + 2], "z", line_no)?,
        ];
        if let Some(col) = cols.charge {
            self.charges
                .push(parse_number(tokens[col], "the charge", line_no)?);
        }
        if let Some(col) = cols.molecule {
            self.molecule_ids
                .push(parse_number(tokens[col], "the molecule id", line_no)?);
        }
        self.ids.push(id);
        self.types.push(atom_type);
        self.positions.push(position);
        Ok(())
    }
}

impl LammpsDataReader {
    /// Parse a LAMMPS data file from a string.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] when a box line is missing or invalid, a number
    /// cannot be parsed, an `Atoms` row is too short, the atom style is
    /// unsupported or cannot be recognised, or a declared count (atoms, bonds,
    /// atom types, velocities) does not match the rows found.
    pub fn from_data_str(data: &str) -> Result<Self> {
        Self::parse(data, None)
    }

    /// Parse a LAMMPS data file, using `style` for the `Atoms` columns instead
    /// of the comment on the `Atoms` keyword line.
    ///
    /// # Errors
    ///
    /// See [`LammpsDataReader::from_data_str`].
    pub fn from_data_str_with_style(data: &str, style: LammpsAtomStyle) -> Result<Self> {
        Self::parse(data, Some(style))
    }

    fn parse(data: &str, forced_style: Option<LammpsAtomStyle>) -> Result<Self> {
        let mut header = Header::default();
        let mut section: Option<Section> = None;
        let mut atoms = AtomTable::default();
        let mut mass_rows: Vec<(u32, f64)> = Vec::new();
        let mut velocity_rows: Vec<(u64, [f64; 3])> = Vec::new();
        let mut bonds: Vec<(u32, u32, u32)> = Vec::new();
        let mut style: Option<LammpsAtomStyle> = None;
        let mut saw_masses = false;
        // The first line is the title, whatever it says.
        for (idx, raw) in data.lines().enumerate().skip(1) {
            let line_no = idx + 1;
            let (content, comment) = match raw.split_once('#') {
                Some((content, comment)) => (content.trim(), comment.trim()),
                None => (raw.trim(), ""),
            };
            if content.is_empty() {
                continue;
            }
            let tokens: Vec<&str> = content.split_whitespace().collect();
            if is_keyword_line(&tokens) {
                section = Some(match content {
                    "Masses" => {
                        saw_masses = true;
                        Section::Masses
                    }
                    "Atoms" => Section::Atoms {
                        hint: match forced_style {
                            Some(forced) => Some(forced),
                            None => style_from_comment(comment, line_no)?,
                        },
                        columns: None,
                    },
                    "Velocities" => Section::Velocities,
                    "Bonds" => Section::Bonds,
                    _ => Section::Skipped,
                });
                continue;
            }
            match section.as_mut() {
                None => header.accept(&tokens, line_no)?,
                Some(Section::Masses) => {
                    if tokens.len() < 2 {
                        return Err(parse_error(line_no, "a Masses row is 'type mass'"));
                    }
                    let atom_type: u32 = parse_number(tokens[0], "the atom type", line_no)?;
                    let mass: f64 = parse_number(tokens[1], "the mass", line_no)?;
                    if atom_type == 0 || !(mass.is_finite() && mass > 0.0) {
                        return Err(parse_error(
                            line_no,
                            format!(
                                "invalid mass entry '{content}': types start at 1 and masses are positive"
                            ),
                        ));
                    }
                    mass_rows.push((atom_type, mass));
                }
                Some(Section::Atoms { hint, columns }) => {
                    let cols = match *columns {
                        Some(cols) => cols,
                        None => {
                            let cols = match *hint {
                                Some(hinted) => AtomColumns::for_style(hinted),
                                None => infer_columns(tokens.len(), line_no)?,
                            };
                            *columns = Some(cols);
                            style = Some(cols.style);
                            cols
                        }
                    };
                    atoms.push_row(&tokens, &cols, line_no)?;
                }
                Some(Section::Velocities) => {
                    if tokens.len() < 4 {
                        return Err(parse_error(line_no, "a Velocities row is 'id vx vy vz'"));
                    }
                    velocity_rows.push((
                        parse_number(tokens[0], "the atom id", line_no)?,
                        [
                            parse_number(tokens[1], "vx", line_no)?,
                            parse_number(tokens[2], "vy", line_no)?,
                            parse_number(tokens[3], "vz", line_no)?,
                        ],
                    ));
                }
                Some(Section::Bonds) => {
                    if tokens.len() < 4 {
                        return Err(parse_error(line_no, "a Bonds row is 'id type atom1 atom2'"));
                    }
                    bonds.push((
                        parse_number(tokens[1], "the bond type", line_no)?,
                        parse_number(tokens[2], "the first bond atom", line_no)?,
                        parse_number(tokens[3], "the second bond atom", line_no)?,
                    ));
                }
                Some(Section::Skipped) => {}
            }
        }

        let (box_lo, box_hi) = header.bounds()?;
        if let Some(expected) = header.n_atoms
            && expected != atoms.ids.len()
        {
            return Err(Error::Parse(format!(
                "LAMMPS data file: header declares {expected} atoms but the Atoms section holds {}",
                atoms.ids.len()
            )));
        }
        if let Some(expected) = header.n_bonds
            && expected != bonds.len()
        {
            return Err(Error::Parse(format!(
                "LAMMPS data file: header declares {expected} bonds but the Bonds section holds {}",
                bonds.len()
            )));
        }
        if let Some(n_types) = header.n_atom_types {
            if let Some(&bad) = atoms.types.iter().find(|&&t| t as usize > n_types) {
                return Err(Error::Parse(format!(
                    "LAMMPS data file: atom type {bad} exceeds the {n_types} declared atom types"
                )));
            }
            if saw_masses && mass_rows.len() != n_types {
                return Err(Error::Parse(format!(
                    "LAMMPS data file: header declares {n_types} atom types but the Masses \
                     section holds {} entries",
                    mass_rows.len()
                )));
            }
        }
        let velocities = if velocity_rows.is_empty() {
            None
        } else {
            Some(assign_velocities(&atoms.ids, &velocity_rows)?)
        };
        let type_mass: HashMap<u32, f64> = mass_rows.into_iter().collect();
        let masses = atoms
            .types
            .iter()
            .map(|t| type_mass.get(t).copied().unwrap_or(1.0))
            .collect();
        Ok(Self {
            atom_ids: atoms.ids,
            positions: atoms.positions,
            masses,
            types: atoms.types,
            charges: atoms.charges,
            molecule_ids: atoms.molecule_ids,
            velocities,
            bonds,
            box_lo,
            box_hi,
            tilt: header.tilt,
            style,
        })
    }

    /// Return parsed positions, in file order.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// Return per-atom masses.
    ///
    /// The mass of an atom comes from the `Masses` entry of its type; a type
    /// without an entry (or a file without a `Masses` section) is given the
    /// placeholder mass `1.0`.
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }

    /// Return per-atom type IDs (1-indexed).
    pub fn types(&self) -> &[u32] {
        &self.types
    }

    /// Return the atom ids, in file order.
    pub fn atom_ids(&self) -> &[u64] {
        &self.atom_ids
    }

    /// Return per-atom charges; empty unless the atom style has a charge column
    /// (`charge`, `full`).
    pub fn charges(&self) -> &[f64] {
        &self.charges
    }

    /// Return per-atom molecule ids; empty unless the atom style has a molecule
    /// column (`molecular`, `bond`, `angle`, `full`).
    pub fn molecule_ids(&self) -> &[u64] {
        &self.molecule_ids
    }

    /// Return per-atom velocities (aligned with [`LammpsDataReader::positions`])
    /// when the file has a `Velocities` section.
    pub fn velocities(&self) -> Option<&[[f64; 3]]> {
        self.velocities.as_deref()
    }

    /// Return the bonds as `(bond_type, atom_i, atom_j)`.
    pub fn bonds(&self) -> &[(u32, u32, u32)] {
        &self.bonds
    }

    /// Return the simulation box bounds as `(lo, hi)`.
    ///
    /// For a triclinic cell these are the `xlo xhi`, `ylo yhi`, `zlo zhi` values
    /// of the file; the tilts are in [`LammpsDataReader::tilt_factors`].
    pub fn box_bounds(&self) -> ([f64; 3], [f64; 3]) {
        (self.box_lo, self.box_hi)
    }

    /// Return the `[xy, xz, yz]` tilt factors of a triclinic cell, if the file
    /// declares them.
    pub fn tilt_factors(&self) -> Option<[f64; 3]> {
        self.tilt
    }

    /// Return the atom style the `Atoms` rows were read with, or `None` when
    /// the file has no `Atoms` section.
    pub fn atom_style(&self) -> Option<LammpsAtomStyle> {
        self.style
    }
}

/// Match `Velocities` rows to atoms by id.
fn assign_velocities(ids: &[u64], rows: &[(u64, [f64; 3])]) -> Result<Vec<[f64; 3]>> {
    if rows.len() != ids.len() {
        return Err(Error::Parse(format!(
            "LAMMPS data file: Velocities section holds {} entries for {} atoms",
            rows.len(),
            ids.len()
        )));
    }
    let index_of: HashMap<u64, usize> = ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let mut velocities: Vec<Option<[f64; 3]>> = vec![None; ids.len()];
    for (id, velocity) in rows {
        let Some(&slot) = index_of.get(id) else {
            return Err(Error::Parse(format!(
                "LAMMPS data file: Velocities refers to unknown atom id {id}"
            )));
        };
        if velocities[slot].replace(*velocity).is_some() {
            return Err(Error::Parse(format!(
                "LAMMPS data file: atom id {id} has two Velocities entries"
            )));
        }
    }
    velocities
        .into_iter()
        .map(|v| {
            v.ok_or_else(|| {
                Error::Parse("LAMMPS data file: an atom has no Velocities entry".to_string())
            })
        })
        .collect()
}

impl std::str::FromStr for LammpsDataReader {
    type Err = Error;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        Self::from_data_str(s)
    }
}
