// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Text serialisation of [`PrmtopFile`] in the AMBER `prmtop` section layout.
//!
//! The writer emits the sections that [`PrmtopFile`] models, formatted the way
//! AMBER itself writes them:
//!
//! | `%FLAG` | `%FORMAT` | content |
//! |---|---|---|
//! | `POINTERS` | `10I8` | the 31 topology counters |
//! | `ATOM_NAME` | `20a4` | atom names |
//! | `CHARGE` | `5E16.8` | charges in AMBER units |
//! | `MASS` | `5E16.8` | atomic masses |
//! | `AMBER_ATOM_TYPE` | `20a4` | atom type names |
//! | `BONDS_WITHOUT_HYDROGEN` | `10I8` | `3*i, 3*j, type` triplets |
//! | `BOND_FORCE_CONSTANT`, `BOND_EQUIL_VALUE` | `5E16.8` | deduplicated bond parameter types |
//! | `ANGLES_WITHOUT_HYDROGEN` | `10I8` | `3*i, 3*j, 3*k, type` quadruplets |
//! | `ANGLE_FORCE_CONSTANT`, `ANGLE_EQUIL_VALUE` | `5E16.8` | deduplicated angle parameter types |
//!
//! Atom pointers are three times the zero-based atom index and the parameter
//! type is a one-based index into the parameter tables, exactly as in a
//! `tleap` topology.  Every other `%FLAG` section found while reading is kept
//! verbatim in [`PrmtopFile::extra_sections`] and written back.

use super::{AmberIoError, PrmtopAngle, PrmtopAtom, PrmtopBond, PrmtopFile};
use crate::fortran_fields::{parse_int_record, parse_real_record, parse_text_record};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

/// Number of integers in the `POINTERS` section.
const POINTER_COUNT: usize = 31;
/// Position of `NATOM` in `POINTERS`.
const PTR_NATOM: usize = 0;
/// Position of `NTYPES` in `POINTERS`.
const PTR_NTYPES: usize = 1;
/// Position of `NBONH` (bonds with hydrogen) in `POINTERS`.
const PTR_NBONH: usize = 2;
/// Position of `MBONA` (bonds without hydrogen) in `POINTERS`.
const PTR_MBONA: usize = 3;
/// Position of `NTHETH` (angles with hydrogen) in `POINTERS`.
const PTR_NTHETH: usize = 4;
/// Position of `MTHETA` (angles without hydrogen) in `POINTERS`.
const PTR_MTHETA: usize = 5;
/// Position of `NBONA` in `POINTERS`.
const PTR_NBONA: usize = 12;
/// Position of `NTHETA` in `POINTERS`.
const PTR_NTHETA: usize = 13;
/// Position of `NUMBND` (distinct bond types) in `POINTERS`.
const PTR_NUMBND: usize = 15;
/// Position of `NUMANG` (distinct angle types) in `POINTERS`.
const PTR_NUMANG: usize = 16;
/// Position of `NATYP` in `POINTERS`.
const PTR_NATYP: usize = 18;

/// Column width of an `I8` field.
const INT_WIDTH: usize = 8;
/// Fields per record of a `10I8` section.
const INTS_PER_LINE: usize = 10;
/// Largest value that fits an `I8` field.
const MAX_INT_FIELD: usize = 99_999_999;
/// Column width of an `E16.8` field.
const REAL_WIDTH: usize = 16;
/// Fields per record of a `5E16.8` section.
const REALS_PER_LINE: usize = 5;
/// Column width of an `a4` field.
const TEXT_WIDTH: usize = 4;
/// Fields per record of a `20a4` section.
const TEXTS_PER_LINE: usize = 20;

const FLAG_POINTERS: &str = "POINTERS";
const FLAG_ATOM_NAME: &str = "ATOM_NAME";
const FLAG_CHARGE: &str = "CHARGE";
const FLAG_MASS: &str = "MASS";
const FLAG_ATOM_TYPE: &str = "AMBER_ATOM_TYPE";
const FLAG_BONDS_H: &str = "BONDS_INC_HYDROGEN";
const FLAG_BONDS_NO_H: &str = "BONDS_WITHOUT_HYDROGEN";
const FLAG_BOND_K: &str = "BOND_FORCE_CONSTANT";
const FLAG_BOND_R0: &str = "BOND_EQUIL_VALUE";
const FLAG_ANGLES_H: &str = "ANGLES_INC_HYDROGEN";
const FLAG_ANGLES_NO_H: &str = "ANGLES_WITHOUT_HYDROGEN";
const FLAG_ANGLE_K: &str = "ANGLE_FORCE_CONSTANT";
const FLAG_ANGLE_THETA0: &str = "ANGLE_EQUIL_VALUE";

/// Sections modelled by [`PrmtopFile`]; they never appear in `extra_sections`.
const MODELLED_FLAGS: [&str; 13] = [
    FLAG_POINTERS,
    FLAG_ATOM_NAME,
    FLAG_CHARGE,
    FLAG_MASS,
    FLAG_ATOM_TYPE,
    FLAG_BONDS_H,
    FLAG_BONDS_NO_H,
    FLAG_BOND_K,
    FLAG_BOND_R0,
    FLAG_ANGLES_H,
    FLAG_ANGLES_NO_H,
    FLAG_ANGLE_K,
    FLAG_ANGLE_THETA0,
];

fn parse_error(message: impl Into<String>) -> AmberIoError {
    AmberIoError::ParseError(message.into())
}

fn invalid_data(message: impl Into<String>) -> AmberIoError {
    AmberIoError::InvalidData(message.into())
}

// ──────────────────────────────────────────────────────────────────────────────
// Writing
// ──────────────────────────────────────────────────────────────────────────────

/// Table of distinct `(first, second)` parameter pairs.
///
/// Pairs are compared bit for bit, so two terms share a type only if their
/// parameters are identical.
#[derive(Default)]
struct TypeTable {
    pairs: Vec<(f64, f64)>,
    slots: HashMap<(u64, u64), usize>,
}

impl TypeTable {
    /// Zero-based type of `(first, second)`, inserting the pair if it is new.
    fn intern(&mut self, first: f64, second: f64) -> usize {
        let key = (first.to_bits(), second.to_bits());
        if let Some(&slot) = self.slots.get(&key) {
            return slot;
        }
        let slot = self.pairs.len();
        self.pairs.push((first, second));
        self.slots.insert(key, slot);
        slot
    }

    fn firsts(&self) -> impl Iterator<Item = f64> + '_ {
        self.pairs.iter().map(|&(first, _)| first)
    }

    fn seconds(&self) -> impl Iterator<Item = f64> + '_ {
        self.pairs.iter().map(|&(_, second)| second)
    }
}

/// Atom pointers and deduplicated parameter tables of the valence terms.
struct ValenceLayout {
    bond_records: Vec<usize>,
    bond_types: TypeTable,
    angle_records: Vec<usize>,
    angle_types: TypeTable,
}

/// Pointer (`3 * index`) of atom `index`, which must exist in a topology of
/// `natom` atoms.
fn atom_pointer(
    index: usize,
    natom: usize,
    kind: &str,
    term: usize,
) -> Result<usize, AmberIoError> {
    if index >= natom {
        return Err(invalid_data(format!(
            "{kind} {term} references atom {index}, but the topology has only {natom} atoms"
        )));
    }
    index
        .checked_mul(3)
        .ok_or_else(|| invalid_data(format!("{kind} {term}: atom index {index} is too large")))
}

impl ValenceLayout {
    fn build(file: &PrmtopFile) -> Result<Self, AmberIoError> {
        let natom = file.atoms.len();
        let mut bond_types = TypeTable::default();
        let mut bond_records = Vec::with_capacity(3 * file.bonds.len());
        for (term, bond) in file.bonds.iter().enumerate() {
            let i = atom_pointer(bond.atom_i, natom, "bond", term)?;
            let j = atom_pointer(bond.atom_j, natom, "bond", term)?;
            let bond_type = bond_types.intern(bond.force_constant, bond.eq_length);
            bond_records.extend([i, j, bond_type + 1]);
        }
        let mut angle_types = TypeTable::default();
        let mut angle_records = Vec::with_capacity(4 * file.angles.len());
        for (term, angle) in file.angles.iter().enumerate() {
            let i = atom_pointer(angle.atom_i, natom, "angle", term)?;
            let j = atom_pointer(angle.atom_j, natom, "angle", term)?;
            let k = atom_pointer(angle.atom_k, natom, "angle", term)?;
            let angle_type = angle_types.intern(angle.force_constant, angle.eq_angle);
            angle_records.extend([i, j, k, angle_type + 1]);
        }
        Ok(Self {
            bond_records,
            bond_types,
            angle_records,
            angle_types,
        })
    }
}

/// Reject atom names and types that do not fit their four-character column or
/// that the reader could not tell from the surrounding file structure.
///
/// A `%` is refused because a record that starts with one is a directive
/// (`%FLAG`, `%FORMAT`, ...), so a name such as `%FLA` followed by `G` in the
/// next column would end the section instead of naming an atom.
fn check_text_field(value: &str, what: &str, atom: usize) -> Result<(), AmberIoError> {
    if value.chars().count() > TEXT_WIDTH || value.chars().any(|c| c.is_control() || c == '%') {
        return Err(invalid_data(format!(
            "{what} `{value}` of atom {atom} must have at most {TEXT_WIDTH} characters, no control characters and no `%`"
        )));
    }
    Ok(())
}

/// Reject extra sections that would not survive a write/read round trip.
fn check_extra_sections(extra: &HashMap<String, Vec<String>>) -> Result<(), AmberIoError> {
    for (name, lines) in extra {
        if name.is_empty() || name.chars().any(|c| c.is_whitespace() || c.is_control()) {
            return Err(invalid_data(format!(
                "extra section name `{name}` must be non-empty and free of whitespace"
            )));
        }
        if MODELLED_FLAGS.contains(&name.as_str()) {
            return Err(invalid_data(format!(
                "extra section `{name}` collides with a section modelled by PrmtopFile"
            )));
        }
        for line in lines {
            let directive = line.trim_start();
            if line.contains(['\n', '\r'])
                || directive.starts_with("%FLAG")
                || directive.starts_with("%VERSION")
            {
                return Err(invalid_data(format!(
                    "extra section `{name}` contains a line that would be read as a new section"
                )));
            }
        }
    }
    Ok(())
}

/// Writes fixed-width fields, `per_line` to a record.
struct FieldWriter<'a, W: Write> {
    out: &'a mut W,
    per_line: usize,
    written: usize,
}

impl<'a, W: Write> FieldWriter<'a, W> {
    /// Start a section: write its `%FLAG` and `%FORMAT` lines.
    fn begin(
        out: &'a mut W,
        flag: &str,
        format: &str,
        per_line: usize,
    ) -> Result<Self, AmberIoError> {
        writeln!(out, "%FLAG {flag}")?;
        writeln!(out, "%FORMAT({format})")?;
        Ok(Self {
            out,
            per_line,
            written: 0,
        })
    }

    fn push(&mut self, field: &str) -> Result<(), AmberIoError> {
        self.out.write_all(field.as_bytes())?;
        self.written += 1;
        if self.written.is_multiple_of(self.per_line) {
            writeln!(self.out)?;
        }
        Ok(())
    }

    /// Terminate the last, possibly short, record.
    fn finish(self) -> Result<(), AmberIoError> {
        if !self.written.is_multiple_of(self.per_line) {
            writeln!(self.out)?;
        }
        Ok(())
    }
}

fn write_int_section<W: Write>(
    out: &mut W,
    flag: &str,
    values: &[usize],
) -> Result<(), AmberIoError> {
    let mut section = FieldWriter::begin(out, flag, "10I8", INTS_PER_LINE)?;
    for &value in values {
        if value > MAX_INT_FIELD {
            return Err(invalid_data(format!(
                "{flag}: value {value} does not fit an I{INT_WIDTH} field"
            )));
        }
        section.push(&format!("{value:INT_WIDTH$}"))?;
    }
    section.finish()
}

fn write_real_section<W: Write>(
    out: &mut W,
    flag: &str,
    values: impl IntoIterator<Item = f64>,
) -> Result<(), AmberIoError> {
    let mut section = FieldWriter::begin(out, flag, "5E16.8", REALS_PER_LINE)?;
    for value in values {
        // A finite `f64` never needs more than 16 characters with eight
        // decimals, so the field cannot overflow its column.
        section.push(&format!("{value:REAL_WIDTH$.8E}"))?;
    }
    section.finish()
}

fn write_text_section<'v, W: Write>(
    out: &mut W,
    flag: &str,
    values: impl IntoIterator<Item = &'v str>,
) -> Result<(), AmberIoError> {
    let mut section = FieldWriter::begin(out, flag, "20a4", TEXTS_PER_LINE)?;
    for value in values {
        section.push(&format!("{value:<TEXT_WIDTH$}"))?;
    }
    section.finish()
}

impl PrmtopFile {
    /// Write the topology to any [`Write`] sink in the AMBER `prmtop` layout.
    ///
    /// The modelled sections are `POINTERS`, `ATOM_NAME`, `CHARGE`, `MASS`,
    /// `AMBER_ATOM_TYPE`, the bond list (`BONDS_WITHOUT_HYDROGEN` with the
    /// deduplicated `BOND_FORCE_CONSTANT` / `BOND_EQUIL_VALUE` tables) and the
    /// angle list (`ANGLES_WITHOUT_HYDROGEN` with `ANGLE_FORCE_CONSTANT` /
    /// `ANGLE_EQUIL_VALUE`).  [`PrmtopBond`] and [`PrmtopAngle`] carry no
    /// hydrogen classification, so every term is written to the
    /// `*_WITHOUT_HYDROGEN` section and keeps its position on a round trip.
    /// The sections in [`PrmtopFile::extra_sections`] follow, sorted by name.
    ///
    /// The counters of `POINTERS` that describe the atoms, bonds and angles
    /// are derived from the topology; all others come from
    /// [`PrmtopFile::pointers`] (zero for a topology built in memory, except
    /// `NTYPES` and `NATYP`, which count the distinct atom type names).
    ///
    /// # Errors
    ///
    /// Returns [`AmberIoError::InvalidData`] when the topology cannot be
    /// represented: a bond or angle refers to a missing atom, an atom name or
    /// type is wider than four characters or contains a control character or a
    /// `%`, a value overflows its column, or an extra section would be
    /// misread.  [`AmberIoError::IoError`] reports a failing sink.
    pub fn write<W: Write>(&self, writer: &mut W) -> Result<(), AmberIoError> {
        if self.version.contains(['\n', '\r']) {
            return Err(invalid_data("the version string must be a single line"));
        }
        for (atom, entry) in self.atoms.iter().enumerate() {
            check_text_field(&entry.name, "atom name", atom)?;
            check_text_field(&entry.atom_type, "atom type", atom)?;
        }
        check_extra_sections(&self.extra_sections)?;
        let layout = ValenceLayout::build(self)?;

        // Counters that describe sections this writer does not produce
        // (residues, dihedrals, Lennard-Jones types, box, ...) keep the value
        // read from the file, so the preserved extra sections stay consistent
        // with POINTERS.  Only a topology built in memory has to estimate the
        // atom type count from the type names.
        let mut pointers = self.pointers.clone();
        if pointers.is_empty() {
            let distinct_types: HashSet<&str> =
                self.atoms.iter().map(|a| a.atom_type.as_str()).collect();
            pointers.resize(POINTER_COUNT, 0);
            pointers[PTR_NTYPES] = distinct_types.len();
            pointers[PTR_NATYP] = distinct_types.len();
        }
        if pointers.len() < POINTER_COUNT {
            pointers.resize(POINTER_COUNT, 0);
        }
        pointers[PTR_NATOM] = self.atoms.len();
        pointers[PTR_NBONH] = 0;
        pointers[PTR_MBONA] = self.bonds.len();
        pointers[PTR_NTHETH] = 0;
        pointers[PTR_MTHETA] = self.angles.len();
        pointers[PTR_NBONA] = self.bonds.len();
        pointers[PTR_NTHETA] = self.angles.len();
        pointers[PTR_NUMBND] = layout.bond_types.pairs.len();
        pointers[PTR_NUMANG] = layout.angle_types.pairs.len();

        let mut out = BufWriter::new(writer);
        writeln!(out, "%VERSION  {}", self.version)?;
        write_int_section(&mut out, FLAG_POINTERS, &pointers)?;
        write_text_section(
            &mut out,
            FLAG_ATOM_NAME,
            self.atoms.iter().map(|a| a.name.as_str()),
        )?;
        write_real_section(&mut out, FLAG_CHARGE, self.atoms.iter().map(|a| a.charge))?;
        write_real_section(&mut out, FLAG_MASS, self.atoms.iter().map(|a| a.mass))?;
        write_text_section(
            &mut out,
            FLAG_ATOM_TYPE,
            self.atoms.iter().map(|a| a.atom_type.as_str()),
        )?;
        write_int_section(&mut out, FLAG_BONDS_NO_H, &layout.bond_records)?;
        write_real_section(&mut out, FLAG_BOND_K, layout.bond_types.firsts())?;
        write_real_section(&mut out, FLAG_BOND_R0, layout.bond_types.seconds())?;
        write_int_section(&mut out, FLAG_ANGLES_NO_H, &layout.angle_records)?;
        write_real_section(&mut out, FLAG_ANGLE_K, layout.angle_types.firsts())?;
        write_real_section(&mut out, FLAG_ANGLE_THETA0, layout.angle_types.seconds())?;

        let mut extra_names: Vec<&String> = self.extra_sections.keys().collect();
        extra_names.sort();
        for name in extra_names {
            writeln!(out, "%FLAG {name}")?;
            if let Some(lines) = self.extra_sections.get(name) {
                for line in lines {
                    writeln!(out, "{line}")?;
                }
            }
        }
        out.flush()?;
        Ok(())
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Reading
// ──────────────────────────────────────────────────────────────────────────────

/// The lines of one `%FLAG` section, verbatim.
struct Section {
    lines: Vec<String>,
}

/// Column width declared by a `%FORMAT(...)` line such as `%FORMAT(10I8)`,
/// `%FORMAT(5E16.8)` or `%FORMAT(20a4)`.
fn format_width(line: &str) -> Option<usize> {
    let spec = line.trim().strip_prefix("%FORMAT")?.trim();
    let spec = spec.strip_prefix('(')?.strip_suffix(')')?;
    let mut rest = spec
        .trim_start_matches(|c: char| c.is_ascii_digit())
        .chars();
    if !rest.next()?.is_ascii_alphabetic() {
        return None;
    }
    let digits: String = rest.take_while(char::is_ascii_digit).collect();
    digits.parse().ok().filter(|&width| width > 0)
}

impl Section {
    /// Column width of the section: the declared one, or `default`.
    fn width(&self, default: usize) -> usize {
        self.lines
            .iter()
            .find_map(|line| format_width(line))
            .unwrap_or(default)
    }

    /// Data records: every non-empty line that is not a `%` directive.
    fn records(&self) -> impl Iterator<Item = &str> {
        self.lines
            .iter()
            .map(String::as_str)
            .filter(|line| !line.is_empty() && !line.trim_start().starts_with('%'))
    }
}

/// Split the file into its `%VERSION` text and its `%FLAG` sections.
fn split_sections(lines: &[String]) -> Result<(String, HashMap<String, Section>), AmberIoError> {
    fn close(
        current: Option<(String, Vec<String>)>,
        sections: &mut HashMap<String, Section>,
    ) -> Result<(), AmberIoError> {
        if let Some((name, lines)) = current
            && sections.insert(name.clone(), Section { lines }).is_some()
        {
            return Err(parse_error(format!("duplicate section %FLAG {name}")));
        }
        Ok(())
    }

    let mut version = String::new();
    let mut sections: HashMap<String, Section> = HashMap::new();
    let mut current: Option<(String, Vec<String>)> = None;
    for line in lines {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("%FLAG") {
            close(current.take(), &mut sections)?;
            let name = rest.trim();
            if name.is_empty() {
                return Err(parse_error("a %FLAG line has no section name"));
            }
            current = Some((name.to_string(), Vec::new()));
        } else if let Some(rest) = trimmed.strip_prefix("%VERSION") {
            close(current.take(), &mut sections)?;
            version = rest.trim().to_string();
        } else if let Some((_, body)) = current.as_mut() {
            body.push(line.clone());
        }
    }
    close(current.take(), &mut sections)?;
    Ok((version, sections))
}

fn take_ints(
    sections: &mut HashMap<String, Section>,
    flag: &str,
) -> Result<Option<Vec<i64>>, AmberIoError> {
    let Some(section) = sections.remove(flag) else {
        return Ok(None);
    };
    let width = section.width(INT_WIDTH);
    let mut values = Vec::new();
    for record in section.records() {
        values.extend(
            parse_int_record(record, width).map_err(|e| parse_error(format!("{flag}: {e}")))?,
        );
    }
    Ok(Some(values))
}

fn take_reals(
    sections: &mut HashMap<String, Section>,
    flag: &str,
) -> Result<Option<Vec<f64>>, AmberIoError> {
    let Some(section) = sections.remove(flag) else {
        return Ok(None);
    };
    let width = section.width(REAL_WIDTH);
    let mut values = Vec::new();
    for record in section.records() {
        values.extend(
            parse_real_record(record, width).map_err(|e| parse_error(format!("{flag}: {e}")))?,
        );
    }
    Ok(Some(values))
}

fn take_texts(sections: &mut HashMap<String, Section>, flag: &str) -> Option<Vec<String>> {
    let section = sections.remove(flag)?;
    let width = section.width(TEXT_WIDTH);
    let mut values = Vec::new();
    for record in section.records() {
        values.extend(parse_text_record(record, width));
    }
    Some(values)
}

/// Per-atom data is either absent (and defaulted) or holds exactly `natom`
/// entries.  A missing section is only an error when `required` and there are
/// atoms.
fn per_atom<T: Clone + Default>(
    values: Option<Vec<T>>,
    natom: usize,
    flag: &str,
    required: bool,
) -> Result<Vec<T>, AmberIoError> {
    match values {
        Some(values) if values.len() == natom => Ok(values),
        Some(values) => Err(parse_error(format!(
            "{flag} holds {} entries, but the topology has {natom} atoms",
            values.len()
        ))),
        None if natom == 0 || !required => Ok(vec![T::default(); natom]),
        None => Err(AmberIoError::MissingSection(flag.to_string())),
    }
}

/// Zero-based atom index designated by a pointer, if the pointer designates
/// one of the `natom` atoms.
fn decode_pointer(pointer: i64, natom: usize) -> Option<usize> {
    let pointer = u64::try_from(pointer).ok()?;
    if !pointer.is_multiple_of(3) {
        return None;
    }
    usize::try_from(pointer / 3)
        .ok()
        .filter(|&atom| atom < natom)
}

/// Decode the records of a `*_INC_HYDROGEN` / `*_WITHOUT_HYDROGEN` section
/// pair into `(atoms, zero-based type)` terms.
///
/// Each record is `N` atom pointers followed by a one-based type index in
/// `1..=n_types`.
fn decode_terms<const N: usize>(
    parts: [(&str, &[i64]); 2],
    natom: usize,
    n_types: usize,
) -> Result<Vec<([usize; N], usize)>, AmberIoError> {
    let stride = N + 1;
    let mut terms = Vec::new();
    for (flag, ints) in parts {
        if !ints.len().is_multiple_of(stride) {
            return Err(parse_error(format!(
                "{flag}: {} integers do not form whole records of {stride} values",
                ints.len()
            )));
        }
        for (index, record) in ints.chunks_exact(stride).enumerate() {
            let mut atoms = [0usize; N];
            for (slot, &pointer) in atoms.iter_mut().zip(record) {
                *slot = decode_pointer(pointer, natom).ok_or_else(|| {
                    parse_error(format!(
                        "{flag}: record {index} has atom pointer {pointer}, which does not designate one of the {natom} atoms"
                    ))
                })?;
            }
            let type_index = record.get(N).copied().unwrap_or(0);
            let term_type = usize::try_from(type_index)
                .ok()
                .filter(|t| (1..=n_types).contains(t))
                .ok_or_else(|| {
                    parse_error(format!(
                        "{flag}: record {index} has parameter type {type_index}, expected 1..={n_types}"
                    ))
                })?;
            terms.push((atoms, term_type - 1));
        }
    }
    Ok(terms)
}

/// The parameter tables and the index sections of one kind of valence term.
struct TermSections {
    kind: &'static str,
    table_flags: [&'static str; 2],
    index_flags: [&'static str; 2],
    first: Vec<f64>,
    second: Vec<f64>,
    inc_hydrogen: Option<Vec<i64>>,
    without_hydrogen: Option<Vec<i64>>,
}

impl TermSections {
    /// Remove the sections of one kind of term from `sections`.
    ///
    /// `table_flags` names the two parameter tables, `index_flags` the
    /// `*_INC_HYDROGEN` and `*_WITHOUT_HYDROGEN` index sections.
    fn take(
        sections: &mut HashMap<String, Section>,
        kind: &'static str,
        table_flags: [&'static str; 2],
        index_flags: [&'static str; 2],
    ) -> Result<Self, AmberIoError> {
        Ok(Self {
            kind,
            table_flags,
            index_flags,
            first: take_reals(sections, table_flags[0])?.unwrap_or_default(),
            second: take_reals(sections, table_flags[1])?.unwrap_or_default(),
            inc_hydrogen: take_ints(sections, index_flags[0])?,
            without_hydrogen: take_ints(sections, index_flags[1])?,
        })
    }

    /// Check that the parameter tables and the index sections belong together.
    fn check(&self) -> Result<(), AmberIoError> {
        if self.first.len() != self.second.len() {
            return Err(parse_error(format!(
                "{} holds {} entries but {} holds {}",
                self.table_flags[0],
                self.first.len(),
                self.table_flags[1],
                self.second.len()
            )));
        }
        let has_tables = !self.first.is_empty();
        if has_tables && self.inc_hydrogen.is_none() && self.without_hydrogen.is_none() {
            return Err(parse_error(format!(
                "{} parameter tables are present, but neither {} nor {} lists the {}s they apply to",
                self.kind, self.index_flags[0], self.index_flags[1], self.kind
            )));
        }
        Ok(())
    }

    /// Decode the terms of both index sections, hydrogen terms first.
    fn terms<const N: usize>(
        &self,
        natom: usize,
    ) -> Result<Vec<([usize; N], usize)>, AmberIoError> {
        decode_terms::<N>(
            [
                (
                    self.index_flags[0],
                    self.inc_hydrogen.as_deref().unwrap_or_default(),
                ),
                (
                    self.index_flags[1],
                    self.without_hydrogen.as_deref().unwrap_or_default(),
                ),
            ],
            natom,
            self.first.len(),
        )
    }
}

impl PrmtopFile {
    /// Parse a topology written by [`PrmtopFile::write`] or by `tleap`.
    ///
    /// Atoms come from `ATOM_NAME`, `CHARGE`, `MASS` and `AMBER_ATOM_TYPE`;
    /// bonds and angles are decoded from the `*_INC_HYDROGEN` and
    /// `*_WITHOUT_HYDROGEN` index sections (hydrogen terms first) and resolved
    /// through the parameter tables.  Fixed-width fields that touch each other
    /// are separated by column.  Every other `%FLAG` section is stored
    /// verbatim, starting with its `%FORMAT` line, in
    /// [`PrmtopFile::extra_sections`], and the `POINTERS` counters in
    /// [`PrmtopFile::pointers`].
    ///
    /// # Errors
    ///
    /// Returns [`AmberIoError::ParseError`] for a malformed file: a duplicate
    /// section, a negative `POINTERS` entry, a per-atom section whose length
    /// differs from `NATOM`, an atom pointer that does not designate an atom,
    /// a parameter type outside its table, a trailing incomplete record, or
    /// bond/angle parameter tables without the index section that lists the
    /// terms.  [`AmberIoError::MissingSection`] is returned when atoms are
    /// declared but `ATOM_NAME` is absent.
    pub fn read<R: Read>(reader: R) -> Result<Self, AmberIoError> {
        let lines: Vec<String> = BufReader::new(reader)
            .lines()
            .collect::<Result<_, std::io::Error>>()?;
        let (version, mut sections) = split_sections(&lines)?;
        let mut file = PrmtopFile::new();
        file.version = version;

        let pointers = take_ints(&mut sections, FLAG_POINTERS)?
            .map(|values| {
                values
                    .iter()
                    .enumerate()
                    .map(|(index, &value)| {
                        usize::try_from(value).map_err(|_| {
                            parse_error(format!("POINTERS entry {index} is negative ({value})"))
                        })
                    })
                    .collect::<Result<Vec<usize>, AmberIoError>>()
            })
            .transpose()?;
        let names = take_texts(&mut sections, FLAG_ATOM_NAME);
        let charges = take_reals(&mut sections, FLAG_CHARGE)?;
        let masses = take_reals(&mut sections, FLAG_MASS)?;
        let types = take_texts(&mut sections, FLAG_ATOM_TYPE);

        let natom = match &pointers {
            Some(values) => values
                .get(PTR_NATOM)
                .copied()
                .ok_or_else(|| parse_error("POINTERS holds no NATOM entry"))?,
            None => names.as_ref().map_or(0, Vec::len),
        };
        file.pointers = pointers.unwrap_or_default();
        let names = per_atom(names, natom, FLAG_ATOM_NAME, true)?;
        let charges = per_atom(charges, natom, FLAG_CHARGE, false)?;
        let masses = per_atom(masses, natom, FLAG_MASS, false)?;
        let types = per_atom(types, natom, FLAG_ATOM_TYPE, false)?;
        file.atoms = names
            .into_iter()
            .zip(charges)
            .zip(masses)
            .zip(types)
            .map(|(((name, charge), mass), atom_type)| {
                PrmtopAtom::new(name, charge, mass, atom_type)
            })
            .collect();

        let bond_sections = TermSections::take(
            &mut sections,
            "bond",
            [FLAG_BOND_K, FLAG_BOND_R0],
            [FLAG_BONDS_H, FLAG_BONDS_NO_H],
        )?;
        bond_sections.check()?;
        for ([atom_i, atom_j], bond_type) in bond_sections.terms::<2>(natom)? {
            file.bonds.push(PrmtopBond::new(
                atom_i,
                atom_j,
                bond_sections.first[bond_type],
                bond_sections.second[bond_type],
            ));
        }

        let angle_sections = TermSections::take(
            &mut sections,
            "angle",
            [FLAG_ANGLE_K, FLAG_ANGLE_THETA0],
            [FLAG_ANGLES_H, FLAG_ANGLES_NO_H],
        )?;
        angle_sections.check()?;
        for ([atom_i, atom_j, atom_k], angle_type) in angle_sections.terms::<3>(natom)? {
            file.angles.push(PrmtopAngle::new(
                atom_i,
                atom_j,
                atom_k,
                angle_sections.first[angle_type],
                angle_sections.second[angle_type],
            ));
        }

        for (name, section) in sections {
            file.extra_sections.insert(name, section.lines);
        }
        Ok(file)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn water() -> PrmtopFile {
        let mut file = PrmtopFile::new();
        file.version = "VERSION_STAMP = V0001.000".to_string();
        file.atoms.push(PrmtopAtom::new("O", -15.1977, 16.0, "OW"));
        file.atoms.push(PrmtopAtom::new("H1", 7.59885, 1.008, "HW"));
        file.atoms.push(PrmtopAtom::new("H2", 7.59885, 1.008, "HW"));
        file.bonds.push(PrmtopBond::new(0, 1, 553.0, 0.9572));
        file.bonds.push(PrmtopBond::new(0, 2, 553.0, 0.9572));
        file.angles.push(PrmtopAngle::new(1, 0, 2, 100.0, 1.824));
        file
    }

    fn write_to_string(file: &PrmtopFile) -> String {
        let mut buf = Vec::new();
        file.write(&mut buf).unwrap();
        String::from_utf8(buf).unwrap()
    }

    fn read_str(text: &str) -> Result<PrmtopFile, AmberIoError> {
        PrmtopFile::read(Cursor::new(text.as_bytes()))
    }

    /// Data lines of section `flag` (everything between its `%FORMAT` line and
    /// the next `%FLAG`).
    fn section_lines(text: &str, flag: &str) -> Vec<String> {
        text.lines()
            .skip_while(|line| *line != format!("%FLAG {flag}"))
            .skip(2)
            .take_while(|line| !line.starts_with('%'))
            .map(str::to_string)
            .collect()
    }

    fn int_lines(values: &[i64]) -> String {
        values
            .chunks(10)
            .map(|row| row.iter().map(|v| format!("{v:8}")).collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    }

    fn real_lines(values: &[f64]) -> String {
        values
            .chunks(5)
            .map(|row| row.iter().map(|v| format!("{v:16.8E}")).collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    }

    fn text_lines(values: &[&str]) -> String {
        values
            .chunks(20)
            .map(|row| row.iter().map(|v| format!("{v:<4}")).collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    }

    fn flag(name: &str, format: &str, body: &str) -> String {
        format!("%FLAG {name}\n%FORMAT({format})\n{body}\n")
    }

    /// Four-atom topology text with hydrogen and heavy-atom bond and angle
    /// sections.
    fn four_atom_text() -> String {
        let mut text = String::from("%VERSION  VERSION_STAMP = V0001.000\n");
        text.push_str(&flag("POINTERS", "10I8", &int_lines(&[4, 3])));
        text.push_str(&flag(
            "ATOM_NAME",
            "20a4",
            &text_lines(&["C1", "H1", "C2", "H2"]),
        ));
        text.push_str(&flag(
            "CHARGE",
            "5E16.8",
            &real_lines(&[1.0, -1.0, 2.0, -2.0]),
        ));
        text.push_str(&flag(
            "MASS",
            "5E16.8",
            &real_lines(&[12.01, 1.008, 12.01, 1.008]),
        ));
        text.push_str(&flag(
            "AMBER_ATOM_TYPE",
            "20a4",
            &text_lines(&["CT", "HC", "CT", "HC"]),
        ));
        text.push_str(&flag("BONDS_INC_HYDROGEN", "10I8", &int_lines(&[0, 3, 1])));
        text.push_str(&flag(
            "BONDS_WITHOUT_HYDROGEN",
            "10I8",
            &int_lines(&[6, 9, 2]),
        ));
        text.push_str(&flag(
            "BOND_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[340.0, 310.0]),
        ));
        text.push_str(&flag(
            "BOND_EQUIL_VALUE",
            "5E16.8",
            &real_lines(&[1.09, 1.526]),
        ));
        text.push_str(&flag(
            "ANGLES_INC_HYDROGEN",
            "10I8",
            &int_lines(&[3, 0, 6, 2]),
        ));
        text.push_str(&flag(
            "ANGLES_WITHOUT_HYDROGEN",
            "10I8",
            &int_lines(&[0, 6, 9, 1]),
        ));
        text.push_str(&flag(
            "ANGLE_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[50.0, 63.0]),
        ));
        text.push_str(&flag(
            "ANGLE_EQUIL_VALUE",
            "5E16.8",
            &real_lines(&[1.91, 1.94]),
        ));
        text
    }

    // ── writing ───────────────────────────────────────────────────────────────

    #[test]
    fn water_round_trips_exactly() {
        let original = water();
        let restored = read_str(&write_to_string(&original)).unwrap();
        assert_eq!(restored.version, original.version);
        assert_eq!(restored.atoms, original.atoms);
        assert_eq!(restored.bonds, original.bonds);
        assert_eq!(restored.angles, original.angles);
        assert!(restored.extra_sections.is_empty());
        let bond_pairs: Vec<(usize, usize)> = restored
            .bonds
            .iter()
            .map(|b| (b.atom_i, b.atom_j))
            .collect();
        assert_eq!(bond_pairs, vec![(0, 1), (0, 2)]);
        let angle = &restored.angles[0];
        assert_eq!((angle.atom_i, angle.atom_j, angle.atom_k), (1, 0, 2));
    }

    #[test]
    fn written_text_lists_atom_pointers_and_deduplicated_tables() {
        let text = write_to_string(&water());
        assert!(text.starts_with("%VERSION  VERSION_STAMP = V0001.000\n"));
        assert_eq!(
            section_lines(&text, "BONDS_WITHOUT_HYDROGEN"),
            vec!["       0       3       1       0       6       1"]
        );
        assert_eq!(
            section_lines(&text, "ANGLES_WITHOUT_HYDROGEN"),
            vec!["       3       0       6       1"]
        );
        let bond_k = section_lines(&text, "BOND_FORCE_CONSTANT");
        assert_eq!(bond_k.len(), 1);
        assert_eq!(bond_k[0].trim(), "5.53000000E2");
        let theta0 = section_lines(&text, "ANGLE_EQUIL_VALUE");
        assert_eq!(theta0[0].trim(), "1.82400000E0");
        assert_eq!(section_lines(&text, "ATOM_NAME"), vec!["O   H1  H2  "]);
        assert_eq!(
            section_lines(&text, "AMBER_ATOM_TYPE"),
            vec!["OW  HW  HW  "]
        );
        assert!(
            !text.contains("ATOM_TYPE_INDEX"),
            "string atom types must not be stored under the integer flag"
        );
    }

    #[test]
    fn pointers_section_carries_the_counters() {
        let text = write_to_string(&water());
        let pointers: Vec<i64> = section_lines(&text, "POINTERS")
            .iter()
            .flat_map(|line| line.split_whitespace().map(|t| t.parse::<i64>().unwrap()))
            .collect();
        assert_eq!(pointers.len(), POINTER_COUNT);
        assert_eq!(pointers[PTR_NATOM], 3);
        assert_eq!(pointers[PTR_NTYPES], 2);
        assert_eq!(pointers[PTR_MBONA], 2);
        assert_eq!(pointers[PTR_MTHETA], 1);
        assert_eq!(pointers[PTR_NUMBND], 1);
        assert_eq!(pointers[PTR_NUMANG], 1);
        assert_eq!(pointers[PTR_NATYP], 2);
    }

    #[test]
    fn distinct_parameters_get_distinct_types() {
        let mut file = PrmtopFile::new();
        for name in ["C1", "C2", "C3", "C4"] {
            file.atoms.push(PrmtopAtom::new(name, 0.0, 12.011, "CT"));
        }
        file.bonds.push(PrmtopBond::new(0, 1, 310.0, 1.526));
        file.bonds.push(PrmtopBond::new(1, 2, 340.0, 1.09));
        file.bonds.push(PrmtopBond::new(2, 3, 310.0, 1.526));
        let text = write_to_string(&file);
        assert_eq!(
            section_lines(&text, "BONDS_WITHOUT_HYDROGEN"),
            vec!["       0       3       1       3       6       2       6       9       1"]
        );
        let k = section_lines(&text, "BOND_FORCE_CONSTANT");
        assert_eq!(k.len(), 1);
        assert_eq!(k[0].split_whitespace().count(), 2);
        let restored = read_str(&text).unwrap();
        assert_eq!(restored.bonds, file.bonds);
    }

    #[test]
    fn long_sections_wrap_and_round_trip() {
        let mut file = PrmtopFile::new();
        for i in 0..25 {
            file.atoms
                .push(PrmtopAtom::new(format!("C{i}"), f64::from(i), 12.0, "CT"));
        }
        for i in 0..24usize {
            file.bonds
                .push(PrmtopBond::new(i, i + 1, 300.0 + i as f64, 1.5));
        }
        for i in 0..23usize {
            file.angles
                .push(PrmtopAngle::new(i + 2, i + 1, i, 50.0 + i as f64, 1.9));
        }
        let text = write_to_string(&file);
        assert_eq!(section_lines(&text, "ATOM_NAME").len(), 2);
        assert_eq!(section_lines(&text, "CHARGE").len(), 5);
        assert_eq!(section_lines(&text, "BONDS_WITHOUT_HYDROGEN").len(), 8);
        assert_eq!(section_lines(&text, "BOND_FORCE_CONSTANT").len(), 5);
        assert_eq!(section_lines(&text, "ANGLES_WITHOUT_HYDROGEN").len(), 10);
        let restored = read_str(&text).unwrap();
        assert_eq!(restored.atoms, file.atoms);
        assert_eq!(restored.bonds, file.bonds);
        assert_eq!(restored.angles, file.angles);
    }

    #[test]
    fn blank_atom_types_round_trip() {
        let mut file = PrmtopFile::new();
        for i in 0..25 {
            file.atoms
                .push(PrmtopAtom::new(format!("A{i}"), 0.0, 1.0, ""));
        }
        let restored = read_str(&write_to_string(&file)).unwrap();
        assert_eq!(restored.atoms, file.atoms);
    }

    /// Raw (unparsed) integers of section `flag` in a written topology.
    fn written_ints(text: &str, flag: &str) -> Vec<usize> {
        section_lines(text, flag)
            .iter()
            .flat_map(|line| line.split_whitespace().map(|t| t.parse::<usize>().unwrap()))
            .collect()
    }

    /// A tleap-style topology with preserved Lennard-Jones, dihedral and
    /// residue sections whose sizes are announced in POINTERS.
    fn tleap_style_text() -> String {
        // NATOM=3, NTYPES=3 (LJ types), NBONH=2, NPHIH=1, NRES=1, IFBOX=1
        let mut pointers = vec![0i64; POINTER_COUNT];
        pointers[PTR_NATOM] = 3;
        pointers[PTR_NTYPES] = 3;
        pointers[PTR_NBONH] = 2;
        pointers[6] = 1;
        pointers[11] = 1;
        pointers[27] = 1;
        let mut text = String::from("%VERSION  VERSION_STAMP = V0001.000\n");
        text.push_str(&flag("POINTERS", "10I8", &int_lines(&pointers)));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["O", "H1", "H2"])));
        text.push_str(&flag("CHARGE", "5E16.8", &real_lines(&[-1.0, 0.5, 0.5])));
        text.push_str(&flag("MASS", "5E16.8", &real_lines(&[16.0, 1.008, 1.008])));
        text.push_str(&flag(
            "AMBER_ATOM_TYPE",
            "20a4",
            &text_lines(&["OW", "HW", "HW"]),
        ));
        text.push_str(&flag(
            "BONDS_INC_HYDROGEN",
            "10I8",
            &int_lines(&[0, 3, 1, 0, 6, 1]),
        ));
        text.push_str(&flag(
            "BOND_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[553.0]),
        ));
        text.push_str(&flag("BOND_EQUIL_VALUE", "5E16.8", &real_lines(&[0.9572])));
        text.push_str(&flag("RESIDUE_LABEL", "20a4", &text_lines(&["WAT"])));
        text.push_str(&flag(
            "DIHEDRALS_INC_HYDROGEN",
            "10I8",
            &int_lines(&[3, 0, 6, 3, 1]),
        ));
        text.push_str(&flag(
            "NONBONDED_PARM_INDEX",
            "10I8",
            &int_lines(&[1, 2, 4, 2, 3, 5, 4, 5, 6]),
        ));
        text.push_str(&flag(
            "LENNARD_JONES_ACOEF",
            "5E16.8",
            &real_lines(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
        ));
        text.push_str(&flag(
            "LENNARD_JONES_BCOEF",
            "5E16.8",
            &real_lines(&[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]),
        ));
        text
    }

    #[test]
    fn pointers_of_preserved_sections_survive_a_round_trip() {
        let original = tleap_style_text();
        let file = read_str(&original).unwrap();
        assert_eq!(file.pointers.len(), POINTER_COUNT);
        assert_eq!(file.pointers[PTR_NTYPES], 3);
        let written = write_to_string(&file);

        let pointers = written_ints(&written, "POINTERS");
        assert_eq!(pointers.len(), POINTER_COUNT);
        assert_eq!(pointers[PTR_NATOM], 3);
        // counters of sections that are only preserved
        assert_eq!(
            pointers[PTR_NTYPES], 3,
            "NTYPES must stay the LJ type count"
        );
        assert_eq!(pointers[6], 1, "NPHIH");
        assert_eq!(pointers[11], 1, "NRES");
        assert_eq!(pointers[27], 1, "IFBOX");
        // counters of the sections the writer produces itself: the bonds are
        // now listed without hydrogen
        assert_eq!(pointers[PTR_NBONH], 0);
        assert_eq!(pointers[PTR_MBONA], 2);
        assert_eq!(pointers[PTR_NBONA], 2);
        assert_eq!(pointers[PTR_NUMBND], 1);

        // The preserved NONBONDED_PARM_INDEX is sized by NTYPES^2, so the
        // Lennard-Jones pairs stay resolvable after the round trip.
        let pairs = |text: &str| -> Vec<(usize, usize, f64, f64)> {
            crate::amber::extract_lj_parameters(text)
                .iter()
                .map(|p| (p.type_a, p.type_b, p.a_coeff, p.b_coeff))
                .collect()
        };
        let before = pairs(&original);
        assert_eq!(before.len(), 6);
        assert_eq!(pairs(&written), before);

        // And the file reads back to the same topology.
        let again = read_str(&written).unwrap();
        assert_eq!(again.pointers, pointers);
        assert_eq!(again.atoms, file.atoms);
        assert_eq!(again.bonds, file.bonds);
        assert_eq!(again.extra_sections, file.extra_sections);
    }

    #[test]
    fn pointers_default_to_the_distinct_type_count_for_in_memory_topologies() {
        let file = water();
        assert!(file.pointers.is_empty());
        let pointers = written_ints(&write_to_string(&file), "POINTERS");
        assert_eq!(pointers[PTR_NTYPES], 2);
        assert_eq!(pointers[PTR_NATYP], 2);
        assert_eq!(pointers[11], 0);
    }

    #[test]
    fn pointers_longer_than_the_standard_table_are_kept() {
        let mut file = water();
        file.pointers = vec![0; POINTER_COUNT + 1];
        file.pointers[POINTER_COUNT] = 7;
        let pointers = written_ints(&write_to_string(&file), "POINTERS");
        assert_eq!(pointers.len(), POINTER_COUNT + 1);
        assert_eq!(pointers[POINTER_COUNT], 7);
        assert_eq!(pointers[PTR_NATOM], 3);
    }

    #[test]
    fn empty_topology_round_trips() {
        let text = write_to_string(&PrmtopFile::new());
        let restored = read_str(&text).unwrap();
        assert_eq!(restored.natom(), 0);
        assert!(restored.bonds.is_empty());
        assert!(restored.angles.is_empty());
        assert!(restored.extra_sections.is_empty());
        assert!(section_lines(&text, "POINTERS")[0].starts_with("       0       0"));
    }

    #[test]
    fn writer_rejects_terms_that_reference_missing_atoms() {
        let mut file = water();
        file.bonds.push(PrmtopBond::new(0, 3, 553.0, 0.9572));
        let mut buf = Vec::new();
        let err = file.write(&mut buf).unwrap_err();
        assert!(matches!(err, AmberIoError::InvalidData(_)), "{err}");
        assert!(err.to_string().contains("atom 3"), "{err}");
        assert!(buf.is_empty(), "nothing may be written on failure");

        let mut file = water();
        file.angles.push(PrmtopAngle::new(0, 1, 7, 1.0, 1.0));
        let err = file.write(&mut Vec::new()).unwrap_err();
        assert!(matches!(err, AmberIoError::InvalidData(_)), "{err}");
    }

    #[test]
    fn writer_rejects_text_that_does_not_fit_its_column() {
        let mut file = water();
        file.atoms[1].name = "HYDRO".to_string();
        assert!(matches!(
            file.write(&mut Vec::new()),
            Err(AmberIoError::InvalidData(_))
        ));
        let mut file = water();
        file.atoms[2].atom_type = "H\nW".to_string();
        assert!(matches!(
            file.write(&mut Vec::new()),
            Err(AmberIoError::InvalidData(_))
        ));
        let mut file = water();
        file.version = "two\nlines".to_string();
        assert!(matches!(
            file.write(&mut Vec::new()),
            Err(AmberIoError::InvalidData(_))
        ));
    }

    #[test]
    fn writer_rejects_names_that_would_be_read_as_directives() {
        // `%FLA` + `G   ` would put "%FLAG" at the start of the record.
        let mut file = water();
        file.atoms[0].name = "%FLA".to_string();
        file.atoms[1].name = "G".to_string();
        let mut buf = Vec::new();
        let err = file.write(&mut buf).unwrap_err();
        assert!(matches!(err, AmberIoError::InvalidData(_)), "{err}");
        assert!(err.to_string().contains("no `%`"), "{err}");
        assert!(buf.is_empty());

        let mut file = water();
        file.atoms[2].atom_type = "H%".to_string();
        assert!(matches!(
            file.write(&mut Vec::new()),
            Err(AmberIoError::InvalidData(_))
        ));

        // Names with other punctuation are untouched.
        let mut file = water();
        file.atoms[0].name = "O'".to_string();
        file.atoms[1].name = "H*1".to_string();
        let restored = read_str(&write_to_string(&file)).unwrap();
        assert_eq!(restored.atoms, file.atoms);
    }

    #[test]
    fn writer_rejects_integers_wider_than_their_column() {
        let mut out = Vec::new();
        let err = write_int_section(
            &mut out,
            "BONDS_WITHOUT_HYDROGEN",
            &[99_999_999, 100_000_000],
        )
        .unwrap_err();
        assert!(matches!(err, AmberIoError::InvalidData(_)), "{err}");
        assert!(err.to_string().contains("BONDS_WITHOUT_HYDROGEN"), "{err}");
        let mut out = Vec::new();
        write_int_section(&mut out, "X", &[99_999_999, 1]).unwrap();
        assert!(
            String::from_utf8(out)
                .unwrap()
                .contains("99999999       1\n")
        );
    }

    #[test]
    fn full_width_reals_round_trip() {
        // "-1.00000000E-300" fills its E16.8 column and touches its neighbours.
        let mut file = water();
        file.atoms[0].charge = 1.0;
        file.atoms[1].charge = -1.0e-300;
        file.atoms[2].charge = -2.0e-300;
        let text = write_to_string(&file);
        let charge_line = &section_lines(&text, "CHARGE")[0];
        assert!(charge_line.split_whitespace().any(|t| t.len() > REAL_WIDTH));
        let restored = read_str(&text).unwrap();
        assert_eq!(restored.atoms, file.atoms);
    }

    #[test]
    fn writer_rejects_unsafe_extra_sections() {
        let mut file = water();
        file.extra_sections
            .insert("POINTERS".to_string(), vec!["%FORMAT(10I8)".to_string()]);
        assert!(matches!(
            file.write(&mut Vec::new()),
            Err(AmberIoError::InvalidData(_))
        ));
        let mut file = water();
        file.extra_sections
            .insert("BAD NAME".to_string(), Vec::new());
        assert!(matches!(
            file.write(&mut Vec::new()),
            Err(AmberIoError::InvalidData(_))
        ));
        let mut file = water();
        file.extra_sections
            .insert("NOTES".to_string(), vec!["%FLAG HIJACK".to_string()]);
        assert!(matches!(
            file.write(&mut Vec::new()),
            Err(AmberIoError::InvalidData(_))
        ));
    }

    /// A sink that rejects every write.
    struct ClosedSink;

    impl Write for ClosedSink {
        fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("sink closed"))
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("sink closed"))
        }
    }

    #[test]
    fn writer_reports_a_failing_sink() {
        let err = water().write(&mut ClosedSink).unwrap_err();
        assert!(matches!(err, AmberIoError::IoError(_)), "{err}");
        assert!(err.to_string().contains("sink closed"), "{err}");
    }

    // ── reading ───────────────────────────────────────────────────────────────

    #[test]
    fn read_decodes_hydrogen_sections_before_heavy_atom_sections() {
        let file = read_str(&four_atom_text()).unwrap();
        assert_eq!(file.version, "VERSION_STAMP = V0001.000");
        let names: Vec<&str> = file.atoms.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, vec!["C1", "H1", "C2", "H2"]);
        assert_eq!(file.atoms[2].atom_type, "CT");
        assert_eq!(file.atoms[3].mass, 1.008);
        assert_eq!(
            file.bonds,
            vec![
                PrmtopBond::new(0, 1, 340.0, 1.09),
                PrmtopBond::new(2, 3, 310.0, 1.526),
            ]
        );
        assert_eq!(
            file.angles,
            vec![
                PrmtopAngle::new(1, 0, 2, 63.0, 1.94),
                PrmtopAngle::new(0, 2, 3, 50.0, 1.91),
            ]
        );
    }

    #[test]
    fn read_separates_touching_integer_fields() {
        // NTYPES = 12345678 fills its I8 column and touches NATOM = 3.
        let pointers = int_lines(&[3, 12_345_678, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert!(pointers.contains("312345678"));
        let mut text = flag("POINTERS", "10I8", &pointers);
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B", "C"])));
        let file = read_str(&text).unwrap();
        assert_eq!(file.natom(), 3);
    }

    #[test]
    fn read_accepts_names_whose_trailing_padding_was_stripped() {
        let mut text = flag("POINTERS", "10I8", &int_lines(&[3]));
        text.push_str("%FLAG ATOM_NAME\n%FORMAT(20a4)\nN   CA  C\n");
        let file = read_str(&text).unwrap();
        let names: Vec<&str> = file.atoms.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, vec!["N", "CA", "C"]);
    }

    #[test]
    fn read_defaults_absent_per_atom_sections() {
        let mut text = flag("POINTERS", "10I8", &int_lines(&[2]));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"])));
        let file = read_str(&text).unwrap();
        assert_eq!(file.atoms.len(), 2);
        assert_eq!(file.atoms[1].charge, 0.0);
        assert_eq!(file.atoms[1].mass, 0.0);
        assert_eq!(file.atoms[1].atom_type, "");
    }

    #[test]
    fn read_without_pointers_counts_the_atom_names() {
        let text = flag("ATOM_NAME", "20a4", &text_lines(&["A", "B", "C"]));
        assert_eq!(read_str(&text).unwrap().natom(), 3);
    }

    #[test]
    fn read_rejects_bond_atoms_beyond_natom() {
        let mut text = flag("POINTERS", "10I8", &int_lines(&[2]));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"])));
        text.push_str(&flag(
            "BONDS_WITHOUT_HYDROGEN",
            "10I8",
            &int_lines(&[0, 6, 1]),
        ));
        text.push_str(&flag(
            "BOND_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[300.0]),
        ));
        text.push_str(&flag("BOND_EQUIL_VALUE", "5E16.8", &real_lines(&[1.5])));
        let err = read_str(&text).unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");
        assert!(err.to_string().contains("atom pointer 6"), "{err}");
        assert!(err.to_string().contains("BONDS_WITHOUT_HYDROGEN"), "{err}");
    }

    #[test]
    fn read_rejects_invalid_atom_pointers() {
        for bad in [4_i64, -3, 5] {
            let mut text = flag("POINTERS", "10I8", &int_lines(&[3]));
            text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B", "C"])));
            text.push_str(&flag(
                "ANGLES_INC_HYDROGEN",
                "10I8",
                &int_lines(&[0, 3, bad, 1]),
            ));
            text.push_str(&flag(
                "ANGLE_FORCE_CONSTANT",
                "5E16.8",
                &real_lines(&[50.0]),
            ));
            text.push_str(&flag("ANGLE_EQUIL_VALUE", "5E16.8", &real_lines(&[1.9])));
            let err = read_str(&text).unwrap_err();
            assert!(
                matches!(err, AmberIoError::ParseError(_)),
                "pointer {bad}: {err}"
            );
        }
    }

    #[test]
    fn read_rejects_parameter_types_outside_the_tables() {
        for bad in [0_i64, 2, -1] {
            let mut text = flag("POINTERS", "10I8", &int_lines(&[2]));
            text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"])));
            text.push_str(&flag(
                "BONDS_WITHOUT_HYDROGEN",
                "10I8",
                &int_lines(&[0, 3, bad]),
            ));
            text.push_str(&flag(
                "BOND_FORCE_CONSTANT",
                "5E16.8",
                &real_lines(&[300.0]),
            ));
            text.push_str(&flag("BOND_EQUIL_VALUE", "5E16.8", &real_lines(&[1.5])));
            let err = read_str(&text).unwrap_err();
            assert!(
                matches!(err, AmberIoError::ParseError(_)),
                "type {bad}: {err}"
            );
            assert!(err.to_string().contains("parameter type"), "{err}");
        }
    }

    #[test]
    fn read_rejects_a_trailing_partial_record() {
        let mut text = flag("POINTERS", "10I8", &int_lines(&[2]));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"])));
        text.push_str(&flag(
            "BONDS_WITHOUT_HYDROGEN",
            "10I8",
            &int_lines(&[0, 3, 1, 3]),
        ));
        text.push_str(&flag(
            "BOND_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[300.0]),
        ));
        text.push_str(&flag("BOND_EQUIL_VALUE", "5E16.8", &real_lines(&[1.5])));
        let err = read_str(&text).unwrap_err();
        assert!(err.to_string().contains("whole records"), "{err}");
    }

    #[test]
    fn read_rejects_bond_tables_without_index_sections() {
        // The layout written by earlier versions: parameter tables only.
        let mut text = flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"]));
        text.push_str(&flag(
            "BOND_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[300.0]),
        ));
        text.push_str(&flag("BOND_EQUIL_VALUE", "5E16.8", &real_lines(&[1.5])));
        let err = read_str(&text).unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");
        assert!(err.to_string().contains("BONDS_WITHOUT_HYDROGEN"), "{err}");

        let mut text = flag("ATOM_NAME", "20a4", &text_lines(&["A", "B", "C"]));
        text.push_str(&flag(
            "ANGLE_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[50.0]),
        ));
        text.push_str(&flag("ANGLE_EQUIL_VALUE", "5E16.8", &real_lines(&[1.9])));
        let err = read_str(&text).unwrap_err();
        assert!(err.to_string().contains("ANGLES_WITHOUT_HYDROGEN"), "{err}");
    }

    #[test]
    fn read_rejects_mismatched_parameter_tables() {
        let mut text = flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"]));
        text.push_str(&flag(
            "BONDS_WITHOUT_HYDROGEN",
            "10I8",
            &int_lines(&[0, 3, 1]),
        ));
        text.push_str(&flag(
            "BOND_FORCE_CONSTANT",
            "5E16.8",
            &real_lines(&[300.0, 310.0]),
        ));
        text.push_str(&flag("BOND_EQUIL_VALUE", "5E16.8", &real_lines(&[1.5])));
        let err = read_str(&text).unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");
    }

    #[test]
    fn read_rejects_inconsistent_atom_counts() {
        let mut text = flag("POINTERS", "10I8", &int_lines(&[3]));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"])));
        let err = read_str(&text).unwrap_err();
        assert!(matches!(err, AmberIoError::ParseError(_)), "{err}");

        let mut text = flag("POINTERS", "10I8", &int_lines(&[2]));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"])));
        text.push_str(&flag("CHARGE", "5E16.8", &real_lines(&[1.0])));
        let err = read_str(&text).unwrap_err();
        assert!(err.to_string().contains("CHARGE"), "{err}");

        let text = flag("POINTERS", "10I8", &int_lines(&[2]));
        assert!(matches!(
            read_str(&text),
            Err(AmberIoError::MissingSection(_))
        ));

        let text = flag("POINTERS", "10I8", &int_lines(&[-1]));
        assert!(matches!(read_str(&text), Err(AmberIoError::ParseError(_))));
    }

    #[test]
    fn read_rejects_duplicate_and_malformed_sections() {
        let mut text = flag("ATOM_NAME", "20a4", &text_lines(&["A"]));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["B"])));
        let err = read_str(&text).unwrap_err();
        assert!(err.to_string().contains("duplicate"), "{err}");

        let mut text = flag("ATOM_NAME", "20a4", &text_lines(&["A"]));
        text.push_str("%FLAG CHARGE\n%FORMAT(5E16.8)\n not-a-number\n");
        let err = read_str(&text).unwrap_err();
        assert!(err.to_string().contains("CHARGE"), "{err}");

        assert!(matches!(
            read_str("%FLAG\n"),
            Err(AmberIoError::ParseError(_))
        ));
    }

    #[test]
    fn unknown_sections_are_preserved_verbatim() {
        let mut text = String::from("%VERSION  V\n");
        text.push_str(&flag("TITLE", "20a4", "water box"));
        text.push_str(&flag("POINTERS", "10I8", &int_lines(&[2])));
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B"])));
        text.push_str("%FLAG RESIDUE_LABEL\n%COMMENT residue names\n%FORMAT(20a4)\nWAT \n");
        text.push_str(&flag("ATOM_TYPE_INDEX", "10I8", &int_lines(&[1, 2])));
        let file = read_str(&text).unwrap();
        assert_eq!(
            file.extra_sections["TITLE"],
            vec!["%FORMAT(20a4)".to_string(), "water box".to_string()]
        );
        assert_eq!(
            file.extra_sections["RESIDUE_LABEL"],
            vec![
                "%COMMENT residue names".to_string(),
                "%FORMAT(20a4)".to_string(),
                "WAT ".to_string()
            ]
        );
        assert_eq!(file.extra_sections["ATOM_TYPE_INDEX"].len(), 2);
        assert_eq!(file.extra_sections.len(), 3);
        assert!(!file.extra_sections.contains_key("POINTERS"));
        assert!(!file.extra_sections.contains_key("ATOM_NAME"));

        let again = read_str(&write_to_string(&file)).unwrap();
        assert_eq!(again.extra_sections, file.extra_sections);
        assert_eq!(again.atoms, file.atoms);
    }

    #[test]
    fn extra_sections_are_written_sorted_by_name() {
        let mut file = PrmtopFile::new();
        file.extra_sections
            .insert("ZETA".to_string(), vec!["z".to_string()]);
        file.extra_sections
            .insert("ALPHA".to_string(), vec!["a".to_string()]);
        let text = write_to_string(&file);
        let alpha = text.find("%FLAG ALPHA").unwrap();
        let zeta = text.find("%FLAG ZETA").unwrap();
        assert!(alpha < zeta);
    }

    // ── helpers ───────────────────────────────────────────────────────────────

    #[test]
    fn format_width_reads_the_declared_column_width() {
        assert_eq!(format_width("%FORMAT(10I8)"), Some(8));
        assert_eq!(format_width("%FORMAT(5E16.8)"), Some(16));
        assert_eq!(format_width("  %FORMAT(20a4)  "), Some(4));
        assert_eq!(format_width("%FORMAT(1a80)"), Some(80));
        assert_eq!(format_width("%FORMAT()"), None);
        assert_eq!(format_width("%FORMAT(I)"), None);
        assert_eq!(format_width("10I8"), None);
    }

    #[test]
    fn decode_pointer_accepts_only_multiples_of_three_inside_the_topology() {
        assert_eq!(decode_pointer(0, 3), Some(0));
        assert_eq!(decode_pointer(3, 3), Some(1));
        assert_eq!(decode_pointer(6, 3), Some(2));
        assert_eq!(decode_pointer(9, 3), None);
        assert_eq!(decode_pointer(4, 3), None);
        assert_eq!(decode_pointer(-3, 3), None);
    }

    #[test]
    fn declared_format_width_overrides_the_default() {
        // A 6-wide integer format: the 7-digit run is two touching fields.
        let mut text = String::from("%FLAG POINTERS\n%FORMAT(10I6)\n     3123456\n");
        text.push_str(&flag("ATOM_NAME", "20a4", &text_lines(&["A", "B", "C"])));
        let file = read_str(&text).unwrap();
        assert_eq!(file.natom(), 3);
    }

    #[test]
    fn mutated_topologies_never_panic_and_stay_consistent() {
        use crate::fortran_fields::test_support::{SplitMix64, mutate_text};
        let seeds = [
            write_to_string(&water()),
            four_atom_text(),
            tleap_style_text(),
        ];
        let mut rng = SplitMix64::new(0x9A27_0001);
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        for seed in &seeds {
            for round in 0..1500 {
                let text = mutate_text(seed, &mut rng, 1 + round % 4);
                match read_str(&text) {
                    Ok(file) => {
                        accepted += 1;
                        let natom = file.natom();
                        if let Some(&declared) = file.pointers.first() {
                            assert_eq!(declared, natom, "NATOM disagrees with the atoms:\n{text}");
                        }
                        for bond in &file.bonds {
                            assert!(
                                bond.atom_i < natom && bond.atom_j < natom,
                                "bond outside the topology:\n{text}"
                            );
                        }
                        for angle in &file.angles {
                            assert!(
                                angle.atom_i < natom
                                    && angle.atom_j < natom
                                    && angle.atom_k < natom,
                                "angle outside the topology:\n{text}"
                            );
                        }
                    }
                    Err(_) => rejected += 1,
                }
            }
        }
        assert!(accepted > 0, "no damaged topology was accepted");
        assert!(rejected > 0, "no damaged topology was rejected");
    }
}
