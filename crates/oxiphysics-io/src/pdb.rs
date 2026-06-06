// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PDB (Protein Data Bank) format writer and reader.
//!
//! Supports ATOM, HETATM, CONECT, CRYST1, multi-MODEL PDB files,
//! occupancy/B-factor handling, and chain ID operations.

use crate::{Error, Result};
use oxiphysics_core::math::Vec3;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Represents a single atom record in PDB format.
#[derive(Debug, Clone)]
pub struct PdbAtom {
    /// Atom serial number.
    pub serial: u32,
    /// Atom name (e.g., "CA", "N").
    pub name: String,
    /// Residue name (e.g., "ALA", "GLY").
    pub res_name: String,
    /// Chain identifier (single character, e.g., 'A').
    pub chain_id: char,
    /// Residue sequence number.
    pub res_seq: u32,
    /// X coordinate in angstroms.
    pub x: f64,
    /// Y coordinate in angstroms.
    pub y: f64,
    /// Z coordinate in angstroms.
    pub z: f64,
    /// Occupancy factor.
    pub occupancy: f64,
    /// Temperature (B) factor.
    pub temp_factor: f64,
    /// Element symbol.
    pub element: String,
    /// Whether this is a HETATM record (vs ATOM).
    pub is_hetatm: bool,
}

/// All named fields for constructing a [`PdbAtom`] via [`PdbAtom::new`].
///
/// Groups the 11 record fields so the constructor stays within the
/// argument-count limit.
#[derive(Debug, Clone)]
pub struct PdbAtomData {
    /// Atom serial number.
    pub serial: u32,
    /// Atom name (e.g., `"CA"`, `"N"`).
    pub name: String,
    /// Residue name (e.g., `"ALA"`, `"GLY"`).
    pub res_name: String,
    /// Chain identifier.
    pub chain_id: char,
    /// Residue sequence number.
    pub res_seq: u32,
    /// X coordinate \[Å\]
    pub x: f64,
    /// Y coordinate \[Å\]
    pub y: f64,
    /// Z coordinate \[Å\]
    pub z: f64,
    /// Occupancy factor.
    pub occupancy: f64,
    /// Temperature (B) factor.
    pub temp_factor: f64,
    /// Element symbol.
    pub element: String,
}

impl PdbAtom {
    /// Create a `PdbAtom` from a [`PdbAtomData`] record.
    pub fn new(data: PdbAtomData) -> Self {
        PdbAtom {
            serial: data.serial,
            name: data.name,
            res_name: data.res_name,
            chain_id: data.chain_id,
            res_seq: data.res_seq,
            x: data.x,
            y: data.y,
            z: data.z,
            occupancy: data.occupancy,
            temp_factor: data.temp_factor,
            element: data.element,
            is_hetatm: false,
        }
    }

    /// Return position as a `Vec3`.
    pub fn position(&self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }

    /// Compute distance to another atom.
    pub fn distance_to(&self, other: &PdbAtom) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// A CONECT record: bonds between atoms.
#[derive(Debug, Clone)]
pub struct ConectRecord {
    /// Serial number of the atom.
    pub atom_serial: u32,
    /// Serial numbers of bonded atoms.
    pub bonded_atoms: Vec<u32>,
}

/// A CRYST1 record: unit cell parameters.
#[derive(Debug, Clone)]
pub struct Cryst1Record {
    /// Unit cell lengths a, b, c in angstroms.
    pub a: f64,
    /// Unit cell length b.
    pub b: f64,
    /// Unit cell length c.
    pub c: f64,
    /// Unit cell angle alpha in degrees.
    pub alpha: f64,
    /// Unit cell angle beta in degrees.
    pub beta: f64,
    /// Unit cell angle gamma in degrees.
    pub gamma: f64,
    /// Space group.
    pub space_group: String,
    /// Z value (number of molecules per unit cell).
    pub z_value: u32,
}

impl Default for Cryst1Record {
    fn default() -> Self {
        Cryst1Record {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
            space_group: "P 1".to_string(),
            z_value: 1,
        }
    }
}

/// A PDB model: a collection of atoms with optional remarks, title, and metadata.
#[derive(Debug, Clone)]
pub struct PdbModel {
    /// All ATOM/HETATM records.
    pub atoms: Vec<PdbAtom>,
    /// REMARK lines (text after the "REMARK" keyword).
    pub remarks: Vec<String>,
    /// TITLE line content.
    pub title: String,
    /// CONECT records.
    pub conect_records: Vec<ConectRecord>,
    /// CRYST1 record (unit cell).
    pub cryst1: Option<Cryst1Record>,
    /// Model number (for multi-model PDB files).
    pub model_number: Option<u32>,
}

/// A multi-model PDB file.
#[derive(Debug, Clone)]
pub struct PdbFile {
    /// All models in the PDB file.
    pub models: Vec<PdbModel>,
    /// File-level remarks.
    pub remarks: Vec<String>,
    /// File-level title.
    pub title: String,
    /// File-level CRYST1 record.
    pub cryst1: Option<Cryst1Record>,
}

/// Parse a single ATOM or HETATM line into a [`PdbAtom`].
///
/// Returns `None` if the line is not an ATOM/HETATM record or is too short to parse.
pub fn parse_atom_line(line: &str) -> Option<PdbAtom> {
    let record = &line[..line.len().min(6)];
    let is_hetatm = record == "HETATM";
    if record != "ATOM  " && !is_hetatm {
        return None;
    }
    if line.len() < 54 {
        return None;
    }

    let serial = line.get(6..11)?.trim().parse::<u32>().ok()?;
    let name = line.get(12..16)?.trim().to_string();
    let res_name = line.get(17..20)?.trim().to_string();
    let chain_id = line.as_bytes().get(21).map_or('A', |&b| b as char);
    let res_seq = line
        .get(22..26)
        .and_then(|s| s.trim().parse::<u32>().ok())
        .unwrap_or(1);
    let x = line.get(30..38)?.trim().parse::<f64>().ok()?;
    let y = line.get(38..46)?.trim().parse::<f64>().ok()?;
    let z = line.get(46..54)?.trim().parse::<f64>().ok()?;
    let occupancy = line
        .get(54..60)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(1.0);
    let temp_factor = line
        .get(60..66)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(0.0);
    let element = line
        .get(76..78)
        .map(|s| s.trim().to_string())
        .unwrap_or_default();

    let mut atom = PdbAtom::new(PdbAtomData {
        serial,
        name,
        res_name,
        chain_id,
        res_seq,
        x,
        y,
        z,
        occupancy,
        temp_factor,
        element,
    });
    atom.is_hetatm = is_hetatm;
    Some(atom)
}

/// Parse a CONECT record line.
pub fn parse_conect_line(line: &str) -> Option<ConectRecord> {
    if !line.starts_with("CONECT") {
        return None;
    }
    let parts: Vec<u32> = line[6..]
        .split_whitespace()
        .filter_map(|s| s.parse::<u32>().ok())
        .collect();
    if parts.is_empty() {
        return None;
    }
    Some(ConectRecord {
        atom_serial: parts[0],
        bonded_atoms: parts[1..].to_vec(),
    })
}

/// Parse a CRYST1 record line.
pub fn parse_cryst1_line(line: &str) -> Option<Cryst1Record> {
    if !line.starts_with("CRYST1") {
        return None;
    }
    if line.len() < 54 {
        return None;
    }
    let a = line.get(6..15)?.trim().parse::<f64>().ok()?;
    let b = line.get(15..24)?.trim().parse::<f64>().ok()?;
    let c = line.get(24..33)?.trim().parse::<f64>().ok()?;
    let alpha = line.get(33..40)?.trim().parse::<f64>().ok()?;
    let beta = line.get(40..47)?.trim().parse::<f64>().ok()?;
    let gamma = line.get(47..54)?.trim().parse::<f64>().ok()?;
    let space_group = line.get(55..66).unwrap_or("P 1").trim().to_string();
    let z_value = line
        .get(66..70)
        .and_then(|s| s.trim().parse::<u32>().ok())
        .unwrap_or(1);

    Some(Cryst1Record {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        space_group,
        z_value,
    })
}

/// Parse PDB file content (as a `&str`) into a [`PdbModel`].
///
/// Handles TITLE, REMARK, ATOM, HETATM, CONECT, and CRYST1 records.
pub fn parse_pdb(content: &str) -> Result<PdbModel> {
    let mut atoms = Vec::new();
    let mut remarks = Vec::new();
    let mut title = String::new();
    let mut conect_records = Vec::new();
    let mut cryst1 = None;

    for line in content.lines() {
        if line.starts_with("TITLE") {
            title = line.get(10..).unwrap_or("").trim().to_string();
        } else if line.starts_with("REMARK") {
            remarks.push(line.get(6..).unwrap_or("").trim().to_string());
        } else if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            match parse_atom_line(line) {
                Some(atom) => atoms.push(atom),
                None => return Err(Error::Parse(format!("Failed to parse atom line: {line}"))),
            }
        } else if line.starts_with("CONECT") {
            if let Some(conect) = parse_conect_line(line) {
                conect_records.push(conect);
            }
        } else if line.starts_with("CRYST1") {
            cryst1 = parse_cryst1_line(line);
        }
    }

    Ok(PdbModel {
        atoms,
        remarks,
        title,
        conect_records,
        cryst1,
        model_number: None,
    })
}

/// Parse a multi-model PDB file.
pub fn parse_multi_model_pdb(content: &str) -> Result<PdbFile> {
    let mut models = Vec::new();
    let mut file_remarks = Vec::new();
    let mut file_title = String::new();
    let mut file_cryst1 = None;

    let mut current_atoms: Vec<PdbAtom> = Vec::new();
    let mut current_model_num: Option<u32> = None;
    let mut current_conect: Vec<ConectRecord> = Vec::new();
    let mut in_model = false;

    for line in content.lines() {
        if line.starts_with("TITLE") {
            file_title = line.get(10..).unwrap_or("").trim().to_string();
        } else if line.starts_with("REMARK") {
            file_remarks.push(line.get(6..).unwrap_or("").trim().to_string());
        } else if line.starts_with("CRYST1") {
            file_cryst1 = parse_cryst1_line(line);
        } else if let Some(rest) = line.strip_prefix("MODEL") {
            in_model = true;
            current_model_num = rest.trim().parse::<u32>().ok();
            current_atoms.clear();
            current_conect.clear();
        } else if line.starts_with("ENDMDL") {
            models.push(PdbModel {
                atoms: std::mem::take(&mut current_atoms),
                remarks: Vec::new(),
                title: String::new(),
                conect_records: std::mem::take(&mut current_conect),
                cryst1: file_cryst1.clone(),
                model_number: current_model_num,
            });
            in_model = false;
        } else if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            if let Some(atom) = parse_atom_line(line) {
                current_atoms.push(atom);
            }
        } else if line.starts_with("CONECT")
            && let Some(conect) = parse_conect_line(line)
        {
            current_conect.push(conect);
        }
    }

    // If no MODEL/ENDMDL was found, treat whole file as one model
    if models.is_empty() && !current_atoms.is_empty() {
        models.push(PdbModel {
            atoms: current_atoms,
            remarks: file_remarks.clone(),
            title: file_title.clone(),
            conect_records: current_conect,
            cryst1: file_cryst1.clone(),
            model_number: None,
        });
    } else if !in_model && !current_atoms.is_empty() {
        // Atoms found after last ENDMDL
        models.push(PdbModel {
            atoms: current_atoms,
            remarks: Vec::new(),
            title: String::new(),
            conect_records: current_conect,
            cryst1: file_cryst1.clone(),
            model_number: None,
        });
    }

    Ok(PdbFile {
        models,
        remarks: file_remarks,
        title: file_title,
        cryst1: file_cryst1,
    })
}

/// Write a [`PdbModel`] to a PDB fixed-column format string.
pub fn write_pdb(model: &PdbModel) -> String {
    let mut out = String::new();

    if !model.title.is_empty() {
        out.push_str(&format!("TITLE     {}\n", model.title));
    }
    for remark in &model.remarks {
        out.push_str(&format!("REMARK {remark}\n"));
    }

    if let Some(ref cryst) = model.cryst1 {
        out.push_str(&format!(
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {:<11}{:>4}\n",
            cryst.a,
            cryst.b,
            cryst.c,
            cryst.alpha,
            cryst.beta,
            cryst.gamma,
            cryst.space_group,
            cryst.z_value,
        ));
    }

    for atom in &model.atoms {
        let record_type = if atom.is_hetatm { "HETATM" } else { "ATOM  " };
        let atom_name = if atom.name.len() < 4 {
            format!(" {:<3}", atom.name)
        } else {
            format!("{:<4}", atom.name)
        };
        let element_field = if atom.element.is_empty() {
            "  ".to_string()
        } else {
            format!("{:>2}", atom.element)
        };
        out.push_str(&format!(
            "{}{:>5} {} {:>3} {}{:>4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {}\n",
            record_type,
            atom.serial,
            atom_name,
            atom.res_name,
            atom.chain_id,
            atom.res_seq,
            atom.x,
            atom.y,
            atom.z,
            atom.occupancy,
            atom.temp_factor,
            element_field,
        ));
    }

    for conect in &model.conect_records {
        let bonded: Vec<String> = conect
            .bonded_atoms
            .iter()
            .map(|a| format!("{:>5}", a))
            .collect();
        out.push_str(&format!(
            "CONECT{:>5}{}\n",
            conect.atom_serial,
            bonded.join("")
        ));
    }

    out.push_str("END\n");
    out
}

/// Write a multi-model PDB file.
pub fn write_multi_model_pdb(pdb_file: &PdbFile) -> String {
    let mut out = String::new();

    if !pdb_file.title.is_empty() {
        out.push_str(&format!("TITLE     {}\n", pdb_file.title));
    }
    for remark in &pdb_file.remarks {
        out.push_str(&format!("REMARK {remark}\n"));
    }
    if let Some(ref cryst) = pdb_file.cryst1 {
        out.push_str(&format!(
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {:<11}{:>4}\n",
            cryst.a,
            cryst.b,
            cryst.c,
            cryst.alpha,
            cryst.beta,
            cryst.gamma,
            cryst.space_group,
            cryst.z_value,
        ));
    }

    for (i, model) in pdb_file.models.iter().enumerate() {
        let model_num = model.model_number.unwrap_or((i + 1) as u32);
        out.push_str(&format!("MODEL     {:>4}\n", model_num));

        for atom in &model.atoms {
            let record_type = if atom.is_hetatm { "HETATM" } else { "ATOM  " };
            let atom_name = if atom.name.len() < 4 {
                format!(" {:<3}", atom.name)
            } else {
                format!("{:<4}", atom.name)
            };
            let element_field = if atom.element.is_empty() {
                "  ".to_string()
            } else {
                format!("{:>2}", atom.element)
            };
            out.push_str(&format!(
                "{}{:>5} {} {:>3} {}{:>4}    {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {}\n",
                record_type,
                atom.serial,
                atom_name,
                atom.res_name,
                atom.chain_id,
                atom.res_seq,
                atom.x,
                atom.y,
                atom.z,
                atom.occupancy,
                atom.temp_factor,
                element_field,
            ));
        }
        for conect in &model.conect_records {
            let bonded: Vec<String> = conect
                .bonded_atoms
                .iter()
                .map(|a| format!("{:>5}", a))
                .collect();
            out.push_str(&format!(
                "CONECT{:>5}{}\n",
                conect.atom_serial,
                bonded.join("")
            ));
        }
        out.push_str("ENDMDL\n");
    }
    out.push_str("END\n");
    out
}

/// Filter atoms by chain ID.
pub fn filter_by_chain(atoms: &[PdbAtom], chain: char) -> Vec<PdbAtom> {
    atoms
        .iter()
        .filter(|a| a.chain_id == chain)
        .cloned()
        .collect()
}

/// Get all unique chain IDs from a list of atoms.
pub fn unique_chain_ids(atoms: &[PdbAtom]) -> Vec<char> {
    let mut chains: Vec<char> = atoms.iter().map(|a| a.chain_id).collect();
    chains.sort_unstable();
    chains.dedup();
    chains
}

/// Get all unique residue names from a list of atoms.
pub fn unique_residue_names(atoms: &[PdbAtom]) -> Vec<String> {
    let mut names: Vec<String> = atoms.iter().map(|a| a.res_name.clone()).collect();
    names.sort();
    names.dedup();
    names
}

/// Compute the center of mass of atoms (equal mass).
pub fn center_of_mass(atoms: &[PdbAtom]) -> [f64; 3] {
    if atoms.is_empty() {
        return [0.0; 3];
    }
    let n = atoms.len() as f64;
    let mut com = [0.0; 3];
    for a in atoms {
        com[0] += a.x;
        com[1] += a.y;
        com[2] += a.z;
    }
    [com[0] / n, com[1] / n, com[2] / n]
}

// ---------------------------------------------------------------------------
// Legacy file-based writer / reader (kept for backward compatibility)
// ---------------------------------------------------------------------------

/// Writer for PDB format files.
pub struct PdbWriter;

impl PdbWriter {
    /// Write a single frame of atom data to a PDB file.
    ///
    /// Each atom is written as an ATOM record following PDB column conventions.
    pub fn write_frame(path: &str, atoms: &[PdbAtom]) -> Result<()> {
        let model = PdbModel {
            atoms: atoms.to_vec(),
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let content = write_pdb(&model);
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        w.write_all(content.as_bytes())?;
        w.flush()?;
        Ok(())
    }
}

/// Reader for PDB format files (basic ATOM records only).
pub struct PdbReader;

impl PdbReader {
    /// Read ATOM records from a PDB file.
    pub fn read(path: &str) -> Result<Vec<PdbAtom>> {
        let file = File::open(Path::new(path))?;
        let reader = BufReader::new(file);
        let mut content = String::new();
        for line in reader.lines() {
            content.push_str(&line?);
            content.push('\n');
        }
        let model = parse_pdb(&content)?;
        Ok(model.atoms)
    }
}

// ---------------------------------------------------------------------------
// SEQRES Parsing
// ---------------------------------------------------------------------------

/// A SEQRES record describing the primary sequence of a chain.
#[derive(Debug, Clone)]
pub struct SeqresRecord {
    /// Chain identifier.
    pub chain_id: char,
    /// Total number of residues in the chain.
    pub num_residues: u32,
    /// List of 3-letter residue names.
    pub residues: Vec<String>,
}

impl SeqresRecord {
    /// Return the sequence as a dash-separated string of 3-letter codes.
    pub fn sequence_string(&self) -> String {
        self.residues.join("-")
    }

    /// Return a 1-letter code sequence (standard amino acids only; unknown → '?').
    pub fn one_letter_sequence(&self) -> String {
        self.residues.iter().map(|r| aa3_to_1(r)).collect()
    }
}

/// Parse all SEQRES records from PDB content.
pub fn parse_seqres(content: &str) -> Vec<SeqresRecord> {
    // Group SEQRES lines by chain ID
    let mut chains: std::collections::HashMap<char, (u32, Vec<String>)> =
        std::collections::HashMap::new();

    for line in content.lines() {
        if !line.starts_with("SEQRES") {
            continue;
        }
        if line.len() < 19 {
            continue;
        }

        let chain_id = line.as_bytes().get(11).map_or('A', |&b| b as char);
        let num_res: u32 = line
            .get(13..17)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0);

        // Residue names start at col 19, each is 4 chars (padded with space)
        let residues_part = line.get(19..).unwrap_or("");
        let residues: Vec<String> = residues_part
            .split_whitespace()
            .map(|s| s.to_string())
            .collect();

        let entry = chains.entry(chain_id).or_insert((0, Vec::new()));
        if entry.0 == 0 {
            entry.0 = num_res;
        }
        entry.1.extend(residues);
    }

    let mut records: Vec<SeqresRecord> = chains
        .into_iter()
        .map(|(chain_id, (num_residues, residues))| SeqresRecord {
            chain_id,
            num_residues,
            residues,
        })
        .collect();
    records.sort_by_key(|r| r.chain_id);
    records
}

/// Convert 3-letter amino acid code to 1-letter code.
fn aa3_to_1(aa: &str) -> char {
    match aa {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        _ => '?',
    }
}

// ---------------------------------------------------------------------------
// Secondary Structure Records
// ---------------------------------------------------------------------------

/// Type of secondary structure element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecStructType {
    /// Alpha helix (HELIX record).
    Helix,
    /// Beta sheet (SHEET record).
    Sheet,
    /// Turn (TURN record).
    Turn,
}

/// A secondary structure record (HELIX, SHEET, or TURN).
#[derive(Debug, Clone)]
pub struct SecStructRecord {
    /// Type of secondary structure.
    pub record_type: SecStructType,
    /// Chain identifier.
    pub chain_id: char,
    /// Starting residue sequence number.
    pub start_res: u32,
    /// Ending residue sequence number.
    pub end_res: u32,
    /// Starting residue name.
    pub start_res_name: String,
    /// Ending residue name.
    pub end_res_name: String,
}

/// Parse HELIX, SHEET, and TURN records from PDB content.
pub fn parse_secondary_structure(content: &str) -> Vec<SecStructRecord> {
    let mut records = Vec::new();

    for line in content.lines() {
        if line.starts_with("HELIX") {
            if line.len() < 26 {
                continue;
            }
            let chain_id = line.as_bytes().get(19).map_or('A', |&b| b as char);
            let start_res_name = line.get(15..18).unwrap_or("   ").trim().to_string();
            let start_res: u32 = line
                .get(21..25)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0);
            let end_res_name = line.get(27..30).unwrap_or("   ").trim().to_string();
            let end_res: u32 = line
                .get(33..37)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0);
            records.push(SecStructRecord {
                record_type: SecStructType::Helix,
                chain_id,
                start_res,
                end_res,
                start_res_name,
                end_res_name,
            });
        } else if line.starts_with("SHEET") {
            if line.len() < 22 {
                continue;
            }
            let chain_id = line.as_bytes().get(21).map_or('A', |&b| b as char);
            let start_res_name = line.get(17..20).unwrap_or("   ").trim().to_string();
            let start_res: u32 = line
                .get(22..26)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0);
            let end_res_name = line.get(28..31).unwrap_or("   ").trim().to_string();
            let end_res: u32 = line
                .get(33..37)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0);
            records.push(SecStructRecord {
                record_type: SecStructType::Sheet,
                chain_id,
                start_res,
                end_res,
                start_res_name,
                end_res_name,
            });
        } else if line.starts_with("TURN") {
            if line.len() < 22 {
                continue;
            }
            let chain_id = line.as_bytes().get(19).map_or('A', |&b| b as char);
            let start_res_name = line.get(15..18).unwrap_or("   ").trim().to_string();
            let start_res: u32 = line
                .get(20..24)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0);
            let end_res_name = line.get(26..29).unwrap_or("   ").trim().to_string();
            let end_res: u32 = line
                .get(30..34)
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or(0);
            records.push(SecStructRecord {
                record_type: SecStructType::Turn,
                chain_id,
                start_res,
                end_res,
                start_res_name,
                end_res_name,
            });
        }
    }

    records
}

// ---------------------------------------------------------------------------
// Biological Assembly
// ---------------------------------------------------------------------------

/// A chain of atoms within a biological assembly.
#[derive(Debug, Clone)]
pub struct AssemblyChain {
    /// Chain identifier.
    pub chain_id: char,
    /// All atoms in this chain.
    pub atoms: Vec<PdbAtom>,
}

/// A biological assembly: one or more copies of the asymmetric unit
/// arranged by crystallographic symmetry.
#[derive(Debug, Clone)]
pub struct BiologicalAssembly {
    /// Assembly number (as in REMARK 350).
    pub id: u32,
    /// All chains in this assembly.
    pub chains: Vec<AssemblyChain>,
}

impl BiologicalAssembly {
    /// Build a biological assembly from a `PdbModel` using the identity transform.
    pub fn from_model_with_identity(model: &PdbModel, id: u32) -> Self {
        // Group atoms by chain
        let mut chain_map: std::collections::HashMap<char, Vec<PdbAtom>> =
            std::collections::HashMap::new();
        for atom in &model.atoms {
            chain_map
                .entry(atom.chain_id)
                .or_default()
                .push(atom.clone());
        }
        let mut chains: Vec<AssemblyChain> = chain_map
            .into_iter()
            .map(|(chain_id, atoms)| AssemblyChain { chain_id, atoms })
            .collect();
        chains.sort_by_key(|c| c.chain_id);
        Self { id, chains }
    }

    /// Apply a 4×3 rotation-translation matrix to all atoms.
    ///
    /// Matrix format: `[[r00,r01,r02,tx\],[r10,r11,r12,ty],[r20,r21,r22,tz]]`
    pub fn apply_transform(&self, m: &[[f64; 4]; 3]) -> Self {
        let transform_atom = |atom: &PdbAtom| -> PdbAtom {
            let x = m[0][0] * atom.x + m[0][1] * atom.y + m[0][2] * atom.z + m[0][3];
            let y = m[1][0] * atom.x + m[1][1] * atom.y + m[1][2] * atom.z + m[1][3];
            let z = m[2][0] * atom.x + m[2][1] * atom.y + m[2][2] * atom.z + m[2][3];
            let mut a = atom.clone();
            a.x = x;
            a.y = y;
            a.z = z;
            a
        };
        let chains = self
            .chains
            .iter()
            .map(|c| AssemblyChain {
                chain_id: c.chain_id,
                atoms: c.atoms.iter().map(transform_atom).collect(),
            })
            .collect();
        Self {
            id: self.id,
            chains,
        }
    }

    /// Total number of atoms across all chains.
    pub fn total_atoms(&self) -> usize {
        self.chains.iter().map(|c| c.atoms.len()).sum()
    }

    /// All chain IDs present.
    pub fn chain_ids(&self) -> Vec<char> {
        self.chains.iter().map(|c| c.chain_id).collect()
    }
}

// ---------------------------------------------------------------------------
// Symmetry Operations
// ---------------------------------------------------------------------------

/// A crystallographic symmetry operation: rotation + translation.
///
/// Stored as a 3×4 matrix: `[[r00,r01,r02,tx\],[r10,r11,r12,ty],[r20,r21,r22,tz]]`.
#[derive(Debug, Clone)]
pub struct SymmetryOperation {
    /// The 3×4 transformation matrix.
    pub matrix: [[f64; 4]; 3],
    /// Optional human-readable description (e.g., "X,Y,Z").
    pub description: String,
}

impl SymmetryOperation {
    /// The identity operation.
    pub fn identity() -> Self {
        Self {
            matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
            description: "X,Y,Z".to_string(),
        }
    }

    /// A pure translation.
    pub fn translation(t: [f64; 3]) -> Self {
        Self {
            matrix: [
                [1.0, 0.0, 0.0, t[0]],
                [0.0, 1.0, 0.0, t[1]],
                [0.0, 0.0, 1.0, t[2]],
            ],
            description: format!("T({},{},{})", t[0], t[1], t[2]),
        }
    }

    /// Apply the symmetry operation to a 3D point.
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let m = &self.matrix;
        [
            m[0][0] * p[0] + m[0][1] * p[1] + m[0][2] * p[2] + m[0][3],
            m[1][0] * p[0] + m[1][1] * p[1] + m[1][2] * p[2] + m[1][3],
            m[2][0] * p[0] + m[2][1] * p[1] + m[2][2] * p[2] + m[2][3],
        ]
    }

    /// Apply the symmetry operation to a list of atoms.
    pub fn apply_to_atoms(&self, atoms: &[PdbAtom]) -> Vec<PdbAtom> {
        atoms
            .iter()
            .map(|a| {
                let [x, y, z] = self.apply([a.x, a.y, a.z]);
                let mut new_a = a.clone();
                new_a.x = x;
                new_a.y = y;
                new_a.z = z;
                new_a
            })
            .collect()
    }

    /// Generate all symmetry-equivalent copies of the atoms.
    pub fn generate_copies(ops: &[SymmetryOperation], atoms: &[PdbAtom]) -> Vec<Vec<PdbAtom>> {
        ops.iter().map(|op| op.apply_to_atoms(atoms)).collect()
    }
}

/// Parse REMARK 290 crystallographic symmetry operations from PDB content.
pub fn parse_symmetry_operations(content: &str) -> Vec<SymmetryOperation> {
    let mut ops = Vec::new();
    let mut in_remark_290 = false;
    let mut current_matrix: Vec<[f64; 4]> = Vec::new();

    for line in content.lines() {
        if line.starts_with("REMARK 290") {
            in_remark_290 = true;
            // Look for SMTRY lines
            if line.contains("SMTRY") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 8
                    && let (Ok(r1), Ok(r2), Ok(r3), Ok(t)) = (
                        parts[4].parse::<f64>(),
                        parts[5].parse::<f64>(),
                        parts[6].parse::<f64>(),
                        parts[7].parse::<f64>(),
                    )
                {
                    current_matrix.push([r1, r2, r3, t]);
                    if current_matrix.len() == 3 {
                        ops.push(SymmetryOperation {
                            matrix: [current_matrix[0], current_matrix[1], current_matrix[2]],
                            description: String::new(),
                        });
                        current_matrix.clear();
                    }
                }
            }
        } else if in_remark_290 && !line.starts_with("REMARK") {
            in_remark_290 = false;
        }
    }

    // Always include identity if no ops found
    if ops.is_empty() {
        ops.push(SymmetryOperation::identity());
    }
    ops
}

// ---------------------------------------------------------------------------
// PDB Validation
// ---------------------------------------------------------------------------

/// Result of validating a PDB model.
#[derive(Debug, Clone, Default)]
pub struct PdbValidationResult {
    /// Fatal errors (model should not be used as-is).
    pub errors: Vec<String>,
    /// Warnings (model is usable but has issues).
    pub warnings: Vec<String>,
}

impl PdbValidationResult {
    /// Whether validation passed with no errors.
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }

    /// Whether there are any issues (errors or warnings).
    pub fn has_issues(&self) -> bool {
        !self.errors.is_empty() || !self.warnings.is_empty()
    }
}

/// Validate a `PdbModel` for common errors and warnings.
///
/// Checks:
/// - Duplicate atom serial numbers (error)
/// - Atoms with extreme coordinates (warning)
/// - Missing occupancy (warning)
/// - Empty model (warning)
pub fn validate_pdb_model(model: &PdbModel) -> PdbValidationResult {
    let mut result = PdbValidationResult::default();

    if model.atoms.is_empty() {
        result.warnings.push("model has no atoms".to_string());
        return result;
    }

    // Check duplicate serial numbers
    let mut seen_serials = std::collections::HashSet::new();
    for atom in &model.atoms {
        if !seen_serials.insert(atom.serial) {
            result
                .errors
                .push(format!("duplicate atom serial number: {}", atom.serial));
        }
    }

    // Check extreme coordinates
    for atom in &model.atoms {
        let max_coord = atom.x.abs().max(atom.y.abs()).max(atom.z.abs());
        if max_coord > 9999.0 {
            result.warnings.push(format!(
                "atom {} has extreme coordinate ({:.1}, {:.1}, {:.1})",
                atom.serial, atom.x, atom.y, atom.z
            ));
        }
    }

    // Check occupancy values
    for atom in &model.atoms {
        if atom.occupancy < 0.0 || atom.occupancy > 1.0 {
            result.warnings.push(format!(
                "atom {} has out-of-range occupancy: {}",
                atom.serial, atom.occupancy
            ));
        }
    }

    // Check B-factors
    for atom in &model.atoms {
        if atom.temp_factor < 0.0 {
            result.warnings.push(format!(
                "atom {} has negative B-factor: {}",
                atom.serial, atom.temp_factor
            ));
        }
    }

    result
}

/// Validate a full multi-model PDB file.
pub fn validate_pdb_file(pdb_file: &PdbFile) -> Vec<PdbValidationResult> {
    pdb_file.models.iter().map(validate_pdb_model).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_atom(
        serial: u32,
        name: &str,
        res: &str,
        chain: char,
        x: f64,
        y: f64,
        z: f64,
    ) -> PdbAtom {
        PdbAtom::new(PdbAtomData {
            serial,
            name: name.to_string(),
            res_name: res.to_string(),
            chain_id: chain,
            res_seq: serial,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: String::new(),
        })
    }

    #[test]
    fn test_pdb_write_and_read_roundtrip() {
        let path = std::env::temp_dir().join("oxiphy_test.pdb");
        let atoms = vec![
            make_atom(1, "CA", "ALA", 'A', 1.234, 5.678, 9.012),
            make_atom(2, "N", "GLY", 'A', -1.0, 2.5, 3.75),
        ];
        PdbWriter::write_frame(path.to_str().unwrap_or(""), &atoms).unwrap();
        let read_atoms = PdbReader::read(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(read_atoms.len(), 2);
        assert!((read_atoms[0].x - 1.234).abs() < 0.01);
        assert!((read_atoms[1].z - 3.75).abs() < 0.01);
        assert_eq!(read_atoms[0].name, "CA");
        assert_eq!(read_atoms[1].res_name, "GLY");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_pdb_write_read_roundtrip() {
        let path = std::env::temp_dir().join("oxiphy_test_rt3.pdb");
        let atoms = vec![
            make_atom(1, "N", "ALA", 'A', 10.111, 20.222, 30.333),
            make_atom(2, "CA", "ALA", 'A', -5.555, 0.001, 15.999),
            make_atom(3, "C", "GLY", 'B', 0.0, 0.0, 0.0),
        ];
        PdbWriter::write_frame(path.to_str().unwrap_or(""), &atoms).unwrap();
        let read_atoms = PdbReader::read(path.to_str().unwrap_or("")).unwrap();

        assert_eq!(read_atoms.len(), 3, "expected 3 atoms");

        assert!((read_atoms[0].x - 10.111).abs() < 0.001);
        assert!((read_atoms[0].y - 20.222).abs() < 0.001);
        assert!((read_atoms[0].z - 30.333).abs() < 0.001);
        assert!((read_atoms[1].x - (-5.555)).abs() < 0.001);
        assert!((read_atoms[1].y - 0.001).abs() < 0.001);
        assert!((read_atoms[1].z - 15.999).abs() < 0.001);
        assert!((read_atoms[2].x - 0.0).abs() < 0.001);
        assert!((read_atoms[2].y - 0.0).abs() < 0.001);
        assert!((read_atoms[2].z - 0.0).abs() < 0.001);

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_pdb_column_widths() {
        let path = std::env::temp_dir().join("oxiphy_test_cols.pdb");
        let atoms = vec![make_atom(1, "CA", "ALA", 'A', 1.0, 2.0, 3.0)];
        PdbWriter::write_frame(path.to_str().unwrap_or(""), &atoms).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let line = content.lines().next().unwrap();
        assert!(line.starts_with("ATOM  "));
        let x_str = &line[30..38];
        assert!(x_str.trim().parse::<f64>().is_ok());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_parse_pdb_minimal() {
        let pdb_str = "\
TITLE     Test molecule
REMARK Test remark
ATOM      1  CA  ALA A   1       1.234   5.678   9.012  1.00  0.00           C
ATOM      2  N   GLY A   2      -1.000   2.500   3.750  1.00  0.00           N
END
";
        let model = parse_pdb(pdb_str).unwrap();
        assert_eq!(model.atoms.len(), 2);
        assert_eq!(model.title, "Test molecule");
        assert_eq!(model.remarks.len(), 1);

        let a0 = &model.atoms[0];
        assert_eq!(a0.name, "CA");
        assert_eq!(a0.res_name, "ALA");
        assert!((a0.x - 1.234).abs() < 0.001);
        assert!((a0.y - 5.678).abs() < 0.001);
        assert!((a0.z - 9.012).abs() < 0.001);

        let a1 = &model.atoms[1];
        assert_eq!(a1.name, "N");
        assert!((a1.x - (-1.0)).abs() < 0.001);
        assert!((a1.z - 3.75).abs() < 0.001);
    }

    #[test]
    fn test_parse_write_parse_roundtrip() {
        let pdb_str = "\
TITLE     Roundtrip test
ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00  5.00
ATOM      2  N   GLY B   2      -5.000   0.500  15.000  0.50 10.00
END
";
        let model1 = parse_pdb(pdb_str).unwrap();
        let written = write_pdb(&model1);
        let model2 = parse_pdb(&written).unwrap();

        assert_eq!(model1.atoms.len(), model2.atoms.len());
        for (a, b) in model1.atoms.iter().zip(model2.atoms.iter()) {
            assert!((a.x - b.x).abs() < 0.001);
            assert!((a.y - b.y).abs() < 0.001);
            assert!((a.z - b.z).abs() < 0.001);
        }
    }

    // ─── New tests ────────────────────────────────────────────────────

    #[test]
    fn test_hetatm_parsing() {
        let pdb_str = "\
HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00  0.00           O
END
";
        let model = parse_pdb(pdb_str).unwrap();
        assert_eq!(model.atoms.len(), 1);
        assert!(model.atoms[0].is_hetatm);
        assert_eq!(model.atoms[0].res_name, "HOH");
    }

    #[test]
    fn test_hetatm_write_roundtrip() {
        let mut atom = make_atom(1, "O", "HOH", 'A', 1.0, 2.0, 3.0);
        atom.is_hetatm = true;
        let model = PdbModel {
            atoms: vec![atom],
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let written = write_pdb(&model);
        assert!(written.contains("HETATM"));
        let parsed = parse_pdb(&written).unwrap();
        assert!(parsed.atoms[0].is_hetatm);
    }

    #[test]
    fn test_conect_record_parsing() {
        let pdb_str = "\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C
ATOM      2  N   ALA A   1       4.000   5.000   6.000  1.00  0.00           N
CONECT    1    2
END
";
        let model = parse_pdb(pdb_str).unwrap();
        assert_eq!(model.conect_records.len(), 1);
        assert_eq!(model.conect_records[0].atom_serial, 1);
        assert_eq!(model.conect_records[0].bonded_atoms, vec![2]);
    }

    #[test]
    fn test_conect_multiple_bonds() {
        let line = "CONECT    1    2    3    4";
        let conect = parse_conect_line(line).unwrap();
        assert_eq!(conect.atom_serial, 1);
        assert_eq!(conect.bonded_atoms, vec![2, 3, 4]);
    }

    #[test]
    fn test_conect_write_roundtrip() {
        let model = PdbModel {
            atoms: vec![
                make_atom(1, "CA", "ALA", 'A', 1.0, 0.0, 0.0),
                make_atom(2, "N", "ALA", 'A', 2.0, 0.0, 0.0),
            ],
            remarks: Vec::new(),
            title: String::new(),
            conect_records: vec![ConectRecord {
                atom_serial: 1,
                bonded_atoms: vec![2],
            }],
            cryst1: None,
            model_number: None,
        };
        let written = write_pdb(&model);
        assert!(written.contains("CONECT"));
        let parsed = parse_pdb(&written).unwrap();
        assert_eq!(parsed.conect_records.len(), 1);
        assert_eq!(parsed.conect_records[0].atom_serial, 1);
    }

    #[test]
    fn test_cryst1_parsing() {
        let pdb_str = "\
CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 1           1
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00
END
";
        let model = parse_pdb(pdb_str).unwrap();
        let cryst = model.cryst1.unwrap();
        assert!((cryst.a - 50.0).abs() < 0.01);
        assert!((cryst.b - 60.0).abs() < 0.01);
        assert!((cryst.c - 70.0).abs() < 0.01);
        assert!((cryst.alpha - 90.0).abs() < 0.01);
    }

    #[test]
    fn test_cryst1_write_roundtrip() {
        let model = PdbModel {
            atoms: vec![make_atom(1, "CA", "ALA", 'A', 1.0, 2.0, 3.0)],
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: Some(Cryst1Record {
                a: 50.0,
                b: 60.0,
                c: 70.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
                space_group: "P 1".to_string(),
                z_value: 1,
            }),
            model_number: None,
        };
        let written = write_pdb(&model);
        assert!(written.contains("CRYST1"));
        let parsed = parse_pdb(&written).unwrap();
        let cryst = parsed.cryst1.unwrap();
        assert!((cryst.a - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_multi_model_pdb() {
        let pdb_str = "\
TITLE     Multi model test
MODEL        1
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00
ENDMDL
MODEL        2
ATOM      1  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00
ENDMDL
END
";
        let pdb_file = parse_multi_model_pdb(pdb_str).unwrap();
        assert_eq!(pdb_file.models.len(), 2);
        assert_eq!(pdb_file.title, "Multi model test");
        assert!((pdb_file.models[0].atoms[0].x - 1.0).abs() < 0.01);
        assert!((pdb_file.models[1].atoms[0].x - 4.0).abs() < 0.01);
        assert_eq!(pdb_file.models[0].model_number, Some(1));
        assert_eq!(pdb_file.models[1].model_number, Some(2));
    }

    #[test]
    fn test_multi_model_write_roundtrip() {
        let pdb_file = PdbFile {
            models: vec![
                PdbModel {
                    atoms: vec![make_atom(1, "CA", "ALA", 'A', 1.0, 2.0, 3.0)],
                    remarks: Vec::new(),
                    title: String::new(),
                    conect_records: Vec::new(),
                    cryst1: None,
                    model_number: Some(1),
                },
                PdbModel {
                    atoms: vec![make_atom(1, "CA", "ALA", 'A', 4.0, 5.0, 6.0)],
                    remarks: Vec::new(),
                    title: String::new(),
                    conect_records: Vec::new(),
                    cryst1: None,
                    model_number: Some(2),
                },
            ],
            remarks: Vec::new(),
            title: "Test".to_string(),
            cryst1: None,
        };
        let written = write_multi_model_pdb(&pdb_file);
        assert!(written.contains("MODEL"));
        assert!(written.contains("ENDMDL"));
        let parsed = parse_multi_model_pdb(&written).unwrap();
        assert_eq!(parsed.models.len(), 2);
    }

    #[test]
    fn test_single_model_as_multi() {
        let pdb_str = "\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00
END
";
        let pdb_file = parse_multi_model_pdb(pdb_str).unwrap();
        assert_eq!(pdb_file.models.len(), 1);
    }

    #[test]
    fn test_occupancy_and_bfactor() {
        let pdb_str = "\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  0.75 25.00
END
";
        let model = parse_pdb(pdb_str).unwrap();
        assert!((model.atoms[0].occupancy - 0.75).abs() < 0.01);
        assert!((model.atoms[0].temp_factor - 25.0).abs() < 0.01);
    }

    #[test]
    fn test_occupancy_bfactor_roundtrip() {
        let atom = PdbAtom::new(PdbAtomData {
            serial: 1,
            name: "CA".into(),
            res_name: "ALA".into(),
            chain_id: 'A',
            res_seq: 1,
            x: 1.0,
            y: 2.0,
            z: 3.0,
            occupancy: 0.50,
            temp_factor: 30.0,
            element: "C".into(),
        });
        let model = PdbModel {
            atoms: vec![atom],
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let written = write_pdb(&model);
        let parsed = parse_pdb(&written).unwrap();
        assert!((parsed.atoms[0].occupancy - 0.50).abs() < 0.01);
        assert!((parsed.atoms[0].temp_factor - 30.0).abs() < 0.01);
    }

    #[test]
    fn test_filter_by_chain() {
        let atoms = vec![
            make_atom(1, "CA", "ALA", 'A', 1.0, 0.0, 0.0),
            make_atom(2, "CA", "GLY", 'B', 2.0, 0.0, 0.0),
            make_atom(3, "N", "ALA", 'A', 3.0, 0.0, 0.0),
        ];
        let chain_a = filter_by_chain(&atoms, 'A');
        assert_eq!(chain_a.len(), 2);
        let chain_b = filter_by_chain(&atoms, 'B');
        assert_eq!(chain_b.len(), 1);
        let chain_c = filter_by_chain(&atoms, 'C');
        assert!(chain_c.is_empty());
    }

    #[test]
    fn test_unique_chain_ids() {
        let atoms = vec![
            make_atom(1, "CA", "ALA", 'B', 0.0, 0.0, 0.0),
            make_atom(2, "CA", "ALA", 'A', 0.0, 0.0, 0.0),
            make_atom(3, "CA", "ALA", 'B', 0.0, 0.0, 0.0),
        ];
        let chains = unique_chain_ids(&atoms);
        assert_eq!(chains, vec!['A', 'B']);
    }

    #[test]
    fn test_unique_residue_names() {
        let atoms = vec![
            make_atom(1, "CA", "GLY", 'A', 0.0, 0.0, 0.0),
            make_atom(2, "CA", "ALA", 'A', 0.0, 0.0, 0.0),
            make_atom(3, "CA", "GLY", 'A', 0.0, 0.0, 0.0),
        ];
        let names = unique_residue_names(&atoms);
        assert_eq!(names, vec!["ALA", "GLY"]);
    }

    #[test]
    fn test_center_of_mass() {
        let atoms = vec![
            make_atom(1, "CA", "ALA", 'A', 0.0, 0.0, 0.0),
            make_atom(2, "CA", "ALA", 'A', 10.0, 0.0, 0.0),
        ];
        let com = center_of_mass(&atoms);
        assert!((com[0] - 5.0).abs() < 1e-10);
        assert!(com[1].abs() < 1e-10);
    }

    #[test]
    fn test_center_of_mass_empty() {
        let com = center_of_mass(&[]);
        assert_eq!(com, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_atom_distance() {
        let a = make_atom(1, "CA", "ALA", 'A', 0.0, 0.0, 0.0);
        let b = make_atom(2, "CA", "ALA", 'A', 3.0, 4.0, 0.0);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_cryst1_default() {
        let cryst = Cryst1Record::default();
        assert!((cryst.a - 1.0).abs() < 1e-10);
        assert!((cryst.alpha - 90.0).abs() < 1e-10);
        assert_eq!(cryst.space_group, "P 1");
    }

    #[test]
    fn test_conect_not_conect() {
        assert!(parse_conect_line("ATOM  stuff").is_none());
    }

    #[test]
    fn test_cryst1_not_cryst1() {
        assert!(parse_cryst1_line("ATOM  stuff").is_none());
    }

    #[test]
    fn test_parse_atom_line_too_short() {
        assert!(parse_atom_line("ATOM  short").is_none());
    }

    // ── SEQRES tests ─────────────────────────────────────────────────────

    #[test]
    fn test_seqres_parsing_basic() {
        let pdb_str = "\
SEQRES   1 A   10  ALA GLY VAL LEU ILE PRO PHE TRP MET SER
END
";
        let records = parse_seqres(pdb_str);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chain_id, 'A');
        assert_eq!(records[0].num_residues, 10);
        assert_eq!(records[0].residues.len(), 10);
        assert_eq!(records[0].residues[0], "ALA");
    }

    #[test]
    fn test_seqres_multi_chain() {
        let pdb_str = "\
SEQRES   1 A    5  ALA GLY VAL LEU ILE
SEQRES   1 B    3  TRP MET SER
END
";
        let records = parse_seqres(pdb_str);
        assert_eq!(records.len(), 2);
        let chains: Vec<char> = records.iter().map(|r| r.chain_id).collect();
        assert!(chains.contains(&'A'));
        assert!(chains.contains(&'B'));
    }

    #[test]
    fn test_seqres_sequence_string() {
        let pdb_str = "\
SEQRES   1 A    3  ALA GLY VAL
END
";
        let records = parse_seqres(pdb_str);
        let seq = records[0].sequence_string();
        assert_eq!(seq, "ALA-GLY-VAL");
    }

    // ── Secondary structure tests ──────────────────────────────────────

    #[test]
    fn test_helix_record_parsing() {
        let pdb_str = "\
HELIX    1   1 ALA A    1  GLY A   10  1                                  10
END
";
        let records = parse_secondary_structure(pdb_str);
        let helices: Vec<_> = records
            .iter()
            .filter(|r| r.record_type == SecStructType::Helix)
            .collect();
        assert_eq!(helices.len(), 1);
        assert_eq!(helices[0].chain_id, 'A');
        assert_eq!(helices[0].start_res, 1);
        assert_eq!(helices[0].end_res, 10);
    }

    #[test]
    fn test_sheet_record_parsing() {
        let pdb_str = "\
SHEET    1   A 2 ALA A   1  GLY A   5  0
END
";
        let records = parse_secondary_structure(pdb_str);
        let sheets: Vec<_> = records
            .iter()
            .filter(|r| r.record_type == SecStructType::Sheet)
            .collect();
        assert_eq!(sheets.len(), 1);
    }

    #[test]
    fn test_turn_record_parsing() {
        let pdb_str = "\
TURN     1 T1  ALA A    1  GLY A    4
END
";
        let records = parse_secondary_structure(pdb_str);
        let turns: Vec<_> = records
            .iter()
            .filter(|r| r.record_type == SecStructType::Turn)
            .collect();
        assert_eq!(turns.len(), 1);
    }

    // ── Biological assembly tests ──────────────────────────────────────

    #[test]
    fn test_biological_assembly_construction() {
        let atoms = vec![
            make_atom(1, "CA", "ALA", 'A', 1.0, 0.0, 0.0),
            make_atom(2, "CA", "GLY", 'A', 2.0, 0.0, 0.0),
        ];
        let model = PdbModel {
            atoms: atoms.clone(),
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let assembly = BiologicalAssembly::from_model_with_identity(&model, 1);
        assert_eq!(assembly.id, 1);
        assert_eq!(assembly.chains.len(), 1); // chain 'A'
        assert_eq!(assembly.chains[0].atoms.len(), 2);
    }

    #[test]
    fn test_biological_assembly_apply_transform() {
        let atoms = vec![make_atom(1, "CA", "ALA", 'A', 1.0, 0.0, 0.0)];
        let model = PdbModel {
            atoms,
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let assembly = BiologicalAssembly::from_model_with_identity(&model, 1);
        // Apply translation [1,0,0] → atom moves to [2,0,0]
        let transform = [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        let transformed = assembly.apply_transform(&transform);
        assert!((transformed.chains[0].atoms[0].x - 2.0).abs() < 1e-6);
    }

    // ── Symmetry operations tests ──────────────────────────────────────

    #[test]
    fn test_symmetry_operation_identity() {
        let sym = SymmetryOperation::identity();
        let pt = [1.0, 2.0, 3.0];
        let result = sym.apply(pt);
        assert!((result[0] - 1.0).abs() < 1e-10);
        assert!((result[1] - 2.0).abs() < 1e-10);
        assert!((result[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_symmetry_operation_translation() {
        let sym = SymmetryOperation::translation([5.0, 0.0, 0.0]);
        let pt = [1.0, 2.0, 3.0];
        let result = sym.apply(pt);
        assert!((result[0] - 6.0).abs() < 1e-10);
        assert!((result[1] - 2.0).abs() < 1e-10);
        assert!((result[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_symmetry_apply_to_atoms() {
        let atoms = vec![make_atom(1, "CA", "ALA", 'A', 1.0, 2.0, 3.0)];
        let sym = SymmetryOperation::translation([10.0, 0.0, 0.0]);
        let moved = sym.apply_to_atoms(&atoms);
        assert!((moved[0].x - 11.0).abs() < 1e-6);
    }

    // ── PDB validation tests ───────────────────────────────────────────

    #[test]
    fn test_pdb_validate_valid_model() {
        let atoms = vec![
            make_atom(1, "CA", "ALA", 'A', 1.0, 2.0, 3.0),
            make_atom(2, "N", "ALA", 'A', 2.0, 3.0, 4.0),
        ];
        let model = PdbModel {
            atoms,
            remarks: Vec::new(),
            title: "Test".to_string(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let result = validate_pdb_model(&model);
        assert!(result.is_valid(), "valid model should pass validation");
    }

    #[test]
    fn test_pdb_validate_duplicate_serial() {
        let mut a1 = make_atom(1, "CA", "ALA", 'A', 1.0, 0.0, 0.0);
        let mut a2 = make_atom(1, "N", "ALA", 'A', 2.0, 0.0, 0.0); // duplicate serial
        a1.serial = 1;
        a2.serial = 1;
        let model = PdbModel {
            atoms: vec![a1, a2],
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let result = validate_pdb_model(&model);
        assert!(
            !result.is_valid(),
            "duplicate serials should fail validation"
        );
        assert!(!result.errors.is_empty());
    }

    #[test]
    fn test_pdb_validate_extreme_coordinates() {
        let atom = PdbAtom {
            serial: 1,
            name: "CA".to_string(),
            res_name: "ALA".to_string(),
            chain_id: 'A',
            res_seq: 1,
            x: 1e10, // unreasonably large
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            is_hetatm: false,
        };
        let model = PdbModel {
            atoms: vec![atom],
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let result = validate_pdb_model(&model);
        assert!(
            result.warnings.iter().any(|w| w.contains("coordinate")),
            "extreme coordinate should generate a warning"
        );
    }

    #[test]
    fn test_pdb_validate_empty_model() {
        let model = PdbModel {
            atoms: Vec::new(),
            remarks: Vec::new(),
            title: String::new(),
            conect_records: Vec::new(),
            cryst1: None,
            model_number: None,
        };
        let result = validate_pdb_model(&model);
        assert!(
            result.warnings.iter().any(|w| w.contains("no atoms")),
            "empty model should warn about missing atoms"
        );
    }
}
