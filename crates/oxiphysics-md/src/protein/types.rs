//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// A protein chain composed of residues.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ProteinChain {
    /// Ordered list of residues in the chain.
    pub residues: Vec<Residue>,
}
impl ProteinChain {
    /// Create an empty protein chain.
    pub fn new() -> Self {
        ProteinChain {
            residues: Vec::new(),
        }
    }
    /// Append a residue with the given amino acid type and Cα position.
    pub fn add_residue(&mut self, aa: AminoAcid, ca_pos: [f64; 3]) {
        self.residues.push(Residue {
            aa,
            ca_position: ca_pos,
        });
    }
    /// Number of residues in the chain.
    pub fn length(&self) -> usize {
        self.residues.len()
    }
    /// Return the sequence of amino acids.
    pub fn sequence(&self) -> Vec<AminoAcid> {
        self.residues.iter().map(|r| r.aa).collect()
    }
    /// Total molecular weight (sum of residue MW + one water molecule).
    pub fn molecular_weight(&self) -> f64 {
        let sum: f64 = self.residues.iter().map(|r| r.aa.molecular_weight()).sum();
        let n = self.residues.len() as f64;
        if n < 1.0 {
            return 0.0;
        }
        sum - (n - 1.0) * 18.015 + 18.015
    }
    /// Radius of gyration of the Cα atoms in Angstroms.
    pub fn radius_of_gyration(&self) -> f64 {
        let n = self.residues.len();
        if n == 0 {
            return 0.0;
        }
        let cx: f64 = self.residues.iter().map(|r| r.ca_position[0]).sum::<f64>() / n as f64;
        let cy: f64 = self.residues.iter().map(|r| r.ca_position[1]).sum::<f64>() / n as f64;
        let cz: f64 = self.residues.iter().map(|r| r.ca_position[2]).sum::<f64>() / n as f64;
        let msd: f64 = self
            .residues
            .iter()
            .map(|r| {
                let dx = r.ca_position[0] - cx;
                let dy = r.ca_position[1] - cy;
                let dz = r.ca_position[2] - cz;
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>()
            / n as f64;
        msd.sqrt()
    }
    /// End-to-end distance between the first and last Cα atoms in Angstroms.
    pub fn end_to_end_distance(&self) -> f64 {
        if self.residues.len() < 2 {
            return 0.0;
        }
        let first = &self.residues[0].ca_position;
        let last = &self.residues[self.residues.len() - 1].ca_position;
        let dx = last[0] - first[0];
        let dy = last[1] - first[1];
        let dz = last[2] - first[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Fraction of hydrophobic residues in the chain.
    pub fn hydrophobic_fraction(&self) -> f64 {
        if self.residues.is_empty() {
            return 0.0;
        }
        let count = self
            .residues
            .iter()
            .filter(|r| r.aa.is_hydrophobic())
            .count();
        count as f64 / self.residues.len() as f64
    }
}
/// DSSP-like secondary structure classification.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecondaryStructure {
    /// Alpha-helix (H)
    Helix,
    /// Beta-strand (E)
    Strand,
    /// Turn (T)
    Turn,
    /// Coil / undefined (C)
    Coil,
}
/// A single residue with its amino acid type and Cα coordinates.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Residue {
    /// The amino acid type.
    pub aa: AminoAcid,
    /// Cα position in Angstroms `[x, y, z]`.
    pub ca_position: [f64; 3],
}
/// Extended residue with backbone dihedral angles.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ResidueExtended {
    /// Amino acid type.
    pub aa: AminoAcid,
    /// Cα position \[x, y, z\] in Angstrom.
    pub ca_position: [f64; 3],
    /// N backbone position (or None if unavailable).
    pub n_position: Option<[f64; 3]>,
    /// C backbone position (or None if unavailable).
    pub c_position: Option<[f64; 3]>,
    /// O backbone position (or None if unavailable).
    pub o_position: Option<[f64; 3]>,
}
/// Ramachandran plot statistics for a chain.
#[allow(dead_code)]
pub struct RamachandranStats {
    /// Fraction in alpha-helix region.
    pub helix_fraction: f64,
    /// Fraction in beta-sheet region.
    pub beta_fraction: f64,
    /// Fraction in left-helix region.
    pub left_helix_fraction: f64,
    /// Fraction in allowed region.
    pub allowed_fraction: f64,
    /// Fraction in disallowed region.
    pub disallowed_fraction: f64,
    /// Total number of classified residues.
    pub n_classified: usize,
}
/// The 20 standard amino acids.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AminoAcid {
    /// Alanine
    Ala,
    /// Glycine
    Gly,
    /// Valine
    Val,
    /// Leucine
    Leu,
    /// Isoleucine
    Ile,
    /// Proline
    Pro,
    /// Phenylalanine
    Phe,
    /// Tryptophan
    Trp,
    /// Methionine
    Met,
    /// Serine
    Ser,
    /// Threonine
    Thr,
    /// Cysteine
    Cys,
    /// Tyrosine
    Tyr,
    /// Histidine
    His,
    /// Aspartate
    Asp,
    /// Glutamate
    Glu,
    /// Asparagine
    Asn,
    /// Glutamine
    Gln,
    /// Lysine
    Lys,
    /// Arginine
    Arg,
}
impl AminoAcid {
    /// Parse from three-letter code (case-insensitive). Returns `None` if unknown.
    pub fn from_three_letter(s: &str) -> Option<AminoAcid> {
        match s.to_uppercase().as_str() {
            "ALA" => Some(AminoAcid::Ala),
            "GLY" => Some(AminoAcid::Gly),
            "VAL" => Some(AminoAcid::Val),
            "LEU" => Some(AminoAcid::Leu),
            "ILE" => Some(AminoAcid::Ile),
            "PRO" => Some(AminoAcid::Pro),
            "PHE" => Some(AminoAcid::Phe),
            "TRP" => Some(AminoAcid::Trp),
            "MET" => Some(AminoAcid::Met),
            "SER" => Some(AminoAcid::Ser),
            "THR" => Some(AminoAcid::Thr),
            "CYS" => Some(AminoAcid::Cys),
            "TYR" => Some(AminoAcid::Tyr),
            "HIS" => Some(AminoAcid::His),
            "ASP" => Some(AminoAcid::Asp),
            "GLU" => Some(AminoAcid::Glu),
            "ASN" => Some(AminoAcid::Asn),
            "GLN" => Some(AminoAcid::Gln),
            "LYS" => Some(AminoAcid::Lys),
            "ARG" => Some(AminoAcid::Arg),
            _ => None,
        }
    }
    /// Parse from one-letter code. Returns `None` if unknown.
    pub fn from_one_letter(c: char) -> Option<AminoAcid> {
        match c.to_uppercase().next().unwrap_or(' ') {
            'A' => Some(AminoAcid::Ala),
            'G' => Some(AminoAcid::Gly),
            'V' => Some(AminoAcid::Val),
            'L' => Some(AminoAcid::Leu),
            'I' => Some(AminoAcid::Ile),
            'P' => Some(AminoAcid::Pro),
            'F' => Some(AminoAcid::Phe),
            'W' => Some(AminoAcid::Trp),
            'M' => Some(AminoAcid::Met),
            'S' => Some(AminoAcid::Ser),
            'T' => Some(AminoAcid::Thr),
            'C' => Some(AminoAcid::Cys),
            'Y' => Some(AminoAcid::Tyr),
            'H' => Some(AminoAcid::His),
            'D' => Some(AminoAcid::Asp),
            'E' => Some(AminoAcid::Glu),
            'N' => Some(AminoAcid::Asn),
            'Q' => Some(AminoAcid::Gln),
            'K' => Some(AminoAcid::Lys),
            'R' => Some(AminoAcid::Arg),
            _ => None,
        }
    }
    /// Residue molecular weight in Daltons (residue mass, not free amino acid).
    pub fn molecular_weight(&self) -> f64 {
        match self {
            AminoAcid::Ala => 89.094,
            AminoAcid::Gly => 75.032,
            AminoAcid::Val => 117.148,
            AminoAcid::Leu => 131.175,
            AminoAcid::Ile => 131.175,
            AminoAcid::Pro => 115.132,
            AminoAcid::Phe => 165.192,
            AminoAcid::Trp => 204.228,
            AminoAcid::Met => 149.208,
            AminoAcid::Ser => 105.093,
            AminoAcid::Thr => 119.120,
            AminoAcid::Cys => 121.158,
            AminoAcid::Tyr => 181.191,
            AminoAcid::His => 155.156,
            AminoAcid::Asp => 133.104,
            AminoAcid::Glu => 147.130,
            AminoAcid::Asn => 132.119,
            AminoAcid::Gln => 146.146,
            AminoAcid::Lys => 146.189,
            AminoAcid::Arg => 174.203,
        }
    }
    /// Kyte-Doolittle hydrophobicity scale (range roughly -4.5 to +4.5).
    pub fn hydrophobicity(&self) -> f64 {
        match self {
            AminoAcid::Ala => 1.8,
            AminoAcid::Gly => -0.4,
            AminoAcid::Val => 4.2,
            AminoAcid::Leu => 3.8,
            AminoAcid::Ile => 4.5,
            AminoAcid::Pro => -1.6,
            AminoAcid::Phe => 2.8,
            AminoAcid::Trp => -0.9,
            AminoAcid::Met => 1.9,
            AminoAcid::Ser => -0.8,
            AminoAcid::Thr => -0.7,
            AminoAcid::Cys => 2.5,
            AminoAcid::Tyr => -1.3,
            AminoAcid::His => -3.2,
            AminoAcid::Asp => -3.5,
            AminoAcid::Glu => -3.5,
            AminoAcid::Asn => -3.5,
            AminoAcid::Gln => -3.5,
            AminoAcid::Lys => -3.9,
            AminoAcid::Arg => -4.5,
        }
    }
    /// Returns `true` if the Kyte-Doolittle hydrophobicity > 0.
    pub fn is_hydrophobic(&self) -> bool {
        self.hydrophobicity() > 0.0
    }
    /// Approximate formal charge at pH 7.
    pub fn charge_at_ph7(&self) -> f64 {
        match self {
            AminoAcid::Asp => -1.0,
            AminoAcid::Glu => -1.0,
            AminoAcid::Lys => 1.0,
            AminoAcid::Arg => 1.0,
            AminoAcid::His => 0.1,
            _ => 0.0,
        }
    }
}
/// Ramachandran plot region classification.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RamachandranRegion {
    /// Favoured alpha-helix region (phi ≈ -60°, psi ≈ -40°)
    AlphaHelix,
    /// Favoured beta-sheet region (phi ≈ -120°, psi ≈ 130°)
    BetaSheet,
    /// Left-handed helix region (phi ≈ 60°, psi ≈ 40°), typical for Gly
    LeftHelix,
    /// Generously allowed region
    AllowedRegion,
    /// Disallowed (outlier) region
    Disallowed,
}
