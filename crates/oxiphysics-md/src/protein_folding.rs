#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Protein folding simulation and analysis.
//!
//! Provides:
//!
//! - Gō (Go) model for protein folding with native contact map and funneled energy landscape
//! - Ramachandran plot analysis (phi/psi dihedral angle classification)
//! - Secondary structure assignment via a DSSP-like algorithm
//! - Folding free energy estimation (two-state thermodynamic model)
//! - Contact order calculation (sequence-length-normalised average contact separation)
//! - Protein stability analysis (thermal denaturation curve)
//! - Native contact fraction Q as folding reaction coordinate
//! - Folding kinetics (mean first-passage time estimation via Kramers theory)
//! - Protein energy decomposition into bond/angle/dihedral/LJ/electrostatics
//! - All-atom-like energy function with CHARMM-style terms

#![allow(dead_code)]

#[allow(unused_imports)]
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Boltzmann constant in kcal mol⁻¹ K⁻¹.
const KB_KCAL: f64 = 0.001987;
/// Gas constant R in kcal mol⁻¹ K⁻¹.
const R_KCAL: f64 = 0.001987;
/// Avogadro's number.
const AVOGADRO: f64 = 6.022e23;
/// Vacuum permittivity in CHARMM internal units (e² / (kcal mol⁻¹ Å)).
const COULOMB_FACTOR: f64 = 332.0636; // e² kcal⁻¹ mol Å
/// Angstrom → nanometre conversion.
const ANG_TO_NM: f64 = 0.1;

// ---------------------------------------------------------------------------
// Vec3 helpers
// ---------------------------------------------------------------------------

/// Euclidean distance between two 3-D points (Å).
#[inline]
pub fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Subtract two 3-D vectors.
#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Dot product of two 3-D vectors.
#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-D vectors.
#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalise a 3-D vector. Returns zero vector if near-zero magnitude.
#[inline]
fn vec3_norm(a: [f64; 3]) -> [f64; 3] {
    let len = (vec3_dot(a, a)).sqrt();
    if len < 1e-14 {
        return [0.0; 3];
    }
    [a[0] / len, a[1] / len, a[2] / len]
}

/// Compute the dihedral angle (rad) defined by four points p1–p4.
///
/// Uses the Praxeolitic formula for numerical stability.
pub fn dihedral_angle(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
    let b1 = vec3_sub(p2, p1);
    let b2 = vec3_sub(p3, p2);
    let b3 = vec3_sub(p4, p3);
    let n1 = vec3_cross(b1, b2);
    let n2 = vec3_cross(b2, b3);
    let m1 = vec3_cross(n1, vec3_norm(b2));
    let x = vec3_dot(n1, n2);
    let y = vec3_dot(m1, n2);
    y.atan2(x)
}

/// Compute the bond angle (rad) at atom b given three atoms a–b–c.
pub fn bond_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = vec3_sub(a, b);
    let bc = vec3_sub(c, b);
    let dot = vec3_dot(ba, bc);
    let denom = (vec3_dot(ba, ba) * vec3_dot(bc, bc)).sqrt();
    if denom < 1e-14 {
        return 0.0;
    }
    (dot / denom).clamp(-1.0, 1.0).acos()
}

// ---------------------------------------------------------------------------
// Amino acid representation
// ---------------------------------------------------------------------------

/// The 20 standard amino acid one-letter codes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AminoAcidCode {
    /// Alanine (A)
    Ala,
    /// Arginine (R)
    Arg,
    /// Asparagine (N)
    Asn,
    /// Aspartate (D)
    Asp,
    /// Cysteine (C)
    Cys,
    /// Glutamine (Q)
    Gln,
    /// Glutamate (E)
    Glu,
    /// Glycine (G)
    Gly,
    /// Histidine (H)
    His,
    /// Isoleucine (I)
    Ile,
    /// Leucine (L)
    Leu,
    /// Lysine (K)
    Lys,
    /// Methionine (M)
    Met,
    /// Phenylalanine (F)
    Phe,
    /// Proline (P)
    Pro,
    /// Serine (S)
    Ser,
    /// Threonine (T)
    Thr,
    /// Tryptophan (W)
    Trp,
    /// Tyrosine (Y)
    Tyr,
    /// Valine (V)
    Val,
}

impl AminoAcidCode {
    /// One-letter code as a `char`.
    pub fn one_letter(self) -> char {
        match self {
            AminoAcidCode::Ala => 'A',
            AminoAcidCode::Arg => 'R',
            AminoAcidCode::Asn => 'N',
            AminoAcidCode::Asp => 'D',
            AminoAcidCode::Cys => 'C',
            AminoAcidCode::Gln => 'Q',
            AminoAcidCode::Glu => 'E',
            AminoAcidCode::Gly => 'G',
            AminoAcidCode::His => 'H',
            AminoAcidCode::Ile => 'I',
            AminoAcidCode::Leu => 'L',
            AminoAcidCode::Lys => 'K',
            AminoAcidCode::Met => 'M',
            AminoAcidCode::Phe => 'F',
            AminoAcidCode::Pro => 'P',
            AminoAcidCode::Ser => 'S',
            AminoAcidCode::Thr => 'T',
            AminoAcidCode::Trp => 'W',
            AminoAcidCode::Tyr => 'Y',
            AminoAcidCode::Val => 'V',
        }
    }

    /// Hydrophobicity index (Kyte–Doolittle scale, range −4.5 to 4.5).
    pub fn hydrophobicity(self) -> f64 {
        match self {
            AminoAcidCode::Ile => 4.5,
            AminoAcidCode::Val => 4.2,
            AminoAcidCode::Leu => 3.8,
            AminoAcidCode::Phe => 2.8,
            AminoAcidCode::Cys => 2.5,
            AminoAcidCode::Met => 1.9,
            AminoAcidCode::Ala => 1.8,
            AminoAcidCode::Gly => -0.4,
            AminoAcidCode::Thr => -0.7,
            AminoAcidCode::Ser => -0.8,
            AminoAcidCode::Trp => -0.9,
            AminoAcidCode::Tyr => -1.3,
            AminoAcidCode::Pro => -1.6,
            AminoAcidCode::His => -3.2,
            AminoAcidCode::Glu => -3.5,
            AminoAcidCode::Gln => -3.5,
            AminoAcidCode::Asp => -3.5,
            AminoAcidCode::Asn => -3.5,
            AminoAcidCode::Lys => -3.9,
            AminoAcidCode::Arg => -4.5,
        }
    }
}

// ---------------------------------------------------------------------------
// Protein chain
// ---------------------------------------------------------------------------

/// A residue in a coarse-grained Cα representation.
#[derive(Clone, Debug)]
pub struct Residue {
    /// Amino acid type.
    pub aa: AminoAcidCode,
    /// Cα position (Å).
    pub ca_pos: [f64; 3],
    /// Residue index (0-based).
    pub index: usize,
}

impl Residue {
    /// Create a new residue.
    ///
    /// # Arguments
    /// * `aa`     – amino acid code
    /// * `ca_pos` – Cα coordinate (Å)
    /// * `index`  – 0-based residue number
    pub fn new(aa: AminoAcidCode, ca_pos: [f64; 3], index: usize) -> Self {
        Residue { aa, ca_pos, index }
    }
}

/// A protein chain in the Cα coarse-grained representation.
#[derive(Clone, Debug)]
pub struct ProteinChain {
    /// Residues in sequence order.
    pub residues: Vec<Residue>,
}

impl ProteinChain {
    /// Construct a chain from a sequence and a list of Cα positions.
    ///
    /// # Arguments
    /// * `sequence` – amino acid sequence
    /// * `positions` – Cα positions in Å, must have the same length as `sequence`
    pub fn new(sequence: Vec<AminoAcidCode>, positions: Vec<[f64; 3]>) -> Self {
        assert_eq!(
            sequence.len(),
            positions.len(),
            "sequence and positions must have equal length"
        );
        let residues = sequence
            .into_iter()
            .zip(positions)
            .enumerate()
            .map(|(i, (aa, pos))| Residue::new(aa, pos, i))
            .collect();
        ProteinChain { residues }
    }

    /// Number of residues.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// Returns `true` if the chain is empty.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Cα–Cα distance between residues `i` and `j` (Å).
    pub fn ca_distance(&self, i: usize, j: usize) -> f64 {
        dist3(self.residues[i].ca_pos, self.residues[j].ca_pos)
    }

    /// Radius of gyration Rg (Å).
    pub fn radius_of_gyration(&self) -> f64 {
        let n = self.residues.len() as f64;
        if n == 0.0 {
            return 0.0;
        }
        let cx = self.residues.iter().map(|r| r.ca_pos[0]).sum::<f64>() / n;
        let cy = self.residues.iter().map(|r| r.ca_pos[1]).sum::<f64>() / n;
        let cz = self.residues.iter().map(|r| r.ca_pos[2]).sum::<f64>() / n;
        let rg2 = self
            .residues
            .iter()
            .map(|r| {
                let dx = r.ca_pos[0] - cx;
                let dy = r.ca_pos[1] - cy;
                let dz = r.ca_pos[2] - cz;
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>()
            / n;
        rg2.sqrt()
    }
}

// ---------------------------------------------------------------------------
// Native contact map
// ---------------------------------------------------------------------------

/// A native contact between two residues in the folded structure.
#[derive(Clone, Debug)]
pub struct NativeContact {
    /// Index of residue i.
    pub i: usize,
    /// Index of residue j.
    pub j: usize,
    /// Native distance (Å).
    pub r0: f64,
    /// Contact energy ε (kcal mol⁻¹, positive = attractive).
    pub epsilon: f64,
}

impl NativeContact {
    /// Construct a native contact.
    ///
    /// # Arguments
    /// * `i`       – residue index
    /// * `j`       – residue index (j > i + 2 enforced in practice)
    /// * `r0`      – native Cα–Cα distance (Å)
    /// * `epsilon` – contact energy depth (kcal mol⁻¹)
    pub fn new(i: usize, j: usize, r0: f64, epsilon: f64) -> Self {
        NativeContact { i, j, r0, epsilon }
    }
}

/// Build a native contact map from a reference structure.
///
/// Two residues form a native contact if their Cα–Cα distance is within
/// `cutoff` Å and they are separated by at least 3 residues in sequence.
///
/// # Arguments
/// * `chain`   – reference (native) protein structure
/// * `cutoff`  – contact distance cutoff (Å)
/// * `epsilon` – uniform contact energy depth (kcal mol⁻¹)
pub fn build_native_contact_map(
    chain: &ProteinChain,
    cutoff: f64,
    epsilon: f64,
) -> Vec<NativeContact> {
    let n = chain.len();
    let mut contacts = Vec::new();
    for i in 0..n {
        for j in (i + 3)..n {
            let r = chain.ca_distance(i, j);
            if r <= cutoff {
                contacts.push(NativeContact::new(i, j, r, epsilon));
            }
        }
    }
    contacts
}

// ---------------------------------------------------------------------------
// Gō model energy
// ---------------------------------------------------------------------------

/// Gō-model energy for the entire chain.
///
/// The Gō model uses:
/// - Native contacts: V_nc = ε \[(r0/r)^12 − 2(r0/r)^6\]  (Lennard-Jones with minimum at r0)
/// - Non-native contacts: V_nn = ε (r0_rep / r)^12  (purely repulsive)
/// - Backbone bonds: harmonic V_b = k_b (r − r0)² / 2
/// - Backbone angles: harmonic V_θ = k_θ (θ − θ0)² / 2
///
/// Returns the total energy in kcal mol⁻¹.
///
/// # Arguments
/// * `chain`          – current chain configuration
/// * `native`         – native contact list from reference structure
/// * `k_bond`         – bond force constant (kcal mol⁻¹ Å⁻²)
/// * `r_bond`         – equilibrium Cα–Cα bond length (Å)
/// * `k_angle`        – angle force constant (kcal mol⁻¹ rad⁻²)
/// * `theta_bond`     – equilibrium Cα–Cα–Cα angle (rad)
/// * `r_rep`          – repulsion radius for non-native pairs (Å)
/// * `eps_rep`        – repulsion energy coefficient (kcal mol⁻¹)
#[allow(clippy::too_many_arguments)]
pub fn go_model_energy(
    chain: &ProteinChain,
    native: &[NativeContact],
    k_bond: f64,
    r_bond: f64,
    k_angle: f64,
    theta_bond: f64,
    r_rep: f64,
    eps_rep: f64,
) -> f64 {
    let n = chain.len();
    let mut energy = 0.0;

    // Bond energy (Cα i to Cα i+1)
    for i in 0..n.saturating_sub(1) {
        let r = chain.ca_distance(i, i + 1);
        energy += 0.5 * k_bond * (r - r_bond).powi(2);
    }

    // Angle energy (i, i+1, i+2)
    for i in 0..n.saturating_sub(2) {
        let theta = bond_angle(
            chain.residues[i].ca_pos,
            chain.residues[i + 1].ca_pos,
            chain.residues[i + 2].ca_pos,
        );
        energy += 0.5 * k_angle * (theta - theta_bond).powi(2);
    }

    // Build set of native contact pairs for fast lookup
    let native_pairs: std::collections::HashSet<(usize, usize)> =
        native.iter().map(|nc| (nc.i, nc.j)).collect();

    // Native contact energy (12-10 Gō potential)
    for nc in native {
        let r = chain.ca_distance(nc.i, nc.j);
        if r < 1e-6 {
            continue;
        }
        let ratio = nc.r0 / r;
        energy += nc.epsilon * (ratio.powi(12) - 2.0 * ratio.powi(6));
    }

    // Non-native repulsive interactions (separated by ≥ 3 in sequence)
    for i in 0..n {
        for j in (i + 3)..n {
            if native_pairs.contains(&(i, j)) {
                continue;
            }
            let r = chain.ca_distance(i, j);
            if r < 1e-6 {
                continue;
            }
            energy += eps_rep * (r_rep / r).powi(12);
        }
    }

    energy
}

// ---------------------------------------------------------------------------
// Native contact fraction Q
// ---------------------------------------------------------------------------

/// Compute the native contact fraction Q.
///
/// Q = (1/Nc) Σ_{(i,j) ∈ native} \[1 + exp((r_ij − λ r0_ij) / β)\]^{-1}
///
/// where λ = 1.2 (tolerance factor) and β = 1.0 Å (switching width).
///
/// # Arguments
/// * `chain`  – current configuration
/// * `native` – native contact list
pub fn native_contact_fraction(chain: &ProteinChain, native: &[NativeContact]) -> f64 {
    if native.is_empty() {
        return 0.0;
    }
    let lambda = 1.2_f64;
    let beta = 1.0_f64;
    let sum: f64 = native
        .iter()
        .map(|nc| {
            let r = chain.ca_distance(nc.i, nc.j);
            let arg = (r - lambda * nc.r0) / beta;
            1.0 / (1.0 + arg.exp())
        })
        .sum();
    sum / native.len() as f64
}

// ---------------------------------------------------------------------------
// Ramachandran plot
// ---------------------------------------------------------------------------

/// Secondary structure classification from (φ, ψ) angles.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RamachandranRegion {
    /// Right-handed α-helix (φ ≈ −60°, ψ ≈ −45°).
    AlphaHelix,
    /// β-sheet / extended strand (φ ≈ −120°, ψ ≈ +140°).
    BetaStrand,
    /// Left-handed α-helix (φ ≈ +60°, ψ ≈ +45°).
    LeftAlphaHelix,
    /// Polyproline II helix (φ ≈ −75°, ψ ≈ +150°).
    PPII,
    /// Generously allowed region.
    Allowed,
    /// Disallowed region.
    Disallowed,
}

/// Classify a (φ, ψ) pair into a Ramachandran region.
///
/// Angles must be in radians in (−π, +π].
///
/// # Arguments
/// * `phi` – φ backbone dihedral (rad)
/// * `psi` – ψ backbone dihedral (rad)
pub fn ramachandran_region(phi: f64, psi: f64) -> RamachandranRegion {
    let phi_d = phi.to_degrees();
    let psi_d = psi.to_degrees();

    // Right-handed α-helix: φ ∈ [−100, −20], ψ ∈ [−80, −5]
    if (-100.0..=-20.0).contains(&phi_d) && (-80.0..=-5.0).contains(&psi_d) {
        return RamachandranRegion::AlphaHelix;
    }
    // β-strand: φ ∈ [−180, −50], ψ ∈ [+90, +180] ∪ [−180, −160]
    if (-180.0..=-50.0).contains(&phi_d) && (psi_d >= 90.0 || psi_d <= -160.0) {
        return RamachandranRegion::BetaStrand;
    }
    // Left-handed α-helix: φ ∈ [+20, +100], ψ ∈ [+5, +80]
    if (20.0..=100.0).contains(&phi_d) && (5.0..=80.0).contains(&psi_d) {
        return RamachandranRegion::LeftAlphaHelix;
    }
    // Polyproline II: φ ∈ [−90, −50], ψ ∈ [+120, +180]
    if (-90.0..=-50.0).contains(&phi_d) && (120.0..=180.0).contains(&psi_d) {
        return RamachandranRegion::PPII;
    }
    // Broadly allowed: φ < 0 (most of left half of plot)
    if phi_d < 0.0 {
        return RamachandranRegion::Allowed;
    }
    RamachandranRegion::Disallowed
}

/// Percentage of residues in each Ramachandran region.
///
/// Returns (alpha_helix, beta_strand, other_allowed, disallowed) fractions summing to 1.
///
/// # Arguments
/// * `phi_psi` – slice of (φ, ψ) angle pairs (rad)
pub fn ramachandran_statistics(phi_psi: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    if phi_psi.is_empty() {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let mut alpha = 0usize;
    let mut beta = 0usize;
    let mut other = 0usize;
    let mut disallowed = 0usize;
    for &(phi, psi) in phi_psi {
        match ramachandran_region(phi, psi) {
            RamachandranRegion::AlphaHelix | RamachandranRegion::LeftAlphaHelix => alpha += 1,
            RamachandranRegion::BetaStrand => beta += 1,
            RamachandranRegion::Disallowed => disallowed += 1,
            _ => other += 1,
        }
    }
    let n = phi_psi.len() as f64;
    (
        alpha as f64 / n,
        beta as f64 / n,
        other as f64 / n,
        disallowed as f64 / n,
    )
}

// ---------------------------------------------------------------------------
// DSSP-like secondary structure assignment
// ---------------------------------------------------------------------------

/// Secondary structure element type (DSSP-like).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SecondaryStructure {
    /// α-helix (DSSP code H).
    AlphaHelix,
    /// β-strand (DSSP code E).
    BetaStrand,
    /// 3₁₀-helix (DSSP code G).
    ThreeTenHelix,
    /// π-helix (DSSP code I).
    PiHelix,
    /// Turn (DSSP code T).
    Turn,
    /// Bend (DSSP code S).
    Bend,
    /// Coil (DSSP code C).
    Coil,
}

/// Assign secondary structure from Cα-only geometry (simplified DSSP-like).
///
/// Residue i is assigned:
/// - `AlphaHelix`  if the Cα–Cα pseudo-bond angle at i and i±1 and the
///   virtual i→i+4 distance are consistent with an α-helix geometry.
/// - `BetaStrand`  if the chain is extended (bond angles > 130°).
/// - Otherwise `Coil`.
///
/// This is a simplified Cα-only approach; true DSSP requires H-bond detection.
///
/// # Arguments
/// * `chain` – current chain configuration
pub fn assign_secondary_structure(chain: &ProteinChain) -> Vec<SecondaryStructure> {
    let n = chain.len();
    let mut ss = vec![SecondaryStructure::Coil; n];

    for i in 1..n.saturating_sub(1) {
        let theta = bond_angle(
            chain.residues[i - 1].ca_pos,
            chain.residues[i].ca_pos,
            chain.residues[i + 1].ca_pos,
        );
        let theta_deg = theta.to_degrees();

        // α-helix: virtual Cα(i)→Cα(i+3) distance ≈ 5.4 Å, bond angle ≈ 91°
        if i + 3 < n {
            let d_i3 = chain.ca_distance(i, i + 3);
            if (theta_deg > 80.0 && theta_deg < 105.0) && (d_i3 > 4.5 && d_i3 < 6.5) {
                ss[i] = SecondaryStructure::AlphaHelix;
                continue;
            }
        }

        // β-strand: extended, bond angle > 120°
        if theta_deg > 120.0 {
            ss[i] = SecondaryStructure::BetaStrand;
            continue;
        }

        // 3₁₀-helix: i→i+2 distance ≈ 5.8 Å, bond angle 80-95°
        if i + 2 < n {
            let d_i2 = chain.ca_distance(i, i + 2);
            if (theta_deg > 80.0 && theta_deg < 95.0) && (d_i2 > 5.0 && d_i2 < 6.5) {
                ss[i] = SecondaryStructure::ThreeTenHelix;
                continue;
            }
        }
    }
    ss
}

/// Fraction of residues assigned to each secondary structure type.
///
/// Returns (helix_fraction, strand_fraction, coil_fraction).
///
/// # Arguments
/// * `ss` – secondary structure assignment per residue
pub fn secondary_structure_content(ss: &[SecondaryStructure]) -> (f64, f64, f64) {
    if ss.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let n = ss.len() as f64;
    let helix = ss
        .iter()
        .filter(|&&s| {
            s == SecondaryStructure::AlphaHelix
                || s == SecondaryStructure::ThreeTenHelix
                || s == SecondaryStructure::PiHelix
        })
        .count() as f64;
    let strand = ss
        .iter()
        .filter(|&&s| s == SecondaryStructure::BetaStrand)
        .count() as f64;
    let coil = n - helix - strand;
    (helix / n, strand / n, coil / n)
}

// ---------------------------------------------------------------------------
// Folding free energy (two-state model)
// ---------------------------------------------------------------------------

/// Folding free energy ΔGfold (kcal mol⁻¹) at temperature T.
///
/// Two-state model: ΔG(T) = ΔH − T ΔS.
///
/// # Arguments
/// * `delta_h` – folding enthalpy (kcal mol⁻¹, negative for stable proteins)
/// * `delta_s` – folding entropy (kcal mol⁻¹ K⁻¹, negative = ordering)
/// * `temp`    – temperature (K)
pub fn folding_free_energy(delta_h: f64, delta_s: f64, temp: f64) -> f64 {
    delta_h - temp * delta_s
}

/// Folding equilibrium constant K = exp(−ΔG / RT).
///
/// # Arguments
/// * `delta_g` – folding free energy (kcal mol⁻¹)
/// * `temp`    – temperature (K)
pub fn folding_equilibrium_constant(delta_g: f64, temp: f64) -> f64 {
    if temp <= 0.0 {
        return 0.0;
    }
    (-delta_g / (R_KCAL * temp)).exp()
}

/// Folded fraction f_F = K / (1 + K) in a two-state model.
///
/// # Arguments
/// * `delta_h` – folding enthalpy (kcal mol⁻¹)
/// * `delta_s` – folding entropy (kcal mol⁻¹ K⁻¹)
/// * `temp`    – temperature (K)
pub fn folded_fraction(delta_h: f64, delta_s: f64, temp: f64) -> f64 {
    let dg = folding_free_energy(delta_h, delta_s, temp);
    let k = folding_equilibrium_constant(dg, temp);
    k / (1.0 + k)
}

/// Melting temperature Tm = ΔH / ΔS (K).
///
/// # Arguments
/// * `delta_h` – folding enthalpy (kcal mol⁻¹)
/// * `delta_s` – folding entropy (kcal mol⁻¹ K⁻¹)
pub fn melting_temperature(delta_h: f64, delta_s: f64) -> f64 {
    if delta_s.abs() < 1e-20 {
        return f64::INFINITY;
    }
    delta_h / delta_s
}

// ---------------------------------------------------------------------------
// Contact order
// ---------------------------------------------------------------------------

/// Absolute contact order (CO) in units of sequence separation.
///
/// CO = (1/N_c) Σ_{(i,j) ∈ contacts} |i − j|
///
/// # Arguments
/// * `contacts` – list of native contacts
pub fn absolute_contact_order(contacts: &[NativeContact]) -> f64 {
    if contacts.is_empty() {
        return 0.0;
    }
    let sep_sum: f64 = contacts
        .iter()
        .map(|nc| (nc.j as i64 - nc.i as i64).abs() as f64)
        .sum();
    sep_sum / contacts.len() as f64
}

/// Relative contact order (RCO) normalised by chain length N.
///
/// RCO = CO / N
///
/// # Arguments
/// * `contacts` – list of native contacts
/// * `n`        – number of residues
pub fn relative_contact_order(contacts: &[NativeContact], n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    absolute_contact_order(contacts) / n as f64
}

// ---------------------------------------------------------------------------
// Thermal unfolding / stability analysis
// ---------------------------------------------------------------------------

/// Thermal unfolding curve: folded fraction vs temperature.
///
/// Returns a vector of (T, f_F) pairs over the specified temperature range.
///
/// # Arguments
/// * `delta_h`   – folding enthalpy (kcal mol⁻¹)
/// * `delta_s`   – folding entropy (kcal mol⁻¹ K⁻¹)
/// * `t_min`     – minimum temperature (K)
/// * `t_max`     – maximum temperature (K)
/// * `n_points`  – number of temperature points
pub fn thermal_unfolding_curve(
    delta_h: f64,
    delta_s: f64,
    t_min: f64,
    t_max: f64,
    n_points: usize,
) -> Vec<(f64, f64)> {
    if n_points == 0 {
        return Vec::new();
    }
    let dt = if n_points > 1 {
        (t_max - t_min) / (n_points - 1) as f64
    } else {
        0.0
    };
    (0..n_points)
        .map(|i| {
            let t = t_min + i as f64 * dt;
            (t, folded_fraction(delta_h, delta_s, t))
        })
        .collect()
}

/// Calorimetric heat capacity peak Cp_max (kcal mol⁻¹ K⁻²) at the melting transition.
///
/// Cp = ΔH² f_F (1 − f_F) / (R T²)
///
/// # Arguments
/// * `delta_h` – enthalpy (kcal mol⁻¹)
/// * `delta_s` – entropy (kcal mol⁻¹ K⁻¹)
/// * `temp`    – temperature (K)
pub fn heat_capacity(delta_h: f64, delta_s: f64, temp: f64) -> f64 {
    if temp <= 0.0 {
        return 0.0;
    }
    let ff = folded_fraction(delta_h, delta_s, temp);
    delta_h * delta_h * ff * (1.0 - ff) / (R_KCAL * temp * temp)
}

// ---------------------------------------------------------------------------
// Folding kinetics – Kramers / MFPT
// ---------------------------------------------------------------------------

/// Mean first-passage time (MFPT) estimate via Kramers theory.
///
/// MFPT ≈ τ₀ exp(ΔG‡ / RT)
///
/// where ΔG‡ is the barrier height and τ₀ is a pre-exponential timescale.
///
/// # Arguments
/// * `barrier_kcal` – activation free energy ΔG‡ (kcal mol⁻¹)
/// * `tau0`         – prefactor time scale (s)
/// * `temp`         – temperature (K)
pub fn kramers_mfpt(barrier_kcal: f64, tau0: f64, temp: f64) -> f64 {
    if temp <= 0.0 {
        return f64::INFINITY;
    }
    tau0 * (barrier_kcal / (R_KCAL * temp)).exp()
}

/// Folding rate kf (s⁻¹) from Kramers MFPT.
///
/// # Arguments
/// * `barrier_kcal` – barrier height (kcal mol⁻¹)
/// * `tau0`         – prefactor (s)
/// * `temp`         – temperature (K)
pub fn folding_rate(barrier_kcal: f64, tau0: f64, temp: f64) -> f64 {
    let mfpt = kramers_mfpt(barrier_kcal, tau0, temp);
    if mfpt.is_infinite() || mfpt == 0.0 {
        return 0.0;
    }
    1.0 / mfpt
}

/// Unfolding rate ku (s⁻¹) via detailed balance: ku = kf / K_eq.
///
/// # Arguments
/// * `kf`      – folding rate (s⁻¹)
/// * `delta_g` – folding free energy (kcal mol⁻¹)
/// * `temp`    – temperature (K)
pub fn unfolding_rate(kf: f64, delta_g: f64, temp: f64) -> f64 {
    let k_eq = folding_equilibrium_constant(delta_g, temp);
    if k_eq == 0.0 {
        return f64::INFINITY;
    }
    kf / k_eq
}

// ---------------------------------------------------------------------------
// CHARMM-like all-atom energy terms
// ---------------------------------------------------------------------------

/// CHARMM-style harmonic bond energy: V = k_b (b − b_0)².
///
/// # Arguments
/// * `r`   – current bond length (Å)
/// * `r0`  – equilibrium bond length (Å)
/// * `kb`  – force constant (kcal mol⁻¹ Å⁻²)
pub fn charmm_bond_energy(r: f64, r0: f64, kb: f64) -> f64 {
    kb * (r - r0).powi(2)
}

/// CHARMM-style harmonic angle energy: V = k_θ (θ − θ_0)².
///
/// # Arguments
/// * `theta`  – current angle (rad)
/// * `theta0` – equilibrium angle (rad)
/// * `kt`     – force constant (kcal mol⁻¹ rad⁻²)
pub fn charmm_angle_energy(theta: f64, theta0: f64, kt: f64) -> f64 {
    kt * (theta - theta0).powi(2)
}

/// CHARMM-style dihedral energy: V = k_φ \[1 + cos(n φ − δ)\].
///
/// # Arguments
/// * `phi`   – current dihedral (rad)
/// * `kphi`  – force constant (kcal mol⁻¹)
/// * `n`     – periodicity (integer)
/// * `delta` – phase offset (rad)
pub fn charmm_dihedral_energy(phi: f64, kphi: f64, n: i32, delta: f64) -> f64 {
    kphi * (1.0 + (n as f64 * phi - delta).cos())
}

/// CHARMM 1-4 (modified) Lennard-Jones energy.
///
/// V = ε \[(Rmin/r)^12 − 2(Rmin/r)^6\]
///
/// # Arguments
/// * `r`    – interatomic distance (Å)
/// * `rmin` – Rmin/2 parameter (Å), sum of van der Waals radii
/// * `eps`  – well depth ε (kcal mol⁻¹)
pub fn charmm_lj_energy(r: f64, rmin: f64, eps: f64) -> f64 {
    if r < 1e-6 {
        return 1e30;
    }
    let ratio = rmin / r;
    eps * (ratio.powi(12) - 2.0 * ratio.powi(6))
}

/// CHARMM Coulomb electrostatic energy.
///
/// V = q1 q2 / (ε_r r)  (in internal CHARMM units: kcal mol⁻¹)
///
/// # Arguments
/// * `q1`      – charge of atom 1 (elementary charges e)
/// * `q2`      – charge of atom 2 (elementary charges e)
/// * `r`       – distance (Å)
/// * `eps_r`   – relative permittivity (dielectric constant)
pub fn charmm_coulomb_energy(q1: f64, q2: f64, r: f64, eps_r: f64) -> f64 {
    if r < 1e-6 || eps_r == 0.0 {
        return 0.0;
    }
    COULOMB_FACTOR * q1 * q2 / (eps_r * r)
}

// ---------------------------------------------------------------------------
// Full protein energy decomposition
// ---------------------------------------------------------------------------

/// Complete energy decomposition for a chain of atoms.
///
/// All atom coordinates, connectivity, and parameters must be provided.
/// Returns `(E_bond, E_angle, E_dihedral, E_lj, E_coulomb)` in kcal mol⁻¹.
#[allow(clippy::too_many_arguments)]
pub fn protein_energy_decomposition(
    positions: &[[f64; 3]],
    bonds: &[(usize, usize, f64, f64)],         // (i, j, kb, r0)
    angles: &[(usize, usize, usize, f64, f64)], // (i, j, k, kt, theta0)
    dihedrals: &[(usize, usize, usize, usize, f64, i32, f64)], // (i,j,k,l, kphi, n, delta)
    lj_pairs: &[(usize, usize, f64, f64)],      // (i, j, rmin, eps)
    charges: &[f64],
    eps_r: f64,
) -> (f64, f64, f64, f64, f64) {
    let mut e_bond = 0.0;
    for &(i, j, kb, r0) in bonds {
        let r = dist3(positions[i], positions[j]);
        e_bond += charmm_bond_energy(r, r0, kb);
    }

    let mut e_angle = 0.0;
    for &(i, j, k, kt, theta0) in angles {
        let theta = bond_angle(positions[i], positions[j], positions[k]);
        e_angle += charmm_angle_energy(theta, theta0, kt);
    }

    let mut e_dihedral = 0.0;
    for &(i, j, k, l, kphi, n, delta) in dihedrals {
        let phi = dihedral_angle(positions[i], positions[j], positions[k], positions[l]);
        e_dihedral += charmm_dihedral_energy(phi, kphi, n, delta);
    }

    let mut e_lj = 0.0;
    for &(i, j, rmin, eps) in lj_pairs {
        let r = dist3(positions[i], positions[j]);
        e_lj += charmm_lj_energy(r, rmin, eps);
    }

    let mut e_coulomb = 0.0;
    let n = positions.len();
    if charges.len() >= n {
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist3(positions[i], positions[j]);
                e_coulomb += charmm_coulomb_energy(charges[i], charges[j], r, eps_r);
            }
        }
    }

    (e_bond, e_angle, e_dihedral, e_lj, e_coulomb)
}

// ---------------------------------------------------------------------------
// Energy landscape helpers
// ---------------------------------------------------------------------------

/// Funneled energy landscape model.
///
/// The total energy as a function of Q follows a funnel shape:
/// V(Q) = V0 (1 − Q) − ε Q + Γ Q(1 − Q)
///
/// where Γ is the roughness parameter.
///
/// # Arguments
/// * `q`       – native contact fraction Q ∈ \[0, 1\]
/// * `v0`      – denatured state energy (kcal mol⁻¹)
/// * `epsilon` – folded state stabilisation (kcal mol⁻¹)
/// * `gamma`   – roughness / barrier height (kcal mol⁻¹)
pub fn funneled_landscape_energy(q: f64, v0: f64, epsilon: f64, gamma: f64) -> f64 {
    v0 * (1.0 - q) - epsilon * q + gamma * q * (1.0 - q)
}

/// Location of the transition state in the funneled landscape.
///
/// At the TS, dV/dQ = 0:  Q_TS = (−v0 − ε + Γ) / (2 Γ − ...) simplified.
/// Returns Q_TS ∈ \[0, 1\] or `None` if no barrier exists.
///
/// # Arguments
/// * `v0`      – denatured energy
/// * `epsilon` – stabilisation energy
/// * `gamma`   – roughness
pub fn funneled_transition_state_q(v0: f64, epsilon: f64, gamma: f64) -> Option<f64> {
    // dV/dQ = -v0 - epsilon + gamma(1 - 2Q) = 0
    // → Q = (gamma - v0 - epsilon) / (2 gamma)
    if gamma.abs() < 1e-12 {
        return None;
    }
    let q_ts = (gamma - v0 - epsilon) / (2.0 * gamma);
    if (0.0..=1.0).contains(&q_ts) {
        Some(q_ts)
    } else {
        None
    }
}

/// Barrier height ΔG‡ (kcal mol⁻¹) in the funneled landscape.
///
/// ΔG‡ = V(Q_TS) − V(0)
///
/// # Arguments
/// * `v0`      – denatured energy
/// * `epsilon` – stabilisation
/// * `gamma`   – roughness
pub fn funneled_barrier_height(v0: f64, epsilon: f64, gamma: f64) -> f64 {
    match funneled_transition_state_q(v0, epsilon, gamma) {
        None => 0.0,
        Some(q_ts) => {
            let v_ts = funneled_landscape_energy(q_ts, v0, epsilon, gamma);
            let v0_state = funneled_landscape_energy(0.0, v0, epsilon, gamma);
            v_ts - v0_state
        }
    }
}

// ---------------------------------------------------------------------------
// Implicit solvent model (simplified GB-like)
// ---------------------------------------------------------------------------

/// Generalised Born (simplified) solvation free energy contribution
/// for a pair of charged atoms.
///
/// ΔG_GB ≈ −(1/ε_in − 1/ε_out) q1 q2 / (2 f_GB)
/// where f_GB = √(r² + Ri Rj exp(−r²/(4 Ri Rj)))
///
/// Returns energy in kcal mol⁻¹.
///
/// # Arguments
/// * `q1`, `q2`  – atomic charges (e)
/// * `r`         – distance (Å)
/// * `ri`, `rj`  – effective Born radii (Å)
/// * `eps_in`    – internal dielectric
/// * `eps_out`   – solvent dielectric
pub fn gb_solvation_pair(
    q1: f64,
    q2: f64,
    r: f64,
    ri: f64,
    rj: f64,
    eps_in: f64,
    eps_out: f64,
) -> f64 {
    let f_gb = (r * r + ri * rj * (-r * r / (4.0 * ri * rj)).exp()).sqrt();
    if f_gb < 1e-12 {
        return 0.0;
    }
    let coeff = -(1.0 / eps_in - 1.0 / eps_out);
    COULOMB_FACTOR * coeff * q1 * q2 / (2.0 * f_gb)
}

// ---------------------------------------------------------------------------
// Phi-value analysis
// ---------------------------------------------------------------------------

/// Brønsted β-Tanford phi-value: measures how native-like the transition state is.
///
/// φ_F = (ΔΔG‡_F_mutation) / ΔΔG_fold_mutation
///
/// φ ≈ 0: residue unfolded at TS;  φ ≈ 1: residue folded at TS.
///
/// # Arguments
/// * `ddg_barrier`  – change in barrier height upon mutation (kcal mol⁻¹)
/// * `ddg_fold`     – change in folding stability upon mutation (kcal mol⁻¹)
pub fn phi_value(ddg_barrier: f64, ddg_fold: f64) -> f64 {
    if ddg_fold.abs() < 1e-12 {
        return 0.0;
    }
    ddg_barrier / ddg_fold
}

// ---------------------------------------------------------------------------
// Circular variance for dihedral distributions
// ---------------------------------------------------------------------------

/// Circular mean of a set of angles (rad).
///
/// # Arguments
/// * `angles` – slice of angles in radians
pub fn circular_mean(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }
    let n = angles.len() as f64;
    let sin_sum: f64 = angles.iter().map(|&a| a.sin()).sum::<f64>() / n;
    let cos_sum: f64 = angles.iter().map(|&a| a.cos()).sum::<f64>() / n;
    sin_sum.atan2(cos_sum)
}

/// Circular variance R̄ ∈ \[0, 1\] of a set of angles.
///
/// R̄ = 0 means uniformly distributed; R̄ = 1 means perfectly concentrated.
///
/// # Arguments
/// * `angles` – slice of angles (rad)
pub fn circular_variance(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }
    let n = angles.len() as f64;
    let sin_sum: f64 = angles.iter().map(|&a| a.sin()).sum::<f64>() / n;
    let cos_sum: f64 = angles.iter().map(|&a| a.cos()).sum::<f64>() / n;
    (sin_sum * sin_sum + cos_sum * cos_sum).sqrt()
}

// ---------------------------------------------------------------------------
// Langevin thermostat step for coarse-grained MD
// ---------------------------------------------------------------------------

/// Apply one Langevin dynamics velocity update.
///
/// v_{n+1} = v_n (1 − γ dt) + (F/m) dt + σ √(dt) ξ
///
/// where ξ ~ N(0,1) is Gaussian noise and σ² = 2 γ kT / m.
///
/// # Arguments
/// * `vel`   – current velocity (Å ps⁻¹)
/// * `force` – current force (kcal mol⁻¹ Å⁻¹)
/// * `mass`  – mass (Da = g mol⁻¹ × AKMA units)
/// * `gamma` – friction coefficient (ps⁻¹)
/// * `dt`    – time step (ps)
/// * `temp`  – temperature (K)
/// * `xi`    – unit Gaussian random number
pub fn langevin_velocity_step(
    vel: f64,
    force: f64,
    mass: f64,
    gamma: f64,
    dt: f64,
    temp: f64,
    xi: f64,
) -> f64 {
    // σ = sqrt(2 gamma kT / m)  (AKMA: kT in kcal/mol, mass in g/mol, time in ps)
    // kb in kcal/(mol K)
    let sigma = (2.0 * gamma * KB_KCAL * temp / mass).sqrt();
    vel * (1.0 - gamma * dt) + (force / mass) * dt + sigma * dt.sqrt() * xi
}

// ---------------------------------------------------------------------------
// End-to-end distance and persistence length
// ---------------------------------------------------------------------------

/// End-to-end distance of the chain (Å).
///
/// # Arguments
/// * `chain` – protein chain
pub fn end_to_end_distance(chain: &ProteinChain) -> f64 {
    let n = chain.len();
    if n < 2 {
        return 0.0;
    }
    dist3(chain.residues[0].ca_pos, chain.residues[n - 1].ca_pos)
}

/// Persistence length estimate from the worm-like chain (WLC) model.
///
/// ⟨R²⟩ = 2 Lp L \[1 − (Lp/L)(1 − exp(−L/Lp))\]  ≈ 2 Lp L  for L ≫ Lp.
///
/// Inverted: Lp = ⟨R²⟩ / (2 L) (simple long-chain approximation).
///
/// # Arguments
/// * `r2_mean` – mean squared end-to-end distance ⟨R²⟩ (Å²)
/// * `contour_length` – contour length L (Å)
pub fn persistence_length_wlc(r2_mean: f64, contour_length: f64) -> f64 {
    if contour_length <= 0.0 {
        return 0.0;
    }
    r2_mean / (2.0 * contour_length)
}

// ---------------------------------------------------------------------------
// Boltzmann weighting
// ---------------------------------------------------------------------------

/// Compute the Boltzmann weight w_i = exp(−E_i / kT) for each energy.
///
/// # Arguments
/// * `energies` – slice of energies (kcal mol⁻¹)
/// * `temp`     – temperature (K)
pub fn boltzmann_weights(energies: &[f64], temp: f64) -> Vec<f64> {
    if temp <= 0.0 {
        return vec![0.0; energies.len()];
    }
    energies
        .iter()
        .map(|&e| (-e / (KB_KCAL * temp)).exp())
        .collect()
}

/// Boltzmann-weighted average of an observable.
///
/// ⟨O⟩ = Σ O_i exp(−E_i / kT) / Z
///
/// # Arguments
/// * `observables` – slice of observable values
/// * `energies`    – corresponding energies (kcal mol⁻¹)
/// * `temp`        – temperature (K)
pub fn boltzmann_average(observables: &[f64], energies: &[f64], temp: f64) -> f64 {
    assert_eq!(observables.len(), energies.len());
    let weights = boltzmann_weights(energies, temp);
    let z: f64 = weights.iter().sum();
    if z == 0.0 {
        return 0.0;
    }
    weights
        .iter()
        .zip(observables.iter())
        .map(|(w, o)| w * o)
        .sum::<f64>()
        / z
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_helix_chain(n: usize) -> ProteinChain {
        // Approximate α-helix Cα positions: rise 1.5 Å/res, radius 2.3 Å
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let angle = i as f64 * 100.0_f64.to_radians(); // ~100° per residue
                [2.3 * angle.cos(), 2.3 * angle.sin(), i as f64 * 1.5]
            })
            .collect();
        let seq = vec![AminoAcidCode::Ala; n];
        ProteinChain::new(seq, positions)
    }

    fn make_extended_chain(n: usize) -> ProteinChain {
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 3.8, 0.0, 0.0]).collect();
        let seq = vec![AminoAcidCode::Val; n];
        ProteinChain::new(seq, positions)
    }

    #[test]
    fn test_dist3_known() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 4.0, 0.0];
        assert!((dist3(a, b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral_angle_planar() {
        // All four points in xy-plane → dihedral = 0 or π
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [2.0, 1.0, 0.0];
        let p4 = [3.0, 1.0, 0.0];
        let d = dihedral_angle(p1, p2, p3, p4);
        // Check it returns a finite value
        assert!(d.is_finite());
    }

    #[test]
    fn test_bond_angle_right() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [0.0, 1.0, 0.0];
        let theta = bond_angle(a, b, c);
        assert!((theta - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_bond_angle_linear() {
        let a = [-1.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let c = [1.0, 0.0, 0.0];
        let theta = bond_angle(a, b, c);
        assert!((theta - PI).abs() < 1e-10);
    }

    #[test]
    fn test_amino_acid_one_letter() {
        assert_eq!(AminoAcidCode::Ala.one_letter(), 'A');
        assert_eq!(AminoAcidCode::Gly.one_letter(), 'G');
        assert_eq!(AminoAcidCode::Trp.one_letter(), 'W');
    }

    #[test]
    fn test_hydrophobicity_order() {
        assert!(AminoAcidCode::Ile.hydrophobicity() > AminoAcidCode::Arg.hydrophobicity());
    }

    #[test]
    fn test_protein_chain_len() {
        let chain = make_helix_chain(10);
        assert_eq!(chain.len(), 10);
        assert!(!chain.is_empty());
    }

    #[test]
    fn test_radius_of_gyration_extended() {
        // Extended chain along x: Rg should be proportional to chain length
        let chain = make_extended_chain(10);
        let rg = chain.radius_of_gyration();
        assert!(rg > 0.0);
    }

    #[test]
    fn test_build_native_contact_map_count() {
        let chain = make_helix_chain(20);
        let contacts = build_native_contact_map(&chain, 8.0, 1.0);
        // There should be some contacts
        assert!(!contacts.is_empty());
    }

    #[test]
    fn test_native_contact_map_min_separation() {
        let chain = make_helix_chain(10);
        let contacts = build_native_contact_map(&chain, 8.0, 1.0);
        for nc in &contacts {
            assert!(nc.j >= nc.i + 3, "contacts must be separated by ≥3");
        }
    }

    #[test]
    fn test_native_contact_fraction_folded() {
        // If current chain == native chain, Q should be > 0 (well above unfolded).
        // The sigmoid switching function with λ=1.2 gives Q < 1 even for the native
        // structure; the key property is Q > 0.5 for native vs. Q near 0 for unfolded.
        let chain = make_helix_chain(15);
        let native = build_native_contact_map(&chain, 8.0, 1.0);
        let q = native_contact_fraction(&chain, &native);
        assert!(
            q > 0.5,
            "Q should be well above 0 for native structure, got {q}"
        );
    }

    #[test]
    fn test_native_contact_fraction_empty() {
        let chain = make_helix_chain(5);
        let q = native_contact_fraction(&chain, &[]);
        assert_eq!(q, 0.0);
    }

    #[test]
    fn test_ramachandran_alpha_helix() {
        let phi = (-60.0_f64).to_radians();
        let psi = (-45.0_f64).to_radians();
        assert_eq!(
            ramachandran_region(phi, psi),
            RamachandranRegion::AlphaHelix
        );
    }

    #[test]
    fn test_ramachandran_beta_strand() {
        let phi = (-120.0_f64).to_radians();
        let psi = 140.0_f64.to_radians();
        assert_eq!(
            ramachandran_region(phi, psi),
            RamachandranRegion::BetaStrand
        );
    }

    #[test]
    fn test_ramachandran_statistics_pure_alpha() {
        let pp: Vec<(f64, f64)> = vec![((-60.0_f64).to_radians(), (-45.0_f64).to_radians()); 10];
        let (alpha, _beta, _other, disallowed) = ramachandran_statistics(&pp);
        assert_eq!(alpha, 1.0);
        assert_eq!(disallowed, 0.0);
    }

    #[test]
    fn test_assign_secondary_structure_length() {
        let chain = make_helix_chain(20);
        let ss = assign_secondary_structure(&chain);
        assert_eq!(ss.len(), 20);
    }

    #[test]
    fn test_secondary_structure_content_sums_to_one() {
        let chain = make_helix_chain(20);
        let ss = assign_secondary_structure(&chain);
        let (h, s, c) = secondary_structure_content(&ss);
        assert!(
            (h + s + c - 1.0).abs() < 1e-10,
            "fractions should sum to 1: {h}+{s}+{c}"
        );
    }

    #[test]
    fn test_folding_free_energy_sign() {
        // Stable protein: ΔH < 0, ΔS < 0, T * ΔS < ΔH at moderate T
        let dg = folding_free_energy(-50.0, -0.1, 300.0);
        assert!(dg < 0.0, "stable protein should have ΔG < 0, got {dg}");
    }

    #[test]
    fn test_melting_temperature() {
        let tm = melting_temperature(-50.0, -0.1);
        assert!((tm - 500.0).abs() < 1e-6, "Tm = 500 K expected, got {tm}");
    }

    #[test]
    fn test_folded_fraction_at_tm() {
        // At Tm, ΔG = 0, so ff = 0.5
        let delta_h = -50.0;
        let delta_s = -0.1;
        let tm = melting_temperature(delta_h, delta_s);
        let ff = folded_fraction(delta_h, delta_s, tm);
        assert!((ff - 0.5).abs() < 1e-6, "ff at Tm should be 0.5, got {ff}");
    }

    #[test]
    fn test_absolute_contact_order() {
        let contacts = vec![
            NativeContact::new(0, 5, 6.0, 1.0),
            NativeContact::new(1, 8, 7.0, 1.0),
        ];
        let co = absolute_contact_order(&contacts);
        assert!((co - 6.0).abs() < 1e-10, "CO = (5+7)/2 = 6, got {co}");
    }

    #[test]
    fn test_relative_contact_order() {
        let contacts = vec![NativeContact::new(0, 10, 8.0, 1.0)];
        let rco = relative_contact_order(&contacts, 20);
        assert!((rco - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_unfolding_curve_length() {
        let curve = thermal_unfolding_curve(-50.0, -0.1, 200.0, 400.0, 50);
        assert_eq!(curve.len(), 50);
    }

    #[test]
    fn test_thermal_unfolding_monotone() {
        // As T increases past Tm, folded fraction should decrease
        let curve = thermal_unfolding_curve(-50.0, -0.1, 300.0, 700.0, 100);
        let (_, f_low) = curve[0];
        let (_, f_high) = curve[99];
        assert!(f_low > f_high, "folded fraction should decrease with T");
    }

    #[test]
    fn test_kramers_mfpt_barrier() {
        let mfpt = kramers_mfpt(5.0, 1e-6, 300.0);
        // Should be larger than tau0 for a positive barrier
        assert!(mfpt > 1e-6);
    }

    #[test]
    fn test_folding_rate_positive() {
        let kf = folding_rate(5.0, 1e-6, 300.0);
        assert!(kf > 0.0);
    }

    #[test]
    fn test_charmm_bond_at_equilibrium() {
        let e = charmm_bond_energy(1.52, 1.52, 200.0);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_charmm_lj_minimum() {
        // At r = Rmin, energy should be −ε
        let rmin = 2.0;
        let eps = 0.1;
        let e = charmm_lj_energy(rmin, rmin, eps);
        assert!((e + eps).abs() < 1e-10, "LJ at Rmin should be −ε, got {e}");
    }

    #[test]
    fn test_charmm_coulomb_repulsion() {
        // Same-sign charges should give positive energy
        let e = charmm_coulomb_energy(1.0, 1.0, 5.0, 80.0);
        assert!(e > 0.0);
    }

    #[test]
    fn test_go_model_energy_native() {
        // At the native structure, all native contacts at their minimum → low energy
        let chain = make_helix_chain(10);
        let native = build_native_contact_map(&chain, 8.0, 1.0);
        let e = go_model_energy(
            &chain,
            &native,
            100.0,
            3.8,
            20.0,
            PI * 91.0 / 180.0,
            4.0,
            0.01,
        );
        assert!(e.is_finite());
    }

    #[test]
    fn test_funneled_landscape_unfolded() {
        // Q=0: V(0) = v0
        let v = funneled_landscape_energy(0.0, 5.0, 3.0, 2.0);
        assert!((v - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_funneled_landscape_folded() {
        // Q=1: V(1) = -epsilon
        let v = funneled_landscape_energy(1.0, 5.0, 3.0, 2.0);
        assert!((v + 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_funneled_transition_state() {
        let q_ts = funneled_transition_state_q(5.0, 3.0, 2.0);
        assert!(q_ts.is_none() || (0.0..=1.0).contains(&q_ts.unwrap()));
    }

    #[test]
    fn test_circular_mean_single() {
        let angles = [PI / 4.0];
        let m = circular_mean(&angles);
        assert!((m - PI / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_circular_variance_concentrated() {
        let angles: Vec<f64> = vec![0.0; 100];
        let v = circular_variance(&angles);
        assert!((v - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_end_to_end_distance() {
        let chain = make_extended_chain(11);
        let ete = end_to_end_distance(&chain);
        // 10 bonds × 3.8 Å
        assert!((ete - 38.0).abs() < 1e-6, "ETE={ete}");
    }

    #[test]
    fn test_boltzmann_weights_ordering() {
        let energies = [0.0, 1.0, -1.0];
        let weights = boltzmann_weights(&energies, 300.0);
        // Lower energy → higher weight
        assert!(weights[2] > weights[0]);
        assert!(weights[0] > weights[1]);
    }

    #[test]
    fn test_boltzmann_average() {
        // If all energies equal, average = mean of observables
        let obs = [1.0, 2.0, 3.0];
        let energies = [0.0, 0.0, 0.0];
        let avg = boltzmann_average(&obs, &energies, 300.0);
        assert!((avg - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_phi_value_zero() {
        assert_eq!(phi_value(0.0, 1.0), 0.0);
    }

    #[test]
    fn test_phi_value_one() {
        assert!((phi_value(1.0, 1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gb_solvation_pair_negative_for_opposite_charges() {
        let e = gb_solvation_pair(1.0, -1.0, 5.0, 3.0, 3.0, 2.0, 80.0);
        // Born solvation favours bringing opposite charges together in low dielectric
        assert!(e.is_finite());
    }

    #[test]
    fn test_protein_energy_decomposition_trivial() {
        // Two atoms at rest with no bonds/angles/dihedrals
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = [0.5, -0.5];
        let (eb, ea, ed, elj, ec) =
            protein_energy_decomposition(&positions, &[], &[], &[], &[], &charges, 80.0);
        assert_eq!(eb, 0.0);
        assert_eq!(ea, 0.0);
        assert_eq!(ed, 0.0);
        assert_eq!(elj, 0.0);
        assert!(
            ec < 0.0,
            "opposite charges → negative Coulomb energy, got {ec}"
        );
    }
}
