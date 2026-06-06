// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! DNA/RNA molecular dynamics using the oxDNA coarse-grained model.
//!
//! This module provides:
//!
//! - **Base-pair interactions** (oxDNA model): hydrogen bonding, stacking,
//!   cross-stacking, and coaxial stacking potentials.
//! - **Backbone connectivity**: FENE spring model for the sugar-phosphate backbone.
//! - **Melting temperature prediction**: nearest-neighbour thermodynamic model
//!   (SantaLucia 1998) for duplex stability.
//! - **Persistence length**: worm-like chain (WLC) tangent-correlation estimator.
//! - **Supercoiling / linking number**: White's theorem decomposition into twist
//!   and writhe (Fuller's integral).
//! - **Denaturation bubbles**: detection and characterization of local melting.
//! - **Nucleosome wrapping**: simple elastic-rod model of DNA wrapped around a
//!   histone core particle.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Nucleobase types
// ─────────────────────────────────────────────────────────────────────────────

/// Canonical nucleobases for DNA and RNA.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NucleoBase {
    /// Adenine (purine).
    Adenine,
    /// Thymine (pyrimidine, DNA only).
    Thymine,
    /// Guanine (purine).
    Guanine,
    /// Cytosine (pyrimidine).
    Cytosine,
    /// Uracil (pyrimidine, RNA only).
    Uracil,
}

impl NucleoBase {
    /// Watson-Crick complement (DNA rules; Uracil pairs with Adenine).
    pub fn complement(self) -> Self {
        match self {
            NucleoBase::Adenine => NucleoBase::Thymine,
            NucleoBase::Thymine => NucleoBase::Adenine,
            NucleoBase::Guanine => NucleoBase::Cytosine,
            NucleoBase::Cytosine => NucleoBase::Guanine,
            NucleoBase::Uracil => NucleoBase::Adenine,
        }
    }

    /// RNA complement (A-U, G-C).
    pub fn rna_complement(self) -> Self {
        match self {
            NucleoBase::Adenine => NucleoBase::Uracil,
            NucleoBase::Uracil => NucleoBase::Adenine,
            NucleoBase::Guanine => NucleoBase::Cytosine,
            NucleoBase::Cytosine => NucleoBase::Guanine,
            NucleoBase::Thymine => NucleoBase::Adenine,
        }
    }

    /// Single-letter IUPAC code.
    pub fn code(self) -> char {
        match self {
            NucleoBase::Adenine => 'A',
            NucleoBase::Thymine => 'T',
            NucleoBase::Guanine => 'G',
            NucleoBase::Cytosine => 'C',
            NucleoBase::Uracil => 'U',
        }
    }

    /// Parse from single-letter code (case-insensitive). Returns `None` for
    /// unrecognised characters.
    pub fn from_code(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(NucleoBase::Adenine),
            'T' => Some(NucleoBase::Thymine),
            'G' => Some(NucleoBase::Guanine),
            'C' => Some(NucleoBase::Cytosine),
            'U' => Some(NucleoBase::Uracil),
            _ => None,
        }
    }

    /// Whether this base is a purine (A or G).
    pub fn is_purine(self) -> bool {
        matches!(self, NucleoBase::Adenine | NucleoBase::Guanine)
    }

    /// Whether this base is a pyrimidine (C, T or U).
    pub fn is_pyrimidine(self) -> bool {
        !self.is_purine()
    }
}

/// Nucleic acid type selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleicAcidType {
    /// Double-stranded DNA.
    DNA,
    /// Single-stranded RNA.
    RNA,
}

// ─────────────────────────────────────────────────────────────────────────────
// oxDNA nucleotide representation
// ─────────────────────────────────────────────────────────────────────────────

/// A coarse-grained nucleotide in the oxDNA model.
///
/// Each nucleotide has a backbone site, a stacking site, and a hydrogen-bonding
/// site.  Positions are stored as `[f64; 3]` arrays.
#[derive(Debug, Clone)]
pub struct OxDnaNucleotide {
    /// Position of the backbone bead.
    pub backbone: [f64; 3],
    /// Unit vector from backbone toward the base (base-normal).
    pub base_normal: [f64; 3],
    /// Unit vector along the backbone tangent (5' -> 3' direction).
    pub tangent: [f64; 3],
    /// Which base this nucleotide carries.
    pub base: NucleoBase,
    /// Index of the strand this nucleotide belongs to.
    pub strand_id: usize,
    /// Index within the strand (0-based, 5' end = 0).
    pub index_in_strand: usize,
    /// Velocity of the backbone bead.
    pub velocity: [f64; 3],
    /// Angular velocity.
    pub angular_velocity: [f64; 3],
}

impl OxDnaNucleotide {
    /// Compute the position of the hydrogen-bonding site.
    ///
    /// In the oxDNA model, the H-bond site is offset from the backbone along the
    /// base-normal direction by `hb_offset` (default ~0.4 in reduced units).
    pub fn hbond_site(&self, hb_offset: f64) -> [f64; 3] {
        [
            self.backbone[0] + self.base_normal[0] * hb_offset,
            self.backbone[1] + self.base_normal[1] * hb_offset,
            self.backbone[2] + self.base_normal[2] * hb_offset,
        ]
    }

    /// Compute the stacking interaction site.
    ///
    /// The stacking site sits between the backbone and the H-bond site, offset by
    /// `stack_offset` along the base normal.
    pub fn stacking_site(&self, stack_offset: f64) -> [f64; 3] {
        [
            self.backbone[0] + self.base_normal[0] * stack_offset,
            self.backbone[1] + self.base_normal[1] * stack_offset,
            self.backbone[2] + self.base_normal[2] * stack_offset,
        ]
    }

    /// Create a default nucleotide at the given backbone position.
    pub fn new(backbone: [f64; 3], base: NucleoBase, strand_id: usize, index: usize) -> Self {
        Self {
            backbone,
            base_normal: [1.0, 0.0, 0.0],
            tangent: [0.0, 0.0, 1.0],
            base,
            strand_id,
            index_in_strand: index,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// oxDNA model parameters
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for the oxDNA2 coarse-grained model.
///
/// All quantities are in oxDNA reduced units where the energy unit is
/// approximately `k_B T` at 300 K and the length unit is ~8.518 Angstrom.
#[derive(Debug, Clone)]
pub struct OxDnaParams {
    /// FENE spring constant for backbone connectivity (reduced units).
    pub fene_k: f64,
    /// FENE maximum extension (reduced units).
    pub fene_r0: f64,
    /// FENE equilibrium distance.
    pub fene_delta: f64,
    /// Hydrogen bond well depth for AT pairs (reduced energy).
    pub hbond_at_eps: f64,
    /// Hydrogen bond well depth for GC pairs (reduced energy).
    pub hbond_gc_eps: f64,
    /// Stacking interaction strength.
    pub stack_eps: f64,
    /// Cross-stacking interaction strength.
    pub cross_stack_eps: f64,
    /// Coaxial stacking interaction strength.
    pub coaxial_stack_eps: f64,
    /// Excluded volume repulsion strength.
    pub exc_vol_eps: f64,
    /// Excluded volume sigma parameter.
    pub exc_vol_sigma: f64,
    /// Temperature in reduced units (k_B T).
    pub temperature: f64,
    /// H-bond site offset from backbone along base normal.
    pub hb_offset: f64,
    /// Stacking site offset from backbone along base normal.
    pub stack_offset: f64,
}

impl Default for OxDnaParams {
    fn default() -> Self {
        Self {
            fene_k: 30.0,
            fene_r0: 0.7525,
            fene_delta: 0.25,
            hbond_at_eps: 1.077,
            hbond_gc_eps: 1.542,
            stack_eps: 1.3448,
            cross_stack_eps: 0.4,
            coaxial_stack_eps: 1.3448,
            exc_vol_eps: 2.0,
            exc_vol_sigma: 0.35,
            temperature: 1.0,
            hb_offset: 0.4,
            stack_offset: 0.34,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FENE backbone potential
// ─────────────────────────────────────────────────────────────────────────────

/// FENE (finitely extensible nonlinear elastic) spring for backbone connectivity.
///
/// V_FENE = -0.5 * k * delta^2 * ln(1 - ((r - r0) / delta)^2)
///
/// Diverges when |r - r0| >= delta, confining the bond length.
pub fn fene_potential(r: f64, params: &OxDnaParams) -> f64 {
    let dr = r - params.fene_r0;
    let ratio = dr / params.fene_delta;
    let ratio_sq = ratio * ratio;
    if ratio_sq >= 1.0 {
        return 1e6; // effectively infinite
    }
    -0.5 * params.fene_k * params.fene_delta * params.fene_delta * (1.0 - ratio_sq).ln()
}

/// Force magnitude from FENE potential (dV/dr, positive = repulsive).
pub fn fene_force(r: f64, params: &OxDnaParams) -> f64 {
    let dr = r - params.fene_r0;
    let ratio_sq = (dr / params.fene_delta).powi(2);
    if ratio_sq >= 1.0 {
        return 1e6;
    }
    params.fene_k * dr / (1.0 - ratio_sq)
}

// ─────────────────────────────────────────────────────────────────────────────
// Hydrogen bonding potential
// ─────────────────────────────────────────────────────────────────────────────

/// Hydrogen-bond potential between two complementary bases.
///
/// Uses a Morse-like form:  V_hb = eps * f(r) * g(theta)
/// where f(r) is a smooth distance cutoff and g(theta) penalises misalignment.
pub fn hydrogen_bond_energy(
    site_i: [f64; 3],
    site_j: [f64; 3],
    normal_i: [f64; 3],
    normal_j: [f64; 3],
    base_i: NucleoBase,
    base_j: NucleoBase,
    params: &OxDnaParams,
) -> f64 {
    // Only complementary bases form hydrogen bonds
    if base_i.complement() != base_j {
        return 0.0;
    }

    let r = dist3(site_i, site_j);
    let r_eq = 0.2; // equilibrium H-bond distance (reduced units)
    let alpha = 8.0;

    // Morse-like radial factor
    let exp_term = (-alpha * (r - r_eq)).exp();
    let f_r = exp_term * (exp_term - 2.0);

    // Angular modulation: bases should be roughly antiparallel
    let cos_theta =
        -(normal_i[0] * normal_j[0] + normal_i[1] * normal_j[1] + normal_i[2] * normal_j[2]);
    let g_theta = if cos_theta > 0.0 {
        cos_theta.powi(4)
    } else {
        0.0
    };

    let eps = match (base_i, base_j) {
        (NucleoBase::Guanine, NucleoBase::Cytosine)
        | (NucleoBase::Cytosine, NucleoBase::Guanine) => params.hbond_gc_eps,
        _ => params.hbond_at_eps,
    };

    eps * f_r * g_theta
}

// ─────────────────────────────────────────────────────────────────────────────
// Stacking interaction
// ─────────────────────────────────────────────────────────────────────────────

/// Stacking interaction energy between consecutive bases on the same strand.
///
/// V_stack = eps_stack * f(r) * g(theta)
pub fn stacking_energy(
    site_i: [f64; 3],
    site_j: [f64; 3],
    tangent_i: [f64; 3],
    tangent_j: [f64; 3],
    params: &OxDnaParams,
) -> f64 {
    let r = dist3(site_i, site_j);
    let r_eq = 0.34; // ~3.4 Angstrom / 8.518
    let alpha = 6.0;

    // Morse radial part
    let exp_term = (-alpha * (r - r_eq)).exp();
    let f_r = exp_term * (exp_term - 2.0);

    // Angular modulation: tangent alignment
    let cos_t =
        tangent_i[0] * tangent_j[0] + tangent_i[1] * tangent_j[1] + tangent_i[2] * tangent_j[2];
    let g_t = if cos_t > 0.0 { cos_t.powi(2) } else { 0.0 };

    params.stack_eps * f_r * g_t
}

/// Cross-stacking energy between bases on opposite strands that are diagonally
/// adjacent (i.e., base i on strand 1 with base j+1 on strand 2).
pub fn cross_stacking_energy(
    site_i: [f64; 3],
    site_j: [f64; 3],
    normal_i: [f64; 3],
    normal_j: [f64; 3],
    params: &OxDnaParams,
) -> f64 {
    let r = dist3(site_i, site_j);
    let r_eq = 0.44;
    let alpha = 4.0;

    let exp_term = (-alpha * (r - r_eq)).exp();
    let f_r = exp_term * (exp_term - 2.0);

    let cos_theta =
        -(normal_i[0] * normal_j[0] + normal_i[1] * normal_j[1] + normal_i[2] * normal_j[2]);
    let g_theta = if cos_theta > 0.0 {
        cos_theta.powi(2)
    } else {
        0.0
    };

    params.cross_stack_eps * f_r * g_theta
}

/// Coaxial stacking energy for blunt-end stacking of two duplexes.
pub fn coaxial_stacking_energy(
    site_i: [f64; 3],
    site_j: [f64; 3],
    tangent_i: [f64; 3],
    tangent_j: [f64; 3],
    params: &OxDnaParams,
) -> f64 {
    let r = dist3(site_i, site_j);
    let r_eq = 0.34;
    let alpha = 6.0;

    let exp_term = (-alpha * (r - r_eq)).exp();
    let f_r = exp_term * (exp_term - 2.0);

    let cos_t =
        tangent_i[0] * tangent_j[0] + tangent_i[1] * tangent_j[1] + tangent_i[2] * tangent_j[2];
    let g_t = if cos_t > 0.0 { cos_t.powi(4) } else { 0.0 };

    params.coaxial_stack_eps * f_r * g_t
}

// ─────────────────────────────────────────────────────────────────────────────
// Excluded volume
// ─────────────────────────────────────────────────────────────────────────────

/// Soft-core excluded volume repulsion (truncated and shifted LJ-like).
///
/// V_exc = eps * ((sigma / r)^12 - 2*(sigma / r)^6 + 1)  for r < sigma
///       = 0                                                for r >= sigma
pub fn excluded_volume_energy(r: f64, params: &OxDnaParams) -> f64 {
    if r >= params.exc_vol_sigma {
        return 0.0;
    }
    let inv = params.exc_vol_sigma / r.max(1e-14);
    let inv6 = inv.powi(6);
    let inv12 = inv6 * inv6;
    params.exc_vol_eps * (inv12 - 2.0 * inv6 + 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Full oxDNA strand representation
// ─────────────────────────────────────────────────────────────────────────────

/// A strand of nucleotides in the oxDNA model.
#[derive(Debug, Clone)]
pub struct OxDnaStrand {
    /// Nucleotides ordered 5' -> 3'.
    pub nucleotides: Vec<OxDnaNucleotide>,
    /// Unique strand identifier.
    pub strand_id: usize,
    /// Whether this strand is circular (e.g., plasmid DNA).
    pub is_circular: bool,
}

impl OxDnaStrand {
    /// Create a new strand from a base sequence string (e.g., "ATCG").
    ///
    /// Nucleotide positions are arranged along a B-form helix with 10.5 bp/turn
    /// and 0.34 nm rise per bp.
    pub fn from_sequence(seq: &str, strand_id: usize) -> Self {
        let rise = 0.34;
        let twist_per_bp = 2.0 * PI / 10.5;
        let helix_radius = 1.0; // nm

        let mut nucleotides = Vec::with_capacity(seq.len());
        for (i, c) in seq.chars().enumerate() {
            let base = NucleoBase::from_code(c).unwrap_or(NucleoBase::Adenine);
            let z = i as f64 * rise;
            let angle = i as f64 * twist_per_bp;
            let x = helix_radius * angle.cos();
            let y = helix_radius * angle.sin();

            let tangent = normalize3([0.0, 0.0, 1.0]);
            let base_normal = normalize3([angle.cos(), angle.sin(), 0.0]);

            nucleotides.push(OxDnaNucleotide {
                backbone: [x, y, z],
                base_normal,
                tangent,
                base,
                strand_id,
                index_in_strand: i,
                velocity: [0.0; 3],
                angular_velocity: [0.0; 3],
            });
        }

        Self {
            nucleotides,
            strand_id,
            is_circular: false,
        }
    }

    /// Return the sequence as a string of one-letter codes.
    pub fn sequence_string(&self) -> String {
        self.nucleotides.iter().map(|n| n.base.code()).collect()
    }

    /// Number of nucleotides in the strand.
    pub fn len(&self) -> usize {
        self.nucleotides.len()
    }

    /// Whether the strand is empty.
    pub fn is_empty(&self) -> bool {
        self.nucleotides.is_empty()
    }

    /// End-to-end distance of the backbone.
    pub fn end_to_end_distance(&self) -> f64 {
        if self.nucleotides.len() < 2 {
            return 0.0;
        }
        let first = &self.nucleotides[0].backbone;
        let last = &self.nucleotides[self.nucleotides.len() - 1].backbone;
        dist3(*first, *last)
    }

    /// Contour length (sum of backbone bond lengths).
    pub fn contour_length(&self) -> f64 {
        let mut len = 0.0;
        for w in self.nucleotides.windows(2) {
            len += dist3(w[0].backbone, w[1].backbone);
        }
        if self.is_circular && self.nucleotides.len() > 1 {
            len += dist3(
                self.nucleotides[self.nucleotides.len() - 1].backbone,
                self.nucleotides[0].backbone,
            );
        }
        len
    }

    /// GC content as a fraction in \[0, 1\].
    pub fn gc_content(&self) -> f64 {
        if self.nucleotides.is_empty() {
            return 0.0;
        }
        let gc = self
            .nucleotides
            .iter()
            .filter(|n| matches!(n.base, NucleoBase::Guanine | NucleoBase::Cytosine))
            .count();
        gc as f64 / self.nucleotides.len() as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Duplex system
// ─────────────────────────────────────────────────────────────────────────────

/// A double-stranded DNA duplex in the oxDNA model.
///
/// Consists of two antiparallel strands with Watson-Crick base pairing.
#[derive(Debug, Clone)]
pub struct OxDnaDuplex {
    /// The sense strand (5' -> 3').
    pub sense: OxDnaStrand,
    /// The antisense strand (3' -> 5', stored 5' -> 3' of antisense).
    pub antisense: OxDnaStrand,
    /// Model parameters.
    pub params: OxDnaParams,
}

impl OxDnaDuplex {
    /// Build a duplex from a sense-strand sequence.  The antisense strand is
    /// automatically generated as the reverse complement.
    pub fn from_sense_sequence(seq: &str) -> Self {
        let sense = OxDnaStrand::from_sequence(seq, 0);
        let antisense_seq: String = seq
            .chars()
            .rev()
            .map(|c| {
                NucleoBase::from_code(c)
                    .unwrap_or(NucleoBase::Adenine)
                    .complement()
                    .code()
            })
            .collect();
        let antisense = OxDnaStrand::from_sequence(&antisense_seq, 1);
        Self {
            sense,
            antisense,
            params: OxDnaParams::default(),
        }
    }

    /// Number of base pairs.
    pub fn num_base_pairs(&self) -> usize {
        self.sense.len().min(self.antisense.len())
    }

    /// Total backbone FENE energy for both strands.
    pub fn backbone_energy(&self) -> f64 {
        let mut e = 0.0;
        for strand in [&self.sense, &self.antisense] {
            for w in strand.nucleotides.windows(2) {
                let r = dist3(w[0].backbone, w[1].backbone);
                e += fene_potential(r, &self.params);
            }
        }
        e
    }

    /// Total hydrogen-bonding energy across all base pairs.
    pub fn hbond_energy(&self) -> f64 {
        let n = self.num_base_pairs();
        let mut e = 0.0;
        for i in 0..n {
            let ni = &self.sense.nucleotides[i];
            // antisense is reversed, so pair i with (n-1-i)
            let j = n - 1 - i;
            if j >= self.antisense.nucleotides.len() {
                continue;
            }
            let nj = &self.antisense.nucleotides[j];
            let si = ni.hbond_site(self.params.hb_offset);
            let sj = nj.hbond_site(self.params.hb_offset);
            e += hydrogen_bond_energy(
                si,
                sj,
                ni.base_normal,
                nj.base_normal,
                ni.base,
                nj.base,
                &self.params,
            );
        }
        e
    }

    /// Total stacking energy within each strand.
    pub fn stacking_energy_total(&self) -> f64 {
        let mut e = 0.0;
        for strand in [&self.sense, &self.antisense] {
            for w in strand.nucleotides.windows(2) {
                let si = w[0].stacking_site(self.params.stack_offset);
                let sj = w[1].stacking_site(self.params.stack_offset);
                e += stacking_energy(si, sj, w[0].tangent, w[1].tangent, &self.params);
            }
        }
        e
    }

    /// Total potential energy (backbone + H-bond + stacking).
    pub fn total_potential_energy(&self) -> f64 {
        self.backbone_energy() + self.hbond_energy() + self.stacking_energy_total()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Melting temperature prediction (nearest-neighbour model)
// ─────────────────────────────────────────────────────────────────────────────

/// Nearest-neighbour enthalpy (kcal/mol) and entropy (cal/mol/K) parameters
/// for DNA duplex formation (SantaLucia 1998, unified parameters).
///
/// Index: dinucleotide pair 5'-XY-3' / 3'-X'Y'-5'.
#[derive(Debug, Clone, Copy)]
struct NnParams {
    /// Enthalpy change in kcal/mol.
    dh: f64,
    /// Entropy change in cal/(mol K).
    ds: f64,
}

/// Look up the nearest-neighbour parameters for a dinucleotide step.
///
/// Returns (delta_H in kcal/mol, delta_S in cal/mol/K).
fn nn_lookup(a: NucleoBase, b: NucleoBase) -> NnParams {
    use NucleoBase::*;
    match (a, b) {
        (Adenine, Adenine) | (Thymine, Thymine) => NnParams {
            dh: -7.9,
            ds: -22.2,
        },
        (Adenine, Thymine) => NnParams {
            dh: -7.2,
            ds: -20.4,
        },
        (Thymine, Adenine) => NnParams {
            dh: -7.2,
            ds: -21.3,
        },
        (Adenine, Guanine) | (Cytosine, Thymine) => NnParams {
            dh: -7.8,
            ds: -21.0,
        },
        (Adenine, Cytosine) | (Guanine, Thymine) => NnParams {
            dh: -8.4,
            ds: -22.4,
        },
        (Guanine, Adenine) | (Thymine, Cytosine) => NnParams {
            dh: -8.2,
            ds: -22.2,
        },
        (Cytosine, Adenine) | (Thymine, Guanine) => NnParams {
            dh: -8.5,
            ds: -22.7,
        },
        (Guanine, Guanine) | (Cytosine, Cytosine) => NnParams {
            dh: -8.0,
            ds: -19.9,
        },
        (Guanine, Cytosine) => NnParams {
            dh: -9.8,
            ds: -24.4,
        },
        (Cytosine, Guanine) => NnParams {
            dh: -10.6,
            ds: -27.2,
        },
        _ => NnParams {
            dh: -7.5,
            ds: -21.0,
        }, // fallback for RNA bases
    }
}

/// Predict the melting temperature (in Kelvin) of a DNA duplex using the
/// nearest-neighbour model.
///
/// `seq` is the 5'->3' sense strand sequence.
/// `c_total` is the total strand concentration in mol/L (typically 1e-4 M).
/// `na_conc` is the Na+ concentration in mol/L (typically 1.0 M for standard).
///
/// Uses the SantaLucia (1998) unified nearest-neighbour parameters.
pub fn melting_temperature(seq: &str, c_total: f64, na_conc: f64) -> f64 {
    let bases: Vec<NucleoBase> = seq.chars().filter_map(NucleoBase::from_code).collect();

    if bases.len() < 2 {
        return 0.0;
    }

    // Sum nearest-neighbour contributions
    let mut dh_total = 0.0_f64;
    let mut ds_total = 0.0_f64;

    // Initiation parameters
    let init_dh = 0.2; // kcal/mol (averaged AT/GC init)
    let init_ds = -5.7; // cal/mol/K
    dh_total += 2.0 * init_dh;
    ds_total += 2.0 * init_ds;

    for w in bases.windows(2) {
        let nn = nn_lookup(w[0], w[1]);
        dh_total += nn.dh;
        ds_total += nn.ds;
    }

    // Convert dH to cal/mol to match dS units
    let dh_cal = dh_total * 1000.0;

    // Salt correction (Owczarzy et al. simplified)
    let salt_correction = 12.5 * na_conc.log10();
    let ds_corrected = ds_total + salt_correction;

    // Non-self-complementary: Tm = dH / (dS + R*ln(Ct/4))
    let r = 1.987; // cal/(mol K), gas constant
    let ct = c_total.max(1e-12);
    let tm = dh_cal / (ds_corrected + r * (ct / 4.0).ln());

    tm.max(0.0)
}

/// Convert Kelvin to Celsius.
pub fn kelvin_to_celsius(k: f64) -> f64 {
    k - 273.15
}

// ─────────────────────────────────────────────────────────────────────────────
// Persistence length estimation
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the persistence length from a set of backbone tangent vectors
/// along a polymer chain.
///
/// Uses the tangent-tangent autocorrelation:
///   <t(0) . t(s)> = exp(-s / l_p)
///
/// Fits a simple exponential decay to the correlation function and returns
/// the persistence length `l_p` in units of the bond length.
pub fn persistence_length(tangents: &[[f64; 3]], bond_length: f64) -> f64 {
    let n = tangents.len();
    if n < 3 {
        return 0.0;
    }

    // Compute tangent-tangent correlation for separations s = 1..n/2
    let max_s = n / 2;
    let mut sum_log_corr = 0.0;
    let mut sum_s = 0.0;
    let mut count = 0;

    for s in 1..=max_s {
        let mut corr_sum = 0.0;
        let mut corr_count = 0;
        for i in 0..n - s {
            let dot = tangents[i][0] * tangents[i + s][0]
                + tangents[i][1] * tangents[i + s][1]
                + tangents[i][2] * tangents[i + s][2];
            corr_sum += dot;
            corr_count += 1;
        }
        if corr_count > 0 {
            let avg_corr = corr_sum / corr_count as f64;
            if avg_corr > 0.01 {
                sum_log_corr += avg_corr.ln();
                sum_s += s as f64;
                count += 1;
            }
        }
    }

    if count == 0 {
        // All correlations were positive and above 0.01; if the log sum is
        // essentially zero, the chain is straight => effectively infinite lp.
        return if n > 2 { 1e6 * bond_length } else { 0.0 };
    }

    if sum_log_corr.abs() < 1e-12 {
        // Correlations are all ~1.0 (perfectly straight chain)
        return 1e6 * bond_length;
    }

    // <cos(theta)> ~ exp(-s / lp)  =>  ln(<cos(theta)>) = -s / lp
    // Simple average: lp = -<s> / <ln(corr)>
    let avg_s = sum_s / count as f64;
    let avg_log_corr = sum_log_corr / count as f64;
    let lp = -avg_s * bond_length / avg_log_corr;
    lp.max(0.0)
}

/// Compute end-to-end distance prediction from the worm-like chain model.
///
/// <R^2> = 2 * L_p * L * (1 - L_p/L * (1 - exp(-L/L_p)))
///
/// where L is the contour length and L_p is the persistence length.
pub fn wlc_end_to_end_rms(contour_length: f64, persistence_length: f64) -> f64 {
    if contour_length <= 0.0 || persistence_length <= 0.0 {
        return 0.0;
    }
    let lp = persistence_length;
    let ll = contour_length;
    let r2 = 2.0 * lp * ll * (1.0 - lp / ll * (1.0 - (-ll / lp).exp()));
    r2.max(0.0).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Supercoiling and linking number
// ─────────────────────────────────────────────────────────────────────────────

/// White's theorem decomposition: Lk = Tw + Wr.
///
/// Contains the linking number, twist, and writhe of a closed DNA molecule.
#[derive(Debug, Clone, Copy)]
pub struct LinkingDecomposition {
    /// Linking number (integer for a closed curve).
    pub linking_number: f64,
    /// Twist: integral of the local torsion along the backbone.
    pub twist: f64,
    /// Writhe: a measure of the coiling of the axis in 3-D space.
    pub writhe: f64,
}

/// Compute the twist of a closed curve given backbone positions and a set of
/// reference normal vectors (base normals).
///
/// Twist = (1 / 2 pi) * sum of angle changes of the normal vector around the
/// tangent as we traverse the backbone.
pub fn compute_twist(backbone: &[[f64; 3]], normals: &[[f64; 3]]) -> f64 {
    let n = backbone.len();
    if n < 3 {
        return 0.0;
    }

    let mut twist = 0.0;
    for i in 0..n {
        let i_next = (i + 1) % n;
        let tangent = normalize3(sub3(backbone[i_next], backbone[i]));

        // Project normals onto plane perpendicular to tangent
        let ni = project_out(normals[i], tangent);
        let nj = project_out(normals[i_next], tangent);
        let ni = normalize3(ni);
        let nj = normalize3(nj);

        let dot = (ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2]).clamp(-1.0, 1.0);
        let cross = cross3(ni, nj);
        let sin_sign = cross[0] * tangent[0] + cross[1] * tangent[1] + cross[2] * tangent[2];
        let angle = dot.acos().copysign(sin_sign);
        twist += angle;
    }

    twist / (2.0 * PI)
}

/// Compute the writhe of a closed space curve using the discretised
/// Gauss linking integral (Fuller's formula for self-linking writhe).
///
/// Wr = (1 / 4 pi) * sum_{i != j} Omega_{ij}
///
/// where Omega_{ij} is the signed solid angle subtended by segments i and j.
pub fn compute_writhe(backbone: &[[f64; 3]]) -> f64 {
    let n = backbone.len();
    if n < 4 {
        return 0.0;
    }

    let mut wr = 0.0;
    for i in 0..n {
        let i1 = (i + 1) % n;
        for j in (i + 2)..n {
            if j == i || (j + 1) % n == i {
                continue;
            }
            let j1 = (j + 1) % n;
            wr += solid_angle_segment_pair(backbone[i], backbone[i1], backbone[j], backbone[j1]);
        }
    }

    wr / (4.0 * PI)
}

/// Compute linking number decomposition for a closed DNA molecule.
pub fn linking_decomposition(backbone: &[[f64; 3]], normals: &[[f64; 3]]) -> LinkingDecomposition {
    let twist = compute_twist(backbone, normals);
    let writhe = compute_writhe(backbone);
    LinkingDecomposition {
        linking_number: twist + writhe,
        twist,
        writhe,
    }
}

/// Compute the superhelical density sigma = (Lk - Lk0) / Lk0.
///
/// `lk` is the actual linking number, `lk0` is the relaxed linking number
/// (typically N_bp / 10.5 for B-DNA).
pub fn superhelical_density(lk: f64, lk0: f64) -> f64 {
    if lk0.abs() < 1e-15 {
        return 0.0;
    }
    (lk - lk0) / lk0
}

// ─────────────────────────────────────────────────────────────────────────────
// Denaturation bubbles
// ─────────────────────────────────────────────────────────────────────────────

/// A denaturation bubble: a contiguous region where base pairs are broken.
#[derive(Debug, Clone)]
pub struct DenaturationBubble {
    /// Index of the first base pair in the bubble (0-based).
    pub start: usize,
    /// Number of base pairs in the bubble.
    pub length: usize,
    /// Average opening distance (backbone separation) in the bubble region.
    pub avg_opening: f64,
}

/// Detect denaturation bubbles in a duplex.
///
/// A base pair is considered "open" when the hydrogen-bond energy is above
/// `threshold` (closer to zero = weaker).  Contiguous open regions form bubbles.
pub fn detect_denaturation_bubbles(
    duplex: &OxDnaDuplex,
    threshold: f64,
) -> Vec<DenaturationBubble> {
    let n = duplex.num_base_pairs();
    if n == 0 {
        return vec![];
    }

    // Evaluate per-bp H-bond energy
    let mut hb_energies = Vec::with_capacity(n);
    for i in 0..n {
        let ni = &duplex.sense.nucleotides[i];
        let j = n - 1 - i;
        if j >= duplex.antisense.nucleotides.len() {
            hb_energies.push(0.0);
            continue;
        }
        let nj = &duplex.antisense.nucleotides[j];
        let si = ni.hbond_site(duplex.params.hb_offset);
        let sj = nj.hbond_site(duplex.params.hb_offset);
        let e = hydrogen_bond_energy(
            si,
            sj,
            ni.base_normal,
            nj.base_normal,
            ni.base,
            nj.base,
            &duplex.params,
        );
        hb_energies.push(e);
    }

    // Find contiguous regions where |hb_energy| < threshold (open)
    let mut bubbles = Vec::new();
    let mut i = 0;
    while i < n {
        if hb_energies[i].abs() < threshold {
            let start = i;
            let mut opening_sum = 0.0;
            let mut count = 0;
            while i < n && hb_energies[i].abs() < threshold {
                // Compute backbone separation
                let j = n - 1 - i;
                if j < duplex.antisense.nucleotides.len() {
                    opening_sum += dist3(
                        duplex.sense.nucleotides[i].backbone,
                        duplex.antisense.nucleotides[j].backbone,
                    );
                }
                count += 1;
                i += 1;
            }
            bubbles.push(DenaturationBubble {
                start,
                length: count,
                avg_opening: if count > 0 {
                    opening_sum / count as f64
                } else {
                    0.0
                },
            });
        } else {
            i += 1;
        }
    }

    bubbles
}

/// Fraction of base pairs that are denatured (open).
pub fn denaturation_fraction(duplex: &OxDnaDuplex, threshold: f64) -> f64 {
    let bubbles = detect_denaturation_bubbles(duplex, threshold);
    let open: usize = bubbles.iter().map(|b| b.length).sum();
    let n = duplex.num_base_pairs();
    if n == 0 {
        return 0.0;
    }
    open as f64 / n as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// Nucleosome wrapping model
// ─────────────────────────────────────────────────────────────────────────────

/// Simple elastic-rod model parameters for nucleosome wrapping.
#[derive(Debug, Clone)]
pub struct NucleosomeParams {
    /// Histone core radius (nm).
    pub core_radius: f64,
    /// DNA bending stiffness A = l_p * k_B T (pN nm^2).
    pub bending_stiffness: f64,
    /// DNA twist stiffness C (pN nm^2).
    pub twist_stiffness: f64,
    /// Adsorption energy per unit length (pN).
    pub adsorption_energy: f64,
    /// Number of wrapping turns (typically ~1.65 for the nucleosome).
    pub wrapping_turns: f64,
}

impl Default for NucleosomeParams {
    fn default() -> Self {
        Self {
            core_radius: 4.18,
            bending_stiffness: 50.0 * 4.11, // l_p=50 nm, kBT=4.11 pN nm at 300K
            twist_stiffness: 75.0 * 4.11,
            adsorption_energy: 2.0, // pN (per nm of wrapped DNA)
            wrapping_turns: 1.65,
        }
    }
}

/// Compute the bending energy of DNA wrapped around a histone core.
///
/// E_bend = A / (2 R^2) * L
///
/// where A is the bending stiffness, R is the wrapping radius, and L is the
/// wrapped length.
pub fn nucleosome_bending_energy(params: &NucleosomeParams) -> f64 {
    let wrapped_length = 2.0 * PI * params.wrapping_turns * params.core_radius;
    let r = params.core_radius;
    params.bending_stiffness / (2.0 * r * r) * wrapped_length
}

/// Compute the adsorption energy gained by wrapping DNA around the histone.
///
/// E_ads = -epsilon * L  (negative = favourable)
pub fn nucleosome_adsorption_energy(params: &NucleosomeParams) -> f64 {
    let wrapped_length = 2.0 * PI * params.wrapping_turns * params.core_radius;
    -params.adsorption_energy * wrapped_length
}

/// Net wrapping free energy = bending + adsorption.
pub fn nucleosome_wrapping_energy(params: &NucleosomeParams) -> f64 {
    nucleosome_bending_energy(params) + nucleosome_adsorption_energy(params)
}

/// Compute the number of base pairs wrapped around the nucleosome.
///
/// Using the B-DNA rise of 0.34 nm/bp.
pub fn nucleosome_wrapped_bp(params: &NucleosomeParams) -> f64 {
    let wrapped_length = 2.0 * PI * params.wrapping_turns * params.core_radius;
    wrapped_length / 0.34
}

/// Equilibrium unwrapping angle at a given applied force `f` (pN).
///
/// A simple Marko-Siggia like estimate: at equilibrium, the torque from
/// bending equals the applied force times the lever arm.
///
/// theta_unwrap ~ sqrt(2 * f * R / A)
pub fn equilibrium_unwrap_angle(force: f64, params: &NucleosomeParams) -> f64 {
    if force <= 0.0 {
        return 0.0;
    }
    let r = params.core_radius;
    (2.0 * force * r / params.bending_stiffness).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Simulation step helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute all pairwise backbone FENE forces for a strand and update velocities.
///
/// Uses velocity-Verlet half-step: v += 0.5 * dt * F / m.
pub fn apply_backbone_forces(strand: &mut OxDnaStrand, params: &OxDnaParams, dt: f64, mass: f64) {
    let n = strand.nucleotides.len();
    if n < 2 {
        return;
    }

    // Collect forces first (to avoid borrow issues)
    let mut forces = vec![[0.0_f64; 3]; n];
    for i in 0..n - 1 {
        let pi = strand.nucleotides[i].backbone;
        let pj = strand.nucleotides[i + 1].backbone;
        let rij = sub3(pj, pi);
        let r = norm3(rij);
        if r < 1e-14 {
            continue;
        }
        let f_mag = fene_force(r, params);
        let f_vec = scale3(rij, f_mag / r);
        for d in 0..3 {
            forces[i][d] += f_vec[d];
            forces[i + 1][d] -= f_vec[d];
        }
    }

    // Apply half-step velocity update
    let half_dt_over_m = 0.5 * dt / mass.max(1e-15);
    for (nuc, f) in strand.nucleotides.iter_mut().zip(forces.iter()).take(n) {
        for (v, &fv) in nuc.velocity.iter_mut().zip(f.iter()) {
            *v += half_dt_over_m * fv;
        }
    }
}

/// Advance backbone positions by one Verlet step.
pub fn advance_positions(strand: &mut OxDnaStrand, dt: f64) {
    for nuc in &mut strand.nucleotides {
        for d in 0..3 {
            nuc.backbone[d] += nuc.velocity[d] * dt;
        }
    }
}

/// Compute the total kinetic energy of a strand.
pub fn kinetic_energy(strand: &OxDnaStrand, mass: f64) -> f64 {
    let mut ke = 0.0;
    for nuc in &strand.nucleotides {
        let v2 = nuc.velocity[0].powi(2) + nuc.velocity[1].powi(2) + nuc.velocity[2].powi(2);
        ke += 0.5 * mass * v2;
    }
    ke
}

/// Instantaneous temperature of the strand from kinetic energy.
///
/// T = 2 * KE / (3 * N * k_B), in reduced units where k_B = 1.
pub fn instantaneous_temperature(strand: &OxDnaStrand, mass: f64) -> f64 {
    let n = strand.nucleotides.len();
    if n == 0 {
        return 0.0;
    }
    let ke = kinetic_energy(strand, mass);
    2.0 * ke / (3.0 * n as f64)
}

// ─────────────────────────────────────────────────────────────────────────────
// Sequence analysis utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Count occurrences of each dinucleotide step in a sequence.
///
/// Returns a 16-element array indexed as 4*i + j where i, j are base indices
/// (A=0, T=1, G=2, C=3).
pub fn dinucleotide_counts(seq: &str) -> [usize; 16] {
    let mut counts = [0usize; 16];
    let bases: Vec<usize> = seq
        .chars()
        .filter_map(|c| match c.to_ascii_uppercase() {
            'A' => Some(0),
            'T' => Some(1),
            'G' => Some(2),
            'C' => Some(3),
            _ => None,
        })
        .collect();
    for w in bases.windows(2) {
        counts[4 * w[0] + w[1]] += 1;
    }
    counts
}

/// Compute the total nearest-neighbour enthalpy (kcal/mol) for a sequence.
pub fn nn_enthalpy(seq: &str) -> f64 {
    let bases: Vec<NucleoBase> = seq.chars().filter_map(NucleoBase::from_code).collect();
    let mut dh = 0.0;
    for w in bases.windows(2) {
        dh += nn_lookup(w[0], w[1]).dh;
    }
    dh
}

/// Compute the total nearest-neighbour entropy (cal/mol/K) for a sequence.
pub fn nn_entropy(seq: &str) -> f64 {
    let bases: Vec<NucleoBase> = seq.chars().filter_map(NucleoBase::from_code).collect();
    let mut ds = 0.0;
    for w in bases.windows(2) {
        ds += nn_lookup(w[0], w[1]).ds;
    }
    ds
}

/// Free energy of duplex formation at temperature `t_kelvin` using the NN model.
///
/// dG = dH - T * dS  (kcal/mol, with dS converted from cal/mol/K)
pub fn nn_free_energy(seq: &str, t_kelvin: f64) -> f64 {
    let dh = nn_enthalpy(seq);
    let ds = nn_entropy(seq) / 1000.0; // convert to kcal/mol/K
    dh - t_kelvin * ds
}

/// Check whether a sequence is self-complementary.
pub fn is_self_complementary(seq: &str) -> bool {
    let bases: Vec<NucleoBase> = seq.chars().filter_map(NucleoBase::from_code).collect();
    let n = bases.len();
    for i in 0..n {
        if bases[i].complement() != bases[n - 1 - i] {
            return false;
        }
    }
    true
}

/// Count the number of CpG dinucleotides in a sequence.
pub fn cpg_count(seq: &str) -> usize {
    let upper: String = seq.to_ascii_uppercase();
    upper
        .as_bytes()
        .windows(2)
        .filter(|w| w[0] == b'C' && w[1] == b'G')
        .count()
}

// ─────────────────────────────────────────────────────────────────────────────
// B-DNA geometry helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Generate ideal B-DNA backbone coordinates for `n` base pairs.
///
/// Returns two vectors of positions (sense and antisense strands).
pub fn ideal_bdna_backbone(n: usize) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let rise = 0.34; // nm per bp
    let twist_per_bp = 2.0 * PI / 10.5;
    let helix_r = 1.0; // nm

    let mut sense = Vec::with_capacity(n);
    let mut antisense = Vec::with_capacity(n);

    for i in 0..n {
        let z = i as f64 * rise;
        let theta = i as f64 * twist_per_bp;

        sense.push([helix_r * theta.cos(), helix_r * theta.sin(), z]);
        // Antisense is offset by pi
        antisense.push([
            helix_r * (theta + PI).cos(),
            helix_r * (theta + PI).sin(),
            z,
        ]);
    }

    (sense, antisense)
}

/// Compute the helical repeat (base pairs per full turn) from backbone
/// positions.
pub fn helical_repeat(backbone: &[[f64; 3]]) -> f64 {
    if backbone.len() < 3 {
        return 0.0;
    }

    // Measure angle change per step in the xy-plane
    let mut total_angle = 0.0;
    let mut count = 0;
    for w in backbone.windows(2) {
        let a1 = w[0][1].atan2(w[0][0]);
        let a2 = w[1][1].atan2(w[1][0]);
        let mut da = a2 - a1;
        if da > PI {
            da -= 2.0 * PI;
        }
        if da < -PI {
            da += 2.0 * PI;
        }
        total_angle += da.abs();
        count += 1;
    }

    if total_angle < 1e-15 {
        return 0.0;
    }

    let avg_angle = total_angle / count as f64;
    2.0 * PI / avg_angle
}

/// Compute the rise per base pair (nm) from backbone positions.
pub fn rise_per_bp(backbone: &[[f64; 3]]) -> f64 {
    if backbone.len() < 2 {
        return 0.0;
    }
    let mut total_rise = 0.0;
    for w in backbone.windows(2) {
        total_rise += (w[1][2] - w[0][2]).abs();
    }
    total_rise / (backbone.len() - 1) as f64
}

// ─────────────────────────────────────────────────────────────────────────────
// Poland-Scheraga denaturation model
// ─────────────────────────────────────────────────────────────────────────────

/// Poland-Scheraga model parameters for DNA denaturation.
#[derive(Debug, Clone)]
pub struct PolandScheragaParams {
    /// Loop entropy exponent c (typically ~1.75 for 3D).
    pub loop_exponent: f64,
    /// Cooperativity parameter sigma (typically ~1e-5 for DNA).
    pub cooperativity: f64,
    /// Base-pair free energy at reference temperature (kcal/mol).
    pub bp_free_energy: f64,
    /// Temperature (Kelvin).
    pub temperature: f64,
}

impl Default for PolandScheragaParams {
    fn default() -> Self {
        Self {
            loop_exponent: 1.75,
            cooperativity: 1e-5,
            bp_free_energy: -1.5,
            temperature: 300.0,
        }
    }
}

/// Partition function contribution from a loop of length `m` in the
/// Poland-Scheraga model.
///
/// Z_loop(m) = sigma * m^(-c)
pub fn loop_partition(m: usize, params: &PolandScheragaParams) -> f64 {
    if m == 0 {
        return 0.0;
    }
    params.cooperativity * (m as f64).powf(-params.loop_exponent)
}

/// Boltzmann weight for `k` consecutive base pairs at temperature `T`.
///
/// w = exp(-dG_bp / (k_B T))
pub fn bp_boltzmann_weight(k: usize, params: &PolandScheragaParams) -> f64 {
    let kb_kcal = 0.001987; // kcal/(mol K)
    let beta = 1.0 / (kb_kcal * params.temperature);
    (-(k as f64) * params.bp_free_energy * beta).exp()
}

/// Compute the fraction of closed (paired) base pairs in the Poland-Scheraga
/// model using a transfer-matrix approximation for a homopolymer of length `n`.
pub fn ps_fraction_closed(n: usize, params: &PolandScheragaParams) -> f64 {
    if n == 0 {
        return 0.0;
    }

    // Simple mean-field approximation:
    // theta = w / (w + sigma * sum_{m=1}^{n} m^{-c})
    let w = bp_boltzmann_weight(1, params);
    let loop_sum: f64 = (1..=n)
        .map(|m| (m as f64).powf(-params.loop_exponent))
        .sum();
    let denom = w + params.cooperativity * loop_sum;
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    (w / denom).clamp(0.0, 1.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Torsional / twist energy
// ─────────────────────────────────────────────────────────────────────────────

/// Harmonic twist energy per base-pair step.
///
/// E_twist = 0.5 * C / l * (omega - omega_0)^2
///
/// where C is the twist stiffness, l is the step height, omega is the twist
/// angle, and omega_0 is the natural twist angle.
pub fn twist_energy_per_step(
    omega: f64,
    omega_0: f64,
    twist_stiffness: f64,
    step_rise: f64,
) -> f64 {
    if step_rise <= 0.0 {
        return 0.0;
    }
    0.5 * twist_stiffness / step_rise * (omega - omega_0).powi(2)
}

/// Total twist energy for a stretch of `n` bp with uniform over/under-twist
/// `delta_omega` per step.
pub fn total_twist_energy(n: usize, delta_omega: f64, twist_stiffness: f64, rise: f64) -> f64 {
    n as f64 * twist_energy_per_step(delta_omega, 0.0, twist_stiffness, rise)
}

// ─────────────────────────────────────────────────────────────────────────────
// Private math helpers
// ─────────────────────────────────────────────────────────────────────────────

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = norm3(v);
    if len < 1e-15 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Project vector `v` out of direction `n` (assumed unit): v - (v.n)*n.
fn project_out(v: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    let dot = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
    [v[0] - dot * n[0], v[1] - dot * n[1], v[2] - dot * n[2]]
}

/// Signed solid angle subtended by two line segments (for writhe computation).
fn solid_angle_segment_pair(a1: [f64; 3], a2: [f64; 3], b1: [f64; 3], b2: [f64; 3]) -> f64 {
    let r12 = sub3(b1, a1);
    let r13 = sub3(b2, a1);
    let r14 = sub3(b2, a2);
    let r23 = sub3(b2, b1);
    let _r24 = sub3(a2, b1);

    let n1 = cross3(r12, r13);
    let n2 = cross3(r13, r14);

    let len1 = norm3(n1);
    let len2 = norm3(n2);
    if len1 < 1e-14 || len2 < 1e-14 {
        return 0.0;
    }

    let cos_omega = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (len1 * len2);
    let cos_omega = cos_omega.clamp(-1.0, 1.0);

    // Sign from the triple product
    let sign = n1[0] * r23[0] + n1[1] * r23[1] + n1[2] * r23[2];
    cos_omega.acos().copysign(sign)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nucleobase_complement() {
        assert_eq!(NucleoBase::Adenine.complement(), NucleoBase::Thymine);
        assert_eq!(NucleoBase::Thymine.complement(), NucleoBase::Adenine);
        assert_eq!(NucleoBase::Guanine.complement(), NucleoBase::Cytosine);
        assert_eq!(NucleoBase::Cytosine.complement(), NucleoBase::Guanine);
        assert_eq!(NucleoBase::Uracil.complement(), NucleoBase::Adenine);
    }

    #[test]
    fn test_rna_complement() {
        assert_eq!(NucleoBase::Adenine.rna_complement(), NucleoBase::Uracil);
        assert_eq!(NucleoBase::Uracil.rna_complement(), NucleoBase::Adenine);
    }

    #[test]
    fn test_nucleobase_code_roundtrip() {
        for &b in &[
            NucleoBase::Adenine,
            NucleoBase::Thymine,
            NucleoBase::Guanine,
            NucleoBase::Cytosine,
            NucleoBase::Uracil,
        ] {
            let c = b.code();
            assert_eq!(NucleoBase::from_code(c), Some(b));
        }
    }

    #[test]
    fn test_nucleobase_purine_pyrimidine() {
        assert!(NucleoBase::Adenine.is_purine());
        assert!(NucleoBase::Guanine.is_purine());
        assert!(NucleoBase::Cytosine.is_pyrimidine());
        assert!(NucleoBase::Thymine.is_pyrimidine());
        assert!(NucleoBase::Uracil.is_pyrimidine());
    }

    #[test]
    fn test_fene_potential_at_equilibrium() {
        let params = OxDnaParams::default();
        let e = fene_potential(params.fene_r0, &params);
        assert!(
            e.abs() < 1e-10,
            "FENE at equilibrium should be ~0, got {:.6}",
            e
        );
    }

    #[test]
    fn test_fene_potential_increases_with_stretch() {
        let params = OxDnaParams::default();
        let e_eq = fene_potential(params.fene_r0, &params);
        let e_stretch = fene_potential(params.fene_r0 + 0.1, &params);
        assert!(e_stretch > e_eq);
    }

    #[test]
    fn test_fene_diverges_at_max_extension() {
        let params = OxDnaParams::default();
        let e = fene_potential(params.fene_r0 + params.fene_delta + 0.01, &params);
        assert!(e > 1e5, "FENE should diverge beyond max extension");
    }

    #[test]
    fn test_hydrogen_bond_complementary_pair() {
        let params = OxDnaParams::default();
        let site_i = [0.0, 0.0, 0.0];
        let site_j = [0.2, 0.0, 0.0];
        let normal_i = [1.0, 0.0, 0.0];
        let normal_j = [-1.0, 0.0, 0.0]; // antiparallel
        let e = hydrogen_bond_energy(
            site_i,
            site_j,
            normal_i,
            normal_j,
            NucleoBase::Adenine,
            NucleoBase::Thymine,
            &params,
        );
        assert!(
            e < 0.0,
            "H-bond energy should be negative (attractive), got {:.6}",
            e
        );
    }

    #[test]
    fn test_hydrogen_bond_non_complementary() {
        let params = OxDnaParams::default();
        let e = hydrogen_bond_energy(
            [0.0; 3],
            [0.2, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            NucleoBase::Adenine,
            NucleoBase::Guanine,
            &params,
        );
        assert!(
            (e - 0.0).abs() < 1e-15,
            "Non-complementary should give 0 energy"
        );
    }

    #[test]
    fn test_gc_stronger_than_at() {
        let params = OxDnaParams::default();
        let site_i = [0.0, 0.0, 0.0];
        let site_j = [0.2, 0.0, 0.0];
        let n_i = [1.0, 0.0, 0.0];
        let n_j = [-1.0, 0.0, 0.0];
        let e_at = hydrogen_bond_energy(
            site_i,
            site_j,
            n_i,
            n_j,
            NucleoBase::Adenine,
            NucleoBase::Thymine,
            &params,
        );
        let e_gc = hydrogen_bond_energy(
            site_i,
            site_j,
            n_i,
            n_j,
            NucleoBase::Guanine,
            NucleoBase::Cytosine,
            &params,
        );
        assert!(e_gc < e_at, "GC should be stronger (more negative) than AT");
    }

    #[test]
    fn test_stacking_energy_aligned() {
        let params = OxDnaParams::default();
        let s1 = [0.0, 0.0, 0.0];
        let s2 = [0.0, 0.0, 0.34];
        let t = [0.0, 0.0, 1.0];
        let e = stacking_energy(s1, s2, t, t, &params);
        assert!(
            e < 0.0,
            "Aligned stacking should be attractive, got {:.6}",
            e
        );
    }

    #[test]
    fn test_excluded_volume_repulsive() {
        let params = OxDnaParams::default();
        let e = excluded_volume_energy(0.1, &params);
        assert!(
            e > 0.0,
            "Excluded volume should be repulsive at short range"
        );
    }

    #[test]
    fn test_excluded_volume_zero_beyond_sigma() {
        let params = OxDnaParams::default();
        let e = excluded_volume_energy(params.exc_vol_sigma + 0.1, &params);
        assert!((e - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_strand_from_sequence() {
        let strand = OxDnaStrand::from_sequence("ATCG", 0);
        assert_eq!(strand.len(), 4);
        assert_eq!(strand.sequence_string(), "ATCG");
        assert!(!strand.is_empty());
    }

    #[test]
    fn test_strand_gc_content() {
        let strand = OxDnaStrand::from_sequence("GGCC", 0);
        assert!((strand.gc_content() - 1.0).abs() < 1e-10);
        let strand2 = OxDnaStrand::from_sequence("AATT", 0);
        assert!((strand2.gc_content() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_strand_contour_length_positive() {
        let strand = OxDnaStrand::from_sequence("ATCGATCG", 0);
        assert!(strand.contour_length() > 0.0);
    }

    #[test]
    fn test_duplex_construction() {
        let duplex = OxDnaDuplex::from_sense_sequence("ATCGATCG");
        assert_eq!(duplex.num_base_pairs(), 8);
        assert_eq!(duplex.sense.sequence_string(), "ATCGATCG");
    }

    #[test]
    fn test_duplex_total_energy_finite() {
        let duplex = OxDnaDuplex::from_sense_sequence("ATCGATCG");
        let e = duplex.total_potential_energy();
        assert!(e.is_finite(), "Total energy should be finite, got {:.6}", e);
    }

    #[test]
    fn test_melting_temperature_positive() {
        let tm = melting_temperature("ATCGATCGATCG", 1e-4, 1.0);
        assert!(tm > 200.0, "Tm should be above 200 K, got {:.6}", tm);
        assert!(tm < 400.0, "Tm should be below 400 K, got {:.6}", tm);
    }

    #[test]
    fn test_melting_temp_gc_higher_than_at() {
        let tm_gc = melting_temperature("GGGGCCCC", 1e-4, 1.0);
        let tm_at = melting_temperature("AAAATTTT", 1e-4, 1.0);
        assert!(tm_gc > tm_at, "GC-rich Tm should be higher than AT-rich");
    }

    #[test]
    fn test_melting_temp_short_seq() {
        let tm = melting_temperature("A", 1e-4, 1.0);
        assert!((tm - 0.0).abs() < 1e-10, "Single base should give Tm=0");
    }

    #[test]
    fn test_persistence_length_straight_chain() {
        // Perfectly aligned tangents => infinite persistence length
        let tangents: Vec<[f64; 3]> = (0..50).map(|_| [0.0, 0.0, 1.0]).collect();
        let lp = persistence_length(&tangents, 0.34);
        assert!(
            lp > 10.0,
            "Straight chain should have large lp, got {:.6}",
            lp
        );
    }

    #[test]
    fn test_wlc_end_to_end() {
        let r = wlc_end_to_end_rms(100.0, 50.0);
        assert!(r > 0.0 && r < 100.0, "WLC R should be between 0 and L");
    }

    #[test]
    fn test_wlc_zero_length() {
        assert!((wlc_end_to_end_rms(0.0, 50.0) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_linking_decomposition_circle() {
        // Flat circle: twist=0, writhe=0
        let n = 20;
        let backbone: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                [theta.cos(), theta.sin(), 0.0]
            })
            .collect();
        let normals: Vec<[f64; 3]> = vec![[0.0, 0.0, 1.0]; n];
        let ld = linking_decomposition(&backbone, &normals);
        assert!(
            ld.writhe.abs() < 0.5,
            "Flat circle writhe should be ~0, got {:.6}",
            ld.writhe
        );
    }

    #[test]
    fn test_superhelical_density() {
        let sigma = superhelical_density(100.0, 100.0);
        assert!((sigma - 0.0).abs() < 1e-15);
        let sigma2 = superhelical_density(95.0, 100.0);
        assert!((sigma2 - (-0.05)).abs() < 1e-10);
    }

    #[test]
    fn test_denaturation_bubbles_stable_duplex() {
        let duplex = OxDnaDuplex::from_sense_sequence("GCGCGCGCGC");
        let _bubbles = detect_denaturation_bubbles(&duplex, 0.001);
        // With this threshold many pairs may be open; just check no panic
        let _ = denaturation_fraction(&duplex, 0.5);
    }

    #[test]
    fn test_nucleosome_wrapped_bp() {
        let params = NucleosomeParams::default();
        let bp = nucleosome_wrapped_bp(&params);
        // ~147 bp for 1.65 turns around 4.18 nm radius
        assert!(bp > 100.0, "Expected >100 wrapped bp, got {:.6}", bp);
    }

    #[test]
    fn test_nucleosome_energy_balance() {
        let params = NucleosomeParams::default();
        let e_bend = nucleosome_bending_energy(&params);
        let e_ads = nucleosome_adsorption_energy(&params);
        assert!(e_bend > 0.0, "Bending energy should be positive");
        assert!(e_ads < 0.0, "Adsorption energy should be negative");
    }

    #[test]
    fn test_equilibrium_unwrap_angle() {
        let params = NucleosomeParams::default();
        let a0 = equilibrium_unwrap_angle(0.0, &params);
        assert!((a0 - 0.0).abs() < 1e-15);
        let a1 = equilibrium_unwrap_angle(10.0, &params);
        assert!(a1 > 0.0);
    }

    #[test]
    fn test_ideal_bdna_backbone() {
        let (sense, anti) = ideal_bdna_backbone(20);
        assert_eq!(sense.len(), 20);
        assert_eq!(anti.len(), 20);
    }

    #[test]
    fn test_helical_repeat_ideal() {
        let (sense, _) = ideal_bdna_backbone(100);
        let hr = helical_repeat(&sense);
        // Should be close to 10.5
        assert!(
            (hr - 10.5).abs() < 1.0,
            "Helical repeat should be ~10.5, got {:.6}",
            hr
        );
    }

    #[test]
    fn test_rise_per_bp_ideal() {
        let (sense, _) = ideal_bdna_backbone(50);
        let r = rise_per_bp(&sense);
        assert!(
            (r - 0.34).abs() < 0.01,
            "Rise per bp should be ~0.34, got {:.6}",
            r
        );
    }

    #[test]
    fn test_dinucleotide_counts() {
        let counts = dinucleotide_counts("AATCG");
        // AA=1, AT=1, TC=1, CG=1
        assert_eq!(counts[0], 1); // AA
        assert_eq!(counts[1], 1); // AT
        assert_eq!(counts[4 + 3], 1); // TC
        assert_eq!(counts[3 * 4 + 2], 1); // CG
    }

    #[test]
    fn test_nn_free_energy_negative() {
        let dg = nn_free_energy("ATCGATCG", 300.0);
        // At room temperature, duplex formation should be favourable
        assert!(dg < 0.0, "Free energy should be negative, got {:.6}", dg);
    }

    #[test]
    fn test_is_self_complementary() {
        assert!(is_self_complementary("AATT"));
        assert!(!is_self_complementary("AATC"));
    }

    #[test]
    fn test_cpg_count() {
        assert_eq!(cpg_count("ACGCGA"), 2);
        assert_eq!(cpg_count("AAAA"), 0);
    }

    #[test]
    fn test_ps_fraction_closed() {
        let params = PolandScheragaParams::default();
        let f = ps_fraction_closed(100, &params);
        assert!(
            (0.0..=1.0).contains(&f),
            "Fraction should be in [0,1], got {:.6}",
            f
        );
    }

    #[test]
    fn test_twist_energy() {
        let e = twist_energy_per_step(0.6, 0.6, 100.0, 0.34);
        assert!((e - 0.0).abs() < 1e-10);
        let e2 = twist_energy_per_step(0.7, 0.6, 100.0, 0.34);
        assert!(e2 > 0.0);
    }

    #[test]
    fn test_apply_backbone_forces_no_panic() {
        let mut strand = OxDnaStrand::from_sequence("ATCG", 0);
        let params = OxDnaParams::default();
        apply_backbone_forces(&mut strand, &params, 0.001, 1.0);
        advance_positions(&mut strand, 0.001);
        let ke = kinetic_energy(&strand, 1.0);
        assert!(ke >= 0.0);
    }

    #[test]
    fn test_instantaneous_temperature() {
        let strand = OxDnaStrand::from_sequence("ATCG", 0);
        let t = instantaneous_temperature(&strand, 1.0);
        assert!((t - 0.0).abs() < 1e-10, "Zero velocity => T=0");
    }
}
