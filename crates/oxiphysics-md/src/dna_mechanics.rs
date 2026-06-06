// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! DNA and nucleic acid mechanics for molecular dynamics.
//!
//! Provides models for B-form double-stranded DNA, the worm-like chain (WLC)
//! polymer model, DNA origami scaffolding, and Coxeter helix geometry.
//! Includes nearest-neighbour stacking energies and hydrogen-bond potentials.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Base kind
// ─────────────────────────────────────────────────────────────────────────────

/// The four canonical DNA nucleobases.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaseKind {
    /// Adenine (purine, pairs with Thymine).
    Adenine,
    /// Thymine (pyrimidine, pairs with Adenine).
    Thymine,
    /// Guanine (purine, pairs with Cytosine).
    Guanine,
    /// Cytosine (pyrimidine, pairs with Guanine).
    Cytosine,
}

impl BaseKind {
    /// Return the Watson-Crick complement of this base.
    pub fn complement(self) -> Self {
        match self {
            BaseKind::Adenine => BaseKind::Thymine,
            BaseKind::Thymine => BaseKind::Adenine,
            BaseKind::Guanine => BaseKind::Cytosine,
            BaseKind::Cytosine => BaseKind::Guanine,
        }
    }

    /// One-letter code for this base.
    pub fn code(self) -> char {
        match self {
            BaseKind::Adenine => 'A',
            BaseKind::Thymine => 'T',
            BaseKind::Guanine => 'G',
            BaseKind::Cytosine => 'C',
        }
    }

    /// Parse a one-letter code (case-insensitive) into a `BaseKind`.
    ///
    /// Returns `None` for unrecognised characters.
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(BaseKind::Adenine),
            'T' => Some(BaseKind::Thymine),
            'G' => Some(BaseKind::Guanine),
            'C' => Some(BaseKind::Cytosine),
            _ => None,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DNABase
// ─────────────────────────────────────────────────────────────────────────────

/// A single DNA nucleobase with position and orientation.
///
/// `orientation` is stored as a unit quaternion `[w, x, y, z]`.
#[derive(Debug, Clone)]
pub struct DNABase {
    /// Identity of this nucleobase.
    pub kind: BaseKind,
    /// Centre-of-mass position in Å.
    pub position: [f64; 3],
    /// Orientation as a unit quaternion \[w, x, y, z\].
    pub orientation: [f64; 4],
}

impl DNABase {
    /// Create a new base at the given position with identity orientation.
    pub fn new(kind: BaseKind, position: [f64; 3]) -> Self {
        Self {
            kind,
            position,
            orientation: [1.0, 0.0, 0.0, 0.0],
        }
    }

    /// Create a base with an explicit quaternion orientation.
    pub fn with_orientation(kind: BaseKind, position: [f64; 3], orientation: [f64; 4]) -> Self {
        Self {
            kind,
            position,
            orientation,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DNADuplex
// ─────────────────────────────────────────────────────────────────────────────

/// A double-stranded DNA duplex modelled in B-form geometry.
///
/// The duplex is characterised by its base-pair sequence, helical rise per
/// base pair, and helical twist per base pair.
#[derive(Debug, Clone)]
pub struct DNADuplex {
    /// Bases of the sense strand in 5'→3' order.
    pub bases: Vec<DNABase>,
    /// Helical rise per base pair in Å (B-DNA default ≈ 3.4 Å).
    pub rise: f64,
    /// Helical twist per base pair in radians (B-DNA default ≈ 36°).
    pub twist: f64,
    /// Number of base pairs.
    pub n_bp: usize,
}

/// B-DNA rise per base pair in Å.
const B_DNA_RISE: f64 = 3.4;
/// B-DNA twist per base pair in degrees.
const B_DNA_TWIST_DEG: f64 = 36.0;
/// B-DNA helix radius in Å.
const B_DNA_RADIUS: f64 = 10.0;
/// Persistence length of dsDNA in nm.
const DNA_LP_NM: f64 = 50.0;

impl DNADuplex {
    /// Build a B-form duplex from a 5'→3' sequence string.
    ///
    /// Unrecognised characters are silently skipped.  The duplex uses
    /// standard B-DNA parameters: rise = 3.4 Å, twist = 36°/bp.
    pub fn new(sequence: &str) -> Self {
        let rise = B_DNA_RISE;
        let twist = B_DNA_TWIST_DEG.to_radians();
        let helix = CoxeterHelixParameters {
            rise_per_bp: rise,
            twist_per_bp: twist,
            radius: B_DNA_RADIUS,
        };
        let mut bases = Vec::new();
        for (i, c) in sequence.chars().enumerate() {
            if let Some(kind) = BaseKind::from_char(c) {
                let pos = helix.helix_position(i, 0);
                bases.push(DNABase::new(kind, pos));
            }
        }
        let n_bp = bases.len();
        Self {
            bases,
            rise,
            twist,
            n_bp,
        }
    }

    /// Persistence length of dsDNA in Å (50 nm = 500 Å).
    pub fn persistence_length(&self) -> f64 {
        DNA_LP_NM * 10.0 // convert nm → Å
    }

    /// Contour length of the duplex in Å.
    ///
    /// `L₀ = n_bp × rise`
    pub fn contour_length(&self) -> f64 {
        self.n_bp as f64 * self.rise
    }

    /// Fractional extension under a stretching force using the Marko-Siggia WLC.
    ///
    /// Returns the end-to-end extension `z` in Å for an applied force `force`
    /// in pN.  Uses the worm-like chain via [`WormLikeChain::wlc_extension`].
    pub fn wlc_extension(&self, force: f64) -> f64 {
        let wlc = WormLikeChain {
            lp: self.persistence_length(),
            l0: self.contour_length(),
        };
        wlc.wlc_extension(force)
    }

    /// Mean bending angle (in radians) over the full contour.
    ///
    /// Derived from the worm-like chain thermal-fluctuation formula:
    /// `⟨θ²⟩ = L / Lp` for small deflections.
    pub fn bend_angle(&self) -> f64 {
        let l = self.contour_length();
        let lp = self.persistence_length();
        (l / lp).sqrt()
    }

    /// Total twist angle over the full duplex in radians.
    pub fn twist_angle(&self) -> f64 {
        self.n_bp as f64 * self.twist
    }

    /// Width of the major groove in Å (B-DNA ≈ 11.7 Å).
    pub fn groove_width_major(&self) -> f64 {
        11.7
    }

    /// Width of the minor groove in Å (B-DNA ≈ 5.7 Å).
    pub fn groove_width_minor(&self) -> f64 {
        5.7
    }

    /// GC-content of the duplex sequence.
    pub fn gc_content(&self) -> f64 {
        if self.bases.is_empty() {
            return 0.0;
        }
        let gc = self
            .bases
            .iter()
            .filter(|b| b.kind == BaseKind::Guanine || b.kind == BaseKind::Cytosine)
            .count();
        gc as f64 / self.bases.len() as f64
    }

    /// Melting temperature estimate (°C) using the Wallace rule for short duplexes.
    ///
    /// `Tm = 2(A+T) + 4(G+C)` \[°C\].
    pub fn melting_temperature(&self) -> f64 {
        let at = self
            .bases
            .iter()
            .filter(|b| b.kind == BaseKind::Adenine || b.kind == BaseKind::Thymine)
            .count();
        let gc = self
            .bases
            .iter()
            .filter(|b| b.kind == BaseKind::Guanine || b.kind == BaseKind::Cytosine)
            .count();
        2.0 * at as f64 + 4.0 * gc as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// WormLikeChain
// ─────────────────────────────────────────────────────────────────────────────

/// Worm-like chain (WLC) polymer model for semiflexible polymers such as DNA.
///
/// The persistence length `lp` and contour length `l0` are both in Å.
/// Forces are in pN.  Uses kT = 4.11 pN·nm = 41.1 pN·Å at 300 K.
#[derive(Debug, Clone)]
pub struct WormLikeChain {
    /// Persistence length in Å.
    pub lp: f64,
    /// Contour length in Å.
    pub l0: f64,
}

/// Thermal energy kT at 300 K in pN·Å.
const KT_300K_PN_ANG: f64 = 41.1;

impl WormLikeChain {
    /// Create a WLC with given persistence and contour lengths (Å).
    pub fn new(lp: f64, l0: f64) -> Self {
        Self { lp, l0 }
    }

    /// Force (pN) required to extend the chain to fractional extension `z / l0`.
    ///
    /// Uses the exact Marko-Siggia interpolation formula.
    /// `z` is the end-to-end distance in Å; clipped to (0, l0).
    pub fn force_extension(&self, z: f64) -> f64 {
        let x = (z / self.l0).clamp(0.001, 0.999);
        let kt = KT_300K_PN_ANG;
        let lp = self.lp;
        (kt / lp) * (0.25 / (1.0 - x).powi(2) - 0.25 + x)
    }

    /// Marko-Siggia approximate inverse: extension `z` in Å for force `f` (pN).
    ///
    /// Iterates Newton's method to invert `force_extension`.
    pub fn marko_siggia_approx(&self, f: f64) -> f64 {
        if f <= 0.0 {
            return 0.0;
        }
        // Initial guess from the high-force asymptotic: x ≈ 1 - sqrt(kt/(4 lp f))
        let kt = KT_300K_PN_ANG;
        let x0 = 1.0 - (kt / (4.0 * self.lp * f)).sqrt().clamp(0.001, 0.998);
        let mut x = x0.clamp(0.001, 0.998);
        // Newton iterations
        for _ in 0..50 {
            let fx = self.force_extension(x * self.l0) - f;
            // df/dx = (kt/lp)*(0.5/(1-x)^3 + 1)
            let dfx = (kt / self.lp) * (0.5 / (1.0 - x).powi(3) + 1.0);
            let dx = fx / dfx;
            x -= dx;
            x = x.clamp(0.001, 0.998);
            if dx.abs() < 1e-9 {
                break;
            }
        }
        x * self.l0
    }

    /// Convenience alias: extension (Å) for a force `f` (pN).
    pub fn wlc_extension(&self, f: f64) -> f64 {
        self.marko_siggia_approx(f)
    }

    /// Extensible WLC model incorporating backbone stretching.
    ///
    /// `z` is the end-to-end distance (Å), `k_stretch` is the stretch modulus
    /// in pN (typically ~1000 pN for dsDNA).
    ///
    /// The WLC force is augmented by an entropic correction for over-extension.
    pub fn extensible_wlc(&self, z: f64, k_stretch: f64) -> f64 {
        // Base WLC force
        let f_wlc = self.force_extension(z);
        // Backbone stretching: if z > l0, add harmonic correction
        let excess = (z - self.l0).max(0.0);
        f_wlc + k_stretch * excess / self.l0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DNA Origami
// ─────────────────────────────────────────────────────────────────────────────

/// A simplified model of a DNA origami structure.
///
/// Strands are represented as ordered lists of base-pair indices into a
/// scaffold of `n_bp` base pairs.
#[derive(Debug, Clone)]
pub struct DNAOrigami {
    /// Scaffold length in base pairs.
    pub n_bp: usize,
    /// Staple-strand definitions; each inner `Vec`usize` lists the scaffold
    /// positions spanned by one staple.
    pub strands: Vec<Vec<usize>>,
}

impl DNAOrigami {
    /// Create an origami model with a given scaffold length.
    pub fn new(n_bp: usize) -> Self {
        Self {
            n_bp,
            strands: Vec::new(),
        }
    }

    /// Add a staple strand spanning `positions` on the scaffold.
    pub fn add_staple(&mut self, positions: Vec<usize>) {
        self.strands.push(positions);
    }

    /// Generate a simple linear scaffold routing (sequential base pairs).
    ///
    /// Returns the scaffold as a `Vec`usize` of base-pair indices 0..n_bp.
    pub fn scaffold_routing(&self) -> Vec<usize> {
        (0..self.n_bp).collect()
    }

    /// Total number of staple strands.
    pub fn staple_count(&self) -> usize {
        self.strands.len()
    }

    /// Estimate the number of staple strands needed to tile the scaffold.
    ///
    /// Uses the standard rule-of-thumb: one staple every ~32 bp covering
    /// two scaffold helices (i.e., scaffold length / 32).
    pub fn estimate_staple_count(&self) -> usize {
        self.n_bp.div_ceil(32)
    }

    /// Check whether all scaffold positions are covered by at least one staple.
    pub fn is_fully_covered(&self) -> bool {
        let mut covered = vec![false; self.n_bp];
        for strand in &self.strands {
            for &idx in strand {
                if idx < self.n_bp {
                    covered[idx] = true;
                }
            }
        }
        covered.iter().all(|&c| c)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Coxeter Helix Parameters
// ─────────────────────────────────────────────────────────────────────────────

/// Geometric parameters for an ideal nucleic acid helix (Coxeter helical model).
///
/// Positions along each strand are computed from the rise, twist, and radius.
#[derive(Debug, Clone)]
pub struct CoxeterHelixParameters {
    /// Rise per base pair in Å.
    pub rise_per_bp: f64,
    /// Twist per base pair in radians.
    pub twist_per_bp: f64,
    /// Helix radius in Å.
    pub radius: f64,
}

impl CoxeterHelixParameters {
    /// Create standard B-DNA helix parameters.
    pub fn b_dna() -> Self {
        Self {
            rise_per_bp: B_DNA_RISE,
            twist_per_bp: B_DNA_TWIST_DEG.to_radians(),
            radius: B_DNA_RADIUS,
        }
    }

    /// Position of base pair `bp` on strand `strand` (0 = sense, 1 = antisense).
    ///
    /// Returns `[x, y, z]` in Å.  The antisense strand is offset by π radians.
    pub fn helix_position(&self, bp: usize, strand: i32) -> [f64; 3] {
        let phi = bp as f64 * self.twist_per_bp + strand as f64 * PI;
        let z = bp as f64 * self.rise_per_bp;
        [self.radius * phi.cos(), self.radius * phi.sin(), z]
    }

    /// Axial length of a helix containing `n_bp` base pairs in Å.
    pub fn axial_length(&self, n_bp: usize) -> f64 {
        n_bp as f64 * self.rise_per_bp
    }

    /// Total twist angle of a helix containing `n_bp` base pairs in radians.
    pub fn total_twist(&self, n_bp: usize) -> f64 {
        n_bp as f64 * self.twist_per_bp
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Nearest-neighbour stacking energies
// ─────────────────────────────────────────────────────────────────────────────

/// Base-stacking energy (kcal/mol) for two consecutive bases on the same strand.
///
/// Uses the SantaLucia 1998 unified nearest-neighbour thermodynamic parameters.
/// The dinucleotide is written 5′→3′ (base1 then base2).
///
/// # Examples
/// ```no_run
/// use oxiphysics_md::dna_mechanics::{stacking_energy, BaseKind};
/// let e = stacking_energy(BaseKind::Adenine, BaseKind::Thymine);
/// assert!(e < 0.0, "stacking should be stabilising");
/// ```
pub fn stacking_energy(base1: BaseKind, base2: BaseKind) -> f64 {
    // SantaLucia 1998 ΔH° values (kcal/mol) for all 16 dinucleotides.
    // Key: (5' base, 3' base) → ΔH°.
    use BaseKind::*;
    match (base1, base2) {
        (Adenine, Adenine) => -7.9,
        (Adenine, Thymine) => -7.2,
        (Adenine, Guanine) => -7.8,
        (Adenine, Cytosine) => -7.8,
        (Thymine, Adenine) => -7.2,
        (Thymine, Thymine) => -7.9,
        (Thymine, Guanine) => -8.5,
        (Thymine, Cytosine) => -8.2,
        (Guanine, Adenine) => -8.2,
        (Guanine, Thymine) => -8.4,
        (Guanine, Guanine) => -8.0,
        (Guanine, Cytosine) => -10.6,
        (Cytosine, Adenine) => -7.8,
        (Cytosine, Thymine) => -7.8,
        (Cytosine, Guanine) => -9.8,
        (Cytosine, Cytosine) => -8.0,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Hydrogen bond energy
// ─────────────────────────────────────────────────────────────────────────────

/// Hydrogen-bond energy (kcal/mol) as a function of inter-base-pair distance.
///
/// Uses a Morse potential centred at 3.0 Å equilibrium separation with well
/// depth 3.0 kcal/mol (typical Watson-Crick H-bond) and decay constant
/// α = 2.0 Å⁻¹.
///
/// `E(r) = D · (1 − exp(−α · (r − r₀)))² − D`
///
/// * At r = r₀: E = −D (energy minimum, attractive).
/// * At r → ∞: E = 0.
/// * At r ≪ r₀: E > 0 (repulsive wall).
///
/// # Parameters
/// * `bp_distance` – distance between the two hydrogen-bond donor/acceptor
///   centres in Å.
pub fn hydrogen_bond_energy(bp_distance: f64) -> f64 {
    let r0 = 3.0; // equilibrium separation, Å
    let d = 3.0; // well depth (positive), kcal/mol
    let a = 2.0; // Morse α, Å⁻¹
    let x = bp_distance - r0;
    d * (1.0 - (-a * x).exp()).powi(2) - d
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper geometry
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the distance between two positions in Å.
#[cfg(test)]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Dot product of two quaternions \[w, x, y, z\].
fn quat_dot(q: [f64; 4], r: [f64; 4]) -> f64 {
    q[0] * r[0] + q[1] * r[1] + q[2] * r[2] + q[3] * r[3]
}

/// Angular distance between two orientations represented as unit quaternions.
///
/// Returns the angle in radians ∈ \[0, π\].
pub fn quaternion_angle(q: [f64; 4], r: [f64; 4]) -> f64 {
    let dot = quat_dot(q, r).abs().clamp(-1.0, 1.0);
    2.0 * dot.acos()
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::BaseKind::*;
    use super::*;

    // ── BaseKind ─────────────────────────────────────────────────────────────

    #[test]
    fn test_complement_at() {
        assert_eq!(BaseKind::Adenine.complement(), BaseKind::Thymine);
        assert_eq!(BaseKind::Thymine.complement(), BaseKind::Adenine);
    }

    #[test]
    fn test_complement_gc() {
        assert_eq!(BaseKind::Guanine.complement(), BaseKind::Cytosine);
        assert_eq!(BaseKind::Cytosine.complement(), BaseKind::Guanine);
    }

    #[test]
    fn test_base_code() {
        assert_eq!(BaseKind::Adenine.code(), 'A');
        assert_eq!(BaseKind::Thymine.code(), 'T');
        assert_eq!(BaseKind::Guanine.code(), 'G');
        assert_eq!(BaseKind::Cytosine.code(), 'C');
    }

    #[test]
    fn test_base_from_char() {
        assert_eq!(BaseKind::from_char('a'), Some(BaseKind::Adenine));
        assert_eq!(BaseKind::from_char('X'), None);
    }

    // ── DNABase ──────────────────────────────────────────────────────────────

    #[test]
    fn test_dnabase_new() {
        let b = DNABase::new(BaseKind::Guanine, [1.0, 2.0, 3.0]);
        assert_eq!(b.kind, BaseKind::Guanine);
        assert_eq!(b.orientation, [1.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_dnabase_with_orientation() {
        let q = [0.0_f64, 1.0, 0.0, 0.0];
        let b = DNABase::with_orientation(BaseKind::Adenine, [0.0; 3], q);
        assert_eq!(b.orientation, q);
    }

    // ── DNADuplex ────────────────────────────────────────────────────────────

    #[test]
    fn test_duplex_new_sequence() {
        let d = DNADuplex::new("ATCG");
        assert_eq!(d.n_bp, 4);
        assert_eq!(d.bases[0].kind, BaseKind::Adenine);
        assert_eq!(d.bases[3].kind, BaseKind::Guanine);
    }

    #[test]
    fn test_duplex_persistence_length() {
        let d = DNADuplex::new("ATCGATCG");
        assert!((d.persistence_length() - 500.0).abs() < 1.0);
    }

    #[test]
    fn test_duplex_contour_length() {
        let d = DNADuplex::new("AAAAAA"); // 6 bp
        let expected = 6.0 * 3.4;
        assert!((d.contour_length() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_duplex_twist_angle() {
        let d = DNADuplex::new("AAAA"); // 4 bp
        let expected = 4.0 * 36_f64.to_radians();
        assert!((d.twist_angle() - expected).abs() < 1e-10);
    }

    #[test]
    fn test_groove_widths() {
        let d = DNADuplex::new("ATCG");
        assert!((d.groove_width_major() - 11.7).abs() < 0.1);
        assert!((d.groove_width_minor() - 5.7).abs() < 0.1);
    }

    #[test]
    fn test_gc_content() {
        let d = DNADuplex::new("GGCC"); // 100% GC
        assert!((d.gc_content() - 1.0).abs() < 1e-10);
        let d2 = DNADuplex::new("AATT"); // 0% GC
        assert!(d2.gc_content() < 1e-10);
    }

    #[test]
    fn test_melting_temperature_wallace() {
        // AATT: 4 AT pairs → 2*4 = 8°C
        let d = DNADuplex::new("AATT");
        assert!((d.melting_temperature() - 8.0).abs() < 1e-10);
        // GGCC: 4 GC pairs → 4*4 = 16°C
        let d2 = DNADuplex::new("GGCC");
        assert!((d2.melting_temperature() - 16.0).abs() < 1e-10);
    }

    #[test]
    fn test_wlc_extension_positive() {
        let d = DNADuplex::new(&"A".repeat(100));
        let z = d.wlc_extension(0.5); // 0.5 pN
        assert!(z > 0.0);
        assert!(z < d.contour_length());
    }

    #[test]
    fn test_bend_angle_positive() {
        let d = DNADuplex::new(&"A".repeat(50));
        assert!(d.bend_angle() > 0.0);
    }

    // ── WormLikeChain ────────────────────────────────────────────────────────

    #[test]
    fn test_wlc_force_zero_extension() {
        let wlc = WormLikeChain::new(500.0, 3400.0);
        let f = wlc.force_extension(1.0); // very small x
        assert!(f >= 0.0);
    }

    #[test]
    fn test_wlc_force_increases_with_extension() {
        let wlc = WormLikeChain::new(500.0, 3400.0);
        let f1 = wlc.force_extension(0.3 * 3400.0);
        let f2 = wlc.force_extension(0.7 * 3400.0);
        assert!(f2 > f1);
    }

    #[test]
    fn test_marko_siggia_roundtrip() {
        let wlc = WormLikeChain::new(500.0, 3400.0);
        let f = 1.0; // pN
        let z = wlc.marko_siggia_approx(f);
        let f_back = wlc.force_extension(z);
        assert!(
            (f_back - f).abs() / f < 0.01,
            "roundtrip error: got {f_back} pN"
        );
    }

    #[test]
    fn test_wlc_extension_zero_force() {
        let wlc = WormLikeChain::new(500.0, 3400.0);
        let z = wlc.wlc_extension(0.0);
        assert_eq!(z, 0.0);
    }

    #[test]
    fn test_extensible_wlc_stretch() {
        let wlc = WormLikeChain::new(500.0, 3400.0);
        // At exactly l0, extensible WLC should equal base WLC
        let z = 0.9 * wlc.l0;
        let f_base = wlc.force_extension(z);
        let f_ext = wlc.extensible_wlc(z, 1000.0);
        assert!((f_ext - f_base).abs() < 1e-6);
        // Beyond l0, extensible should be larger
        let f_over = wlc.extensible_wlc(wlc.l0 * 1.01, 1000.0);
        assert!(f_over > f_base);
    }

    // ── DNAOrigami ───────────────────────────────────────────────────────────

    #[test]
    fn test_origami_staple_count() {
        let mut o = DNAOrigami::new(100);
        o.add_staple(vec![0, 1, 2]);
        o.add_staple(vec![3, 4, 5]);
        assert_eq!(o.staple_count(), 2);
    }

    #[test]
    fn test_origami_scaffold_routing() {
        let o = DNAOrigami::new(10);
        let r = o.scaffold_routing();
        assert_eq!(r, (0..10).collect::<Vec<_>>());
    }

    #[test]
    fn test_origami_estimate_staples() {
        let o = DNAOrigami::new(64);
        assert_eq!(o.estimate_staple_count(), 2);
    }

    #[test]
    fn test_origami_coverage() {
        let mut o = DNAOrigami::new(3);
        assert!(!o.is_fully_covered());
        o.add_staple(vec![0, 1, 2]);
        assert!(o.is_fully_covered());
    }

    // ── CoxeterHelixParameters ───────────────────────────────────────────────

    #[test]
    fn test_helix_position_strand0() {
        let h = CoxeterHelixParameters::b_dna();
        let p = h.helix_position(0, 0);
        // bp=0, phi=0: x=radius, y=0, z=0
        assert!((p[0] - B_DNA_RADIUS).abs() < 1e-10);
        assert!(p[1].abs() < 1e-10);
        assert!(p[2].abs() < 1e-10);
    }

    #[test]
    fn test_helix_position_z_increases() {
        let h = CoxeterHelixParameters::b_dna();
        let p0 = h.helix_position(0, 0);
        let p1 = h.helix_position(1, 0);
        assert!(p1[2] > p0[2]);
        assert!((p1[2] - p0[2] - B_DNA_RISE).abs() < 1e-10);
    }

    #[test]
    fn test_helix_axial_length() {
        let h = CoxeterHelixParameters::b_dna();
        assert!((h.axial_length(10) - 34.0).abs() < 1e-10);
    }

    #[test]
    fn test_helix_total_twist() {
        let h = CoxeterHelixParameters::b_dna();
        let t = h.total_twist(10);
        let expected = 10.0 * 36_f64.to_radians();
        assert!((t - expected).abs() < 1e-10);
    }

    // ── Stacking energy ──────────────────────────────────────────────────────

    #[test]
    fn test_stacking_energy_all_negative() {
        let bases = [Adenine, Thymine, Guanine, Cytosine];
        for &b1 in &bases {
            for &b2 in &bases {
                assert!(
                    stacking_energy(b1, b2) < 0.0,
                    "stacking energy should be negative for {:?}-{:?}",
                    b1,
                    b2
                );
            }
        }
    }

    #[test]
    fn test_stacking_gc_strongest() {
        // GC/GC should have the most negative stacking energy
        let e_gc = stacking_energy(BaseKind::Guanine, BaseKind::Cytosine);
        let e_at = stacking_energy(BaseKind::Adenine, BaseKind::Thymine);
        assert!(e_gc < e_at);
    }

    // ── Hydrogen bond energy ─────────────────────────────────────────────────

    #[test]
    fn test_hbond_at_equilibrium() {
        // At r = r0 = 3.0 Å, energy should be near zero (Morse minimum at 0)
        // Actually Morse: D*(1 - exp(0))^2 - D = 0 - D → correction: at r=r0 → exp(0)=1 → (1-1)^2=0 → -D
        // E(r0) = D*(0) - D = -D (the well bottom)
        let e = hydrogen_bond_energy(3.0);
        assert!(e < 0.0, "should be attractive at equilibrium distance");
    }

    #[test]
    fn test_hbond_large_distance_zero() {
        // At very large distances the energy should approach 0
        let e = hydrogen_bond_energy(100.0);
        assert!(
            e.abs() < 0.01,
            "energy should vanish at large distance, got {e}"
        );
    }

    #[test]
    fn test_hbond_repulsive_at_short() {
        // Very short distance: energy should be positive (repulsive)
        let e = hydrogen_bond_energy(0.5);
        assert!(e > 0.0, "energy should be repulsive at short distances");
    }

    // ── Quaternion angle ─────────────────────────────────────────────────────

    #[test]
    fn test_quaternion_angle_same() {
        let q = [1.0_f64, 0.0, 0.0, 0.0];
        let angle = quaternion_angle(q, q);
        assert!(
            angle.abs() < 1e-10,
            "angle between identical quaternions should be 0"
        );
    }

    #[test]
    fn test_quaternion_angle_opposite() {
        let q = [1.0_f64, 0.0, 0.0, 0.0];
        let r = [-1.0_f64, 0.0, 0.0, 0.0];
        // |dot| = 1 → acos(1) = 0, same orientation
        let angle = quaternion_angle(q, r);
        assert!(angle.abs() < 1e-10);
    }

    // ── Distance helper ──────────────────────────────────────────────────────

    #[test]
    fn test_dist3() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 4.0, 0.0];
        assert!((dist3(a, b) - 5.0).abs() < 1e-10);
    }
}
