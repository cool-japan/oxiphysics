// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! DNA/RNA molecular simulation: coarse-grained models, base stacking, elasticity.
//!
//! Implements coarse-grained DNA/RNA models including:
//! - Watson-Crick base pairing with hydrogen bond energetics
//! - Base stacking (π-π interactions)
//! - oxDNA coarse-grained potential
//! - Kirchhoff elastic rod model for DNA mechanics
//! - Worm-like chain (WLC) force-extension model
//! - DNA origami scaffold routing
//! - Melting temperature prediction (nearest-neighbor thermodynamics)

/// Physical constants used in DNA simulation.
pub mod constants {
    /// Boltzmann constant in kcal/(mol·K).
    pub const KB_KCAL: f64 = 0.001987;
    /// Gas constant in kcal/(mol·K).
    pub const R_KCAL: f64 = 0.001987;
    /// Avogadro's number.
    pub const AVOGADRO: f64 = 6.022e23;
    /// Rise per base pair in Angstroms.
    pub const RISE_PER_BP: f64 = 3.4;
    /// Twist per base pair in degrees.
    pub const TWIST_PER_BP: f64 = 36.0;
    /// DNA bending persistence length in nm.
    pub const PERSISTENCE_LENGTH_NM: f64 = 50.0;
    /// DNA twist persistence length in nm.
    pub const TWIST_PERSISTENCE_LENGTH_NM: f64 = 75.0;
    /// DNA stretch modulus in pN.
    pub const STRETCH_MODULUS_PN: f64 = 1000.0;
    /// Room temperature in Kelvin.
    pub const ROOM_TEMP_K: f64 = 298.15;
}

/// Type of nucleotide base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NucleotideType {
    /// Adenine (DNA and RNA).
    A,
    /// Thymine (DNA only).
    T,
    /// Guanine (DNA and RNA).
    G,
    /// Cytosine (DNA and RNA).
    C,
    /// Uracil (RNA only).
    U,
}

impl NucleotideType {
    /// Returns true if this nucleotide is a purine (A or G).
    pub fn is_purine(&self) -> bool {
        matches!(self, NucleotideType::A | NucleotideType::G)
    }

    /// Returns true if this nucleotide is a pyrimidine (T, C, or U).
    pub fn is_pyrimidine(&self) -> bool {
        matches!(
            self,
            NucleotideType::T | NucleotideType::C | NucleotideType::U
        )
    }

    /// Returns the Watson-Crick complement of this base.
    /// For RNA context, use U as complement of A.
    pub fn complement_dna(&self) -> NucleotideType {
        match self {
            NucleotideType::A => NucleotideType::T,
            NucleotideType::T => NucleotideType::A,
            NucleotideType::G => NucleotideType::C,
            NucleotideType::C => NucleotideType::G,
            NucleotideType::U => NucleotideType::A,
        }
    }

    /// Returns the GC content contribution (1 for G or C, 0 otherwise).
    pub fn gc_weight(&self) -> f64 {
        match self {
            NucleotideType::G | NucleotideType::C => 1.0,
            _ => 0.0,
        }
    }

    /// Returns a single-character representation.
    pub fn as_char(&self) -> char {
        match self {
            NucleotideType::A => 'A',
            NucleotideType::T => 'T',
            NucleotideType::G => 'G',
            NucleotideType::C => 'C',
            NucleotideType::U => 'U',
        }
    }
}

/// Coarse-grained DNA base with position and sugar-phosphate geometry.
#[derive(Debug, Clone)]
pub struct DnaBase {
    /// 3D position of the base center in Angstroms \[x, y, z\].
    pub position: [f64; 3],
    /// Type of this nucleotide base.
    pub nucleotide_type: NucleotideType,
    /// Position of the sugar C1' atom relative to base center.
    pub sugar_c1_offset: [f64; 3],
    /// Position of the phosphorus atom relative to base center.
    pub phosphorus_offset: [f64; 3],
    /// Unit vector pointing from sugar to base (glycosidic bond direction).
    pub glycosidic_direction: [f64; 3],
    /// Index along the DNA strand (0-based).
    pub strand_index: usize,
}

impl DnaBase {
    /// Creates a new DnaBase at the given position with default geometry.
    pub fn new(position: [f64; 3], nucleotide_type: NucleotideType, strand_index: usize) -> Self {
        Self {
            position,
            nucleotide_type,
            sugar_c1_offset: [0.0, -2.5, 0.0],
            phosphorus_offset: [0.0, -4.0, -1.0],
            glycosidic_direction: [0.0, -1.0, 0.0],
            strand_index,
        }
    }

    /// Returns the absolute position of the sugar C1' atom.
    pub fn sugar_c1_position(&self) -> [f64; 3] {
        [
            self.position[0] + self.sugar_c1_offset[0],
            self.position[1] + self.sugar_c1_offset[1],
            self.position[2] + self.sugar_c1_offset[2],
        ]
    }

    /// Returns the absolute position of the phosphorus atom.
    pub fn phosphorus_position(&self) -> [f64; 3] {
        [
            self.position[0] + self.phosphorus_offset[0],
            self.position[1] + self.phosphorus_offset[1],
            self.position[2] + self.phosphorus_offset[2],
        ]
    }

    /// Returns the molecular mass of this base in Da.
    pub fn mass_da(&self) -> f64 {
        match self.nucleotide_type {
            NucleotideType::A => 313.2,
            NucleotideType::T => 304.2,
            NucleotideType::G => 329.2,
            NucleotideType::C => 289.2,
            NucleotideType::U => 306.2,
        }
    }
}

/// Watson-Crick base pair geometry and hydrogen bond energetics.
#[derive(Debug, Clone)]
pub struct WatsonCrickPairing {
    /// First base in the pair.
    pub base1: DnaBase,
    /// Second base in the pair (complementary strand).
    pub base2: DnaBase,
    /// Number of hydrogen bonds (2 for AT, 3 for GC).
    pub num_hbonds: u32,
    /// Hydrogen bond energy in kcal/mol (negative = stabilizing).
    pub hbond_energy: f64,
    /// Distance between base centers in Angstroms.
    pub pair_distance: f64,
    /// Propeller twist angle in degrees.
    pub propeller_twist: f64,
    /// Buckle angle in degrees.
    pub buckle: f64,
}

impl WatsonCrickPairing {
    /// Creates a Watson-Crick base pair from two DnaBase instances.
    ///
    /// Returns None if the bases are not complementary.
    pub fn new(base1: DnaBase, base2: DnaBase) -> Option<Self> {
        let (num_hbonds, hbond_energy) = match (&base1.nucleotide_type, &base2.nucleotide_type) {
            (NucleotideType::A, NucleotideType::T)
            | (NucleotideType::T, NucleotideType::A)
            | (NucleotideType::A, NucleotideType::U)
            | (NucleotideType::U, NucleotideType::A) => (2, -7.0),
            (NucleotideType::G, NucleotideType::C) | (NucleotideType::C, NucleotideType::G) => {
                (3, -11.0)
            }
            _ => return None,
        };

        let dx = base2.position[0] - base1.position[0];
        let dy = base2.position[1] - base1.position[1];
        let dz = base2.position[2] - base1.position[2];
        let pair_distance = (dx * dx + dy * dy + dz * dz).sqrt();

        Some(Self {
            base1,
            base2,
            num_hbonds,
            hbond_energy,
            pair_distance,
            propeller_twist: -11.4,
            buckle: 0.0,
        })
    }

    /// Returns true if this is a GC base pair.
    pub fn is_gc_pair(&self) -> bool {
        self.num_hbonds == 3
    }

    /// Returns true if this is an AT/AU base pair.
    pub fn is_at_pair(&self) -> bool {
        self.num_hbonds == 2
    }

    /// Returns the total hydrogen bond stabilization energy in kcal/mol.
    pub fn total_hbond_energy(&self) -> f64 {
        self.hbond_energy
    }

    /// Computes the opening angle of the base pair due to thermal fluctuations (approximate).
    pub fn thermal_opening_angle(&self, temperature_k: f64) -> f64 {
        let kt = constants::KB_KCAL * temperature_k;
        // Simple harmonic approximation: <θ²> ~ kT / k_angle
        let k_angle = self.num_hbonds as f64 * 0.5; // kcal/(mol·rad²) approx
        (kt / k_angle).sqrt().to_degrees()
    }
}

/// π-π base stacking energy between adjacent bases.
#[derive(Debug, Clone)]
pub struct BaseStacking {
    /// Stacking energy in kcal/mol (negative = stabilizing).
    pub stacking_energy: f64,
    /// Rise along helix axis in Angstroms (≈ 3.4 Å).
    pub rise_angstrom: f64,
    /// Twist angle between adjacent base pairs in degrees (≈ 36°).
    pub twist_degrees: f64,
    /// Tilt angle in degrees.
    pub tilt_degrees: f64,
    /// Roll angle in degrees.
    pub roll_degrees: f64,
    /// Types of the two stacked bases.
    pub base_types: (NucleotideType, NucleotideType),
}

impl BaseStacking {
    /// Creates a BaseStacking interaction between two adjacent bases.
    ///
    /// The stacking energy depends on the specific dinucleotide step.
    pub fn new(base_above: NucleotideType, base_below: NucleotideType) -> Self {
        // Nearest-neighbor stacking energies (approximate, kcal/mol)
        let stacking_energy = Self::dinucleotide_stacking_energy(base_above, base_below);
        Self {
            stacking_energy,
            rise_angstrom: constants::RISE_PER_BP,
            twist_degrees: constants::TWIST_PER_BP,
            tilt_degrees: 0.0,
            roll_degrees: 0.0,
            base_types: (base_above, base_below),
        }
    }

    /// Returns the nearest-neighbor stacking energy for a dinucleotide step.
    pub fn dinucleotide_stacking_energy(top: NucleotideType, bottom: NucleotideType) -> f64 {
        // Simplified nearest-neighbor stacking energies (kcal/mol)
        match (top, bottom) {
            (NucleotideType::A, NucleotideType::A) | (NucleotideType::T, NucleotideType::T) => -1.0,
            (NucleotideType::A, NucleotideType::T) => -0.88,
            (NucleotideType::T, NucleotideType::A) => -0.58,
            (NucleotideType::A, NucleotideType::G) | (NucleotideType::C, NucleotideType::T) => {
                -1.28
            }
            (NucleotideType::G, NucleotideType::A) | (NucleotideType::T, NucleotideType::C) => {
                -1.30
            }
            (NucleotideType::C, NucleotideType::A) | (NucleotideType::T, NucleotideType::G) => {
                -1.45
            }
            (NucleotideType::A, NucleotideType::C) | (NucleotideType::G, NucleotideType::T) => {
                -1.44
            }
            (NucleotideType::G, NucleotideType::G) | (NucleotideType::C, NucleotideType::C) => {
                -1.84
            }
            (NucleotideType::G, NucleotideType::C) => -2.24,
            (NucleotideType::C, NucleotideType::G) => -2.17,
            _ => -1.0,
        }
    }

    /// Returns the twist angle in radians.
    pub fn twist_radians(&self) -> f64 {
        self.twist_degrees.to_radians()
    }

    /// Returns the rise in nanometers.
    pub fn rise_nm(&self) -> f64 {
        self.rise_angstrom * 0.1
    }

    /// Computes the stacking overlap area (approximate, Å²).
    pub fn overlap_area(&self) -> f64 {
        // Purines overlap more than pyrimidines
        let area_top = if self.base_types.0.is_purine() {
            35.0
        } else {
            25.0
        };
        let area_bottom = if self.base_types.1.is_purine() {
            35.0
        } else {
            25.0
        };
        0.5 * (area_top + area_bottom) * (1.0 - self.twist_degrees.abs() / 90.0).max(0.0)
    }
}

/// oxDNA coarse-grained potential for DNA simulation.
///
/// Implements the oxDNA model (Ouldridge *et al.* 2011) with excluded volume,
/// backbone connectivity, stacking, and hydrogen bonding terms.
#[derive(Debug, Clone)]
pub struct OxDnaPotential {
    /// Temperature in Kelvin.
    pub temperature: f64,
    /// Excluded volume energy scale (kcal/mol).
    pub eps_exc: f64,
    /// Backbone spring constant (kcal/mol/Å²).
    pub k_backbone: f64,
    /// Equilibrium backbone distance in Angstroms.
    pub r0_backbone: f64,
    /// Stacking energy scale (kcal/mol).
    pub eps_stack: f64,
    /// Hydrogen bonding energy scale (kcal/mol).
    pub eps_hbond: f64,
}

impl OxDnaPotential {
    /// Creates a new OxDnaPotential with default oxDNA parameters.
    pub fn new(temperature: f64) -> Self {
        Self {
            temperature,
            eps_exc: 2.0,
            k_backbone: 2.0,
            r0_backbone: 6.5,
            eps_stack: 1.5,
            eps_hbond: 1.0,
        }
    }

    /// Computes the excluded volume energy between two sites at distance r (Å).
    pub fn excluded_volume_energy(&self, r: f64) -> f64 {
        let sigma = 3.4;
        if r < sigma {
            self.eps_exc * (sigma / r).powi(12)
        } else {
            0.0
        }
    }

    /// Computes the backbone connectivity energy for a bond of length r (Å).
    pub fn backbone_energy(&self, r: f64) -> f64 {
        let dr = r - self.r0_backbone;
        0.5 * self.k_backbone * dr * dr
    }

    /// Computes the stacking potential energy as a function of distance and twist.
    pub fn stacking_energy(&self, r: f64, theta: f64) -> f64 {
        let r0 = 3.4;
        let theta0 = constants::TWIST_PER_BP.to_radians();
        let r_term = (-0.5 * ((r - r0) / 0.5).powi(2)).exp();
        let theta_term = ((theta - theta0).cos() + 1.0) * 0.5;
        -self.eps_stack * r_term * theta_term
    }

    /// Computes the hydrogen bonding energy for a base pair at distance r (Å).
    pub fn hbond_energy(&self, r: f64, base1: NucleotideType, base2: NucleotideType) -> f64 {
        let is_pair = matches!(
            (base1, base2),
            (NucleotideType::A, NucleotideType::T)
                | (NucleotideType::T, NucleotideType::A)
                | (NucleotideType::G, NucleotideType::C)
                | (NucleotideType::C, NucleotideType::G)
        );
        if !is_pair {
            return 0.0;
        }
        let r0 = 5.86;
        let scale = if matches!(
            (base1, base2),
            (NucleotideType::G, NucleotideType::C) | (NucleotideType::C, NucleotideType::G)
        ) {
            1.5
        } else {
            1.0
        };
        let r_term = (-0.5 * ((r - r0) / 0.4).powi(2)).exp();
        -self.eps_hbond * scale * r_term
    }

    /// Returns the total energy for a system of N nucleotides (simplified).
    pub fn total_energy(&self, positions: &[[f64; 3]], types: &[NucleotideType]) -> f64 {
        let n = positions.len().min(types.len());
        let mut energy = 0.0;
        // Backbone connectivity
        for i in 0..n.saturating_sub(1) {
            let dx = positions[i + 1][0] - positions[i][0];
            let dy = positions[i + 1][1] - positions[i][1];
            let dz = positions[i + 1][2] - positions[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            energy += self.backbone_energy(r);
        }
        // Pairwise excluded volume (simplified, only nearest neighbors)
        for i in 0..n {
            for j in (i + 2)..n {
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-10);
                energy += self.excluded_volume_energy(r);
            }
        }
        energy
    }
}

/// DNA modeled as an elastic rod using the Kirchhoff rod theory.
///
/// The Kirchhoff rod model captures bending, twisting, and stretching
/// deformations of DNA with experimentally measured elastic constants.
#[derive(Debug, Clone)]
pub struct KirchhoffRod {
    /// Contour length of the rod in nm.
    pub contour_length_nm: f64,
    /// Bending persistence length L_p in nm (≈ 50 nm for dsDNA).
    pub persistence_length_nm: f64,
    /// Twist persistence length C in nm (≈ 75 nm for dsDNA).
    pub twist_persistence_length_nm: f64,
    /// Stretch modulus S in pN (≈ 1000 pN for dsDNA).
    pub stretch_modulus_pn: f64,
    /// Temperature in Kelvin.
    pub temperature_k: f64,
}

impl KirchhoffRod {
    /// Creates a KirchhoffRod with standard dsDNA elastic constants.
    pub fn new(contour_length_nm: f64, temperature_k: f64) -> Self {
        Self {
            contour_length_nm,
            persistence_length_nm: constants::PERSISTENCE_LENGTH_NM,
            twist_persistence_length_nm: constants::TWIST_PERSISTENCE_LENGTH_NM,
            stretch_modulus_pn: constants::STRETCH_MODULUS_PN,
            temperature_k,
        }
    }

    /// Returns kT in pN·nm units.
    pub fn kt_pn_nm(&self) -> f64 {
        // kT at given temperature in pN·nm (1 kcal/mol = 6.9477 pN·nm)
        constants::KB_KCAL * self.temperature_k * 6947.7 / 1000.0
    }

    /// Computes the bending stiffness A = kT * L_p (pN·nm²).
    pub fn bending_stiffness(&self) -> f64 {
        self.kt_pn_nm() * self.persistence_length_nm
    }

    /// Computes the torsional stiffness C_torsion = kT * C (pN·nm²).
    pub fn torsional_stiffness(&self) -> f64 {
        self.kt_pn_nm() * self.twist_persistence_length_nm
    }

    /// Computes the bending energy for a uniform curvature κ (1/nm) over length L.
    pub fn bending_energy(&self, curvature_per_nm: f64) -> f64 {
        0.5 * self.bending_stiffness()
            * curvature_per_nm
            * curvature_per_nm
            * self.contour_length_nm
    }

    /// Computes the twist energy for excess linking number ΔLk.
    pub fn twist_energy(&self, delta_lk: f64) -> f64 {
        // E_twist = (kT * C * 4π² * ΔLk²) / (2 * L)
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        self.kt_pn_nm() * self.twist_persistence_length_nm * 4.0 * pi2 * delta_lk * delta_lk
            / (2.0 * self.contour_length_nm)
    }

    /// Computes the stretch energy for an extension δL (nm).
    pub fn stretch_energy(&self, delta_l_nm: f64) -> f64 {
        0.5 * self.stretch_modulus_pn * delta_l_nm * delta_l_nm / self.contour_length_nm
    }

    /// Computes the mean-square end-to-end distance `R²` using wormlike chain formula.
    pub fn mean_square_end_to_end(&self) -> f64 {
        let l = self.contour_length_nm;
        let lp = self.persistence_length_nm;
        2.0 * lp * l * (1.0 - (lp / l) * (1.0 - (-l / lp).exp()))
    }
}

/// Worm-Like Chain (WLC) model for DNA force-extension.
///
/// Implements the interpolation formula of Marko and Siggia (1995).
#[derive(Debug, Clone)]
pub struct WormLikeChain {
    /// Contour length L in nm.
    pub contour_length_nm: f64,
    /// Persistence length L_p in nm.
    pub persistence_length_nm: f64,
    /// Temperature in Kelvin.
    pub temperature_k: f64,
    /// Stretch modulus S in pN (for extensible WLC).
    pub stretch_modulus_pn: f64,
}

impl WormLikeChain {
    /// Creates a WormLikeChain with standard dsDNA parameters.
    pub fn new(contour_length_nm: f64, temperature_k: f64) -> Self {
        Self {
            contour_length_nm,
            persistence_length_nm: constants::PERSISTENCE_LENGTH_NM,
            temperature_k,
            stretch_modulus_pn: constants::STRETCH_MODULUS_PN,
        }
    }

    /// Returns kT in pN·nm.
    pub fn kt_pn_nm(&self) -> f64 {
        constants::KB_KCAL * self.temperature_k * 6947.7 / 1000.0
    }

    /// Computes the force in pN for a given fractional extension x/L (0 < x/L < 1).
    ///
    /// Uses the Marko-Siggia interpolation formula:
    /// F = (kT/L_p) * (1/(4*(1-x/L)²) - 1/4 + x/L)
    pub fn force_pn(&self, extension_nm: f64) -> f64 {
        let x_over_l = (extension_nm / self.contour_length_nm).clamp(0.0, 0.999);
        let kt = self.kt_pn_nm();
        let lp = self.persistence_length_nm;
        (kt / lp) * (1.0 / (4.0 * (1.0 - x_over_l).powi(2)) - 0.25 + x_over_l)
    }

    /// Computes extension in nm for a given force in pN by numerical inversion.
    ///
    /// Uses bisection search on the WLC force-extension relation.
    pub fn extension_nm(&self, force_pn: f64) -> f64 {
        if force_pn <= 0.0 {
            return 0.0;
        }
        let l = self.contour_length_nm;
        let mut lo = 0.0_f64;
        let mut hi = 0.9999 * l;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.force_pn(mid) < force_pn {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Computes the fractional extension x/L for a given force in pN.
    pub fn fractional_extension(&self, force_pn: f64) -> f64 {
        self.extension_nm(force_pn) / self.contour_length_nm
    }

    /// Computes the effective spring constant at a given extension (pN/nm).
    pub fn spring_constant_at_extension(&self, extension_nm: f64) -> f64 {
        let x_over_l = (extension_nm / self.contour_length_nm).clamp(0.001, 0.999);
        let kt = self.kt_pn_nm();
        let lp = self.persistence_length_nm;
        // dF/dx = (kT/Lp*L) * (1/(2*(1-x/L)^3) + 1)
        (kt / (lp * self.contour_length_nm)) * (0.5 / (1.0 - x_over_l).powi(3) + 1.0)
    }
}

/// DNA origami design with scaffold and staple strand routing.
///
/// Captures the logical structure of a DNA origami nanostructure.
#[derive(Debug, Clone)]
pub struct DnaOrigami {
    /// Scaffold strand (M13mp18 or custom), as nucleotide sequence.
    pub scaffold: Vec<NucleotideType>,
    /// Staple strands, each as a nucleotide sequence.
    pub staples: Vec<Vec<NucleotideType>>,
    /// Helix count in the origami lattice.
    pub num_helices: usize,
    /// Target shape descriptor (e.g., "rectangle", "tube", "custom").
    pub shape_descriptor: String,
    /// Estimated contour length of scaffold in nm.
    pub scaffold_length_nm: f64,
}

impl DnaOrigami {
    /// Creates a new DnaOrigami with the given scaffold.
    pub fn new(scaffold: Vec<NucleotideType>, shape_descriptor: &str) -> Self {
        let n = scaffold.len();
        let length_nm = n as f64 * constants::RISE_PER_BP * 0.1;
        Self {
            scaffold,
            staples: Vec::new(),
            num_helices: 1,
            shape_descriptor: shape_descriptor.to_string(),
            scaffold_length_nm: length_nm,
        }
    }

    /// Adds a staple strand to the origami design.
    pub fn add_staple(&mut self, staple: Vec<NucleotideType>) {
        self.staples.push(staple);
    }

    /// Returns the number of scaffold nucleotides.
    pub fn scaffold_length(&self) -> usize {
        self.scaffold.len()
    }

    /// Returns the number of staple strands.
    pub fn num_staples(&self) -> usize {
        self.staples.len()
    }

    /// Computes the GC content fraction of the scaffold.
    pub fn scaffold_gc_content(&self) -> f64 {
        if self.scaffold.is_empty() {
            return 0.0;
        }
        let gc_count = self.scaffold.iter().filter(|b| b.gc_weight() > 0.5).count();
        gc_count as f64 / self.scaffold.len() as f64
    }

    /// Estimates total scaffold coverage by staples (fraction of scaffold paired).
    pub fn staple_coverage(&self) -> f64 {
        let total_staple_nt: usize = self.staples.iter().map(|s| s.len()).sum();
        if self.scaffold.is_empty() {
            return 0.0;
        }
        (total_staple_nt as f64 / self.scaffold.len() as f64).min(1.0)
    }
}

/// Torque-angle relationship for twisted DNA.
///
/// Implements the linear twist-torque relationship and describes
/// the transition from B-DNA to over/under-wound forms.
#[derive(Debug, Clone)]
pub struct TwistingForce {
    /// Torsional stiffness C in pN·nm² (= kT * C_torsional_persistence).
    pub torsional_stiffness_pn_nm2: f64,
    /// Contour length of the DNA segment in nm.
    pub length_nm: f64,
    /// Equilibrium twist density (rad/bp × bp/nm).
    pub twist_density_rad_per_nm: f64,
}

impl TwistingForce {
    /// Creates a TwistingForce with default dsDNA parameters.
    pub fn new(length_nm: f64, temperature_k: f64) -> Self {
        let kt = constants::KB_KCAL * temperature_k * 6947.7 / 1000.0;
        let torsional_stiffness = kt * constants::TWIST_PERSISTENCE_LENGTH_NM;
        // Twist density: 36°/bp = π/5 rad/bp; 10 bp/turn; rise = 3.4 Å/bp = 0.34 nm/bp
        let twist_density_rad_per_nm =
            constants::TWIST_PER_BP.to_radians() / (constants::RISE_PER_BP * 0.1);
        Self {
            torsional_stiffness_pn_nm2: torsional_stiffness,
            length_nm,
            twist_density_rad_per_nm,
        }
    }

    /// Computes the torque (pN·nm) for a given excess linking number ΔLk.
    pub fn torque_pn_nm(&self, delta_lk: f64) -> f64 {
        // τ = C * 2π * ΔLk / L
        self.torsional_stiffness_pn_nm2 * 2.0 * std::f64::consts::PI * delta_lk / self.length_nm
    }

    /// Computes ΔLk for a given torque (pN·nm).
    pub fn delta_lk_from_torque(&self, torque_pn_nm: f64) -> f64 {
        torque_pn_nm * self.length_nm
            / (self.torsional_stiffness_pn_nm2 * 2.0 * std::f64::consts::PI)
    }

    /// Computes the twist energy for a given ΔLk in kcal/mol.
    pub fn twist_energy_kcal(&self, delta_lk: f64) -> f64 {
        let e_pn_nm =
            0.5 * self.torsional_stiffness_pn_nm2 * (2.0 * std::f64::consts::PI * delta_lk).powi(2)
                / self.length_nm;
        e_pn_nm / 6947.7 * 1000.0 // pN·nm → kcal/mol
    }

    /// Returns the equilibrium twist angle per nm (rad/nm).
    pub fn equilibrium_twist_density(&self) -> f64 {
        self.twist_density_rad_per_nm
    }
}

/// Nearest-neighbor thermodynamics for DNA melting temperature.
///
/// Implements the SantaLucia (1998) unified nearest-neighbor parameters
/// for dsDNA stability and melting temperature prediction.
#[derive(Debug, Clone)]
pub struct MeltingTemperature {
    /// Duplex sequence (5'→3' strand).
    pub sequence: Vec<NucleotideType>,
    /// Total strand concentration C_T in molar.
    pub strand_concentration: f64,
    /// Salt concentration \[Na+\] in molar.
    pub salt_concentration: f64,
    /// Enthalpy ΔH in kcal/mol.
    pub delta_h: f64,
    /// Entropy ΔS in cal/(mol·K).
    pub delta_s: f64,
}

impl MeltingTemperature {
    /// Nearest-neighbor enthalpy (kcal/mol) and entropy (cal/mol/K) parameters.
    /// Based on SantaLucia 1998.
    fn nn_params(b1: NucleotideType, b2: NucleotideType) -> (f64, f64) {
        match (b1, b2) {
            (NucleotideType::A, NucleotideType::A) | (NucleotideType::T, NucleotideType::T) => {
                (-7.9, -22.2)
            }
            (NucleotideType::A, NucleotideType::T) => (-7.2, -20.4),
            (NucleotideType::T, NucleotideType::A) => (-7.2, -21.3),
            (NucleotideType::C, NucleotideType::A) | (NucleotideType::T, NucleotideType::G) => {
                (-8.5, -22.7)
            }
            (NucleotideType::G, NucleotideType::T) | (NucleotideType::A, NucleotideType::C) => {
                (-8.4, -22.4)
            }
            (NucleotideType::C, NucleotideType::T) | (NucleotideType::A, NucleotideType::G) => {
                (-7.8, -21.0)
            }
            (NucleotideType::G, NucleotideType::A) | (NucleotideType::T, NucleotideType::C) => {
                (-7.8, -21.0)
            }
            (NucleotideType::G, NucleotideType::G) | (NucleotideType::C, NucleotideType::C) => {
                (-8.0, -19.9)
            }
            (NucleotideType::G, NucleotideType::C) => (-9.8, -24.4),
            (NucleotideType::C, NucleotideType::G) => (-10.6, -27.2),
            _ => (-8.0, -22.0),
        }
    }

    /// Creates a MeltingTemperature calculator for the given sequence.
    pub fn new(sequence: Vec<NucleotideType>, strand_concentration_m: f64) -> Self {
        let mut delta_h = 0.0_f64;
        let mut delta_s = 0.0_f64;
        // Initiation penalty
        delta_h += 0.2;
        delta_s += -5.7;
        // Nearest-neighbor sum
        for i in 0..sequence.len().saturating_sub(1) {
            let (dh, ds) = Self::nn_params(sequence[i], sequence[i + 1]);
            delta_h += dh;
            delta_s += ds;
        }
        Self {
            sequence,
            strand_concentration: strand_concentration_m,
            salt_concentration: 1.0,
            delta_h,
            delta_s,
        }
    }

    /// Computes the melting temperature T_m in Kelvin.
    ///
    /// Uses: T_m = ΔH / (ΔS + R·ln(C_T/4))
    pub fn melting_temperature_k(&self) -> f64 {
        let r_cal = 1.987; // cal/(mol·K)
        let ds = self.delta_s; // cal/(mol·K)
        let dh = self.delta_h * 1000.0; // kcal/mol → cal/mol
        let ct = self.strand_concentration;
        dh / (ds + r_cal * (ct / 4.0).ln())
    }

    /// Computes the melting temperature T_m in Celsius.
    pub fn melting_temperature_c(&self) -> f64 {
        self.melting_temperature_k() - 273.15
    }

    /// Returns the GC content fraction of this sequence.
    pub fn gc_content(&self) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let gc = self.sequence.iter().filter(|b| b.gc_weight() > 0.5).count();
        gc as f64 / self.sequence.len() as f64
    }

    /// Applies a simplified salt correction (Owen-Richards).
    pub fn salt_corrected_tm_k(&self) -> f64 {
        let tm = self.melting_temperature_k();
        // ΔTm ≈ 16.6 * log10([Na+]) — simplified approximation
        tm + 16.6 * (self.salt_concentration.log10())
    }
}

/// Computes the GC content fraction of a nucleotide sequence.
pub fn gc_content(sequence: &[NucleotideType]) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }
    let gc = sequence.iter().filter(|b| b.gc_weight() > 0.5).count();
    gc as f64 / sequence.len() as f64
}

/// Generates a B-DNA helix with n base pairs starting at origin.
///
/// Returns a vector of (sense, antisense) base pairs along the helix axis.
pub fn generate_bform_dna(sequence: &[NucleotideType]) -> Vec<(DnaBase, DnaBase)> {
    let rise = constants::RISE_PER_BP;
    let twist_deg = constants::TWIST_PER_BP;
    let radius = 10.0_f64; // Å, approx C1'-C1' distance / 2
    sequence
        .iter()
        .enumerate()
        .map(|(i, &nt)| {
            let z = i as f64 * rise;
            let angle = (i as f64 * twist_deg).to_radians();
            let x = radius * angle.cos();
            let y = radius * angle.sin();
            let sense = DnaBase::new([x, y, z], nt, i);
            let comp_nt = nt.complement_dna();
            let anti = DnaBase::new([-x, -y, z], comp_nt, i);
            (sense, anti)
        })
        .collect()
}

/// Computes the number of base pairs in a sequence that form Watson-Crick pairs.
pub fn count_watson_crick_pairs(seq1: &[NucleotideType], seq2: &[NucleotideType]) -> usize {
    seq1.iter()
        .zip(seq2.iter())
        .filter(|(a, b)| {
            matches!(
                (a, b),
                (NucleotideType::A, NucleotideType::T)
                    | (NucleotideType::T, NucleotideType::A)
                    | (NucleotideType::G, NucleotideType::C)
                    | (NucleotideType::C, NucleotideType::G)
                    | (NucleotideType::A, NucleotideType::U)
                    | (NucleotideType::U, NucleotideType::A)
            )
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- NucleotideType tests ----

    #[test]
    fn test_nucleotide_purine_pyrimidine() {
        assert!(NucleotideType::A.is_purine());
        assert!(NucleotideType::G.is_purine());
        assert!(!NucleotideType::T.is_purine());
        assert!(NucleotideType::T.is_pyrimidine());
        assert!(NucleotideType::C.is_pyrimidine());
        assert!(NucleotideType::U.is_pyrimidine());
    }

    #[test]
    fn test_nucleotide_complement() {
        assert_eq!(NucleotideType::A.complement_dna(), NucleotideType::T);
        assert_eq!(NucleotideType::T.complement_dna(), NucleotideType::A);
        assert_eq!(NucleotideType::G.complement_dna(), NucleotideType::C);
        assert_eq!(NucleotideType::C.complement_dna(), NucleotideType::G);
    }

    #[test]
    fn test_nucleotide_gc_weight() {
        assert_eq!(NucleotideType::G.gc_weight(), 1.0);
        assert_eq!(NucleotideType::C.gc_weight(), 1.0);
        assert_eq!(NucleotideType::A.gc_weight(), 0.0);
        assert_eq!(NucleotideType::T.gc_weight(), 0.0);
    }

    #[test]
    fn test_nucleotide_as_char() {
        assert_eq!(NucleotideType::A.as_char(), 'A');
        assert_eq!(NucleotideType::U.as_char(), 'U');
    }

    // ---- DnaBase tests ----

    #[test]
    fn test_dna_base_construction() {
        let base = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::A, 0);
        assert_eq!(base.nucleotide_type, NucleotideType::A);
        assert_eq!(base.strand_index, 0);
    }

    #[test]
    fn test_dna_base_mass() {
        let base_a = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::A, 0);
        let base_g = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::G, 0);
        // Purine should be heavier
        assert!(base_g.mass_da() > base_a.mass_da());
    }

    #[test]
    fn test_dna_base_sugar_position() {
        let base = DnaBase::new([5.0, 5.0, 5.0], NucleotideType::C, 0);
        let sugar = base.sugar_c1_position();
        // Should be offset from position
        assert!((sugar[1] - 2.5).abs() < 1e-10);
    }

    // ---- WatsonCrickPairing tests ----

    #[test]
    fn test_watson_crick_at_pairing() {
        let b1 = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::A, 0);
        let b2 = DnaBase::new([10.0, 0.0, 0.0], NucleotideType::T, 1);
        let pair = WatsonCrickPairing::new(b1, b2).unwrap();
        assert_eq!(pair.num_hbonds, 2);
        assert!(pair.hbond_energy < 0.0);
        assert!(pair.is_at_pair());
    }

    #[test]
    fn test_watson_crick_gc_pairing() {
        let b1 = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::G, 0);
        let b2 = DnaBase::new([10.0, 0.0, 0.0], NucleotideType::C, 1);
        let pair = WatsonCrickPairing::new(b1, b2).unwrap();
        assert_eq!(pair.num_hbonds, 3);
        assert!(pair.hbond_energy < 0.0);
        assert!(pair.is_gc_pair());
    }

    #[test]
    fn test_watson_crick_gc_stronger_than_at() {
        let b1_at = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::A, 0);
        let b2_at = DnaBase::new([10.0, 0.0, 0.0], NucleotideType::T, 1);
        let pair_at = WatsonCrickPairing::new(b1_at, b2_at).unwrap();

        let b1_gc = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::G, 0);
        let b2_gc = DnaBase::new([10.0, 0.0, 0.0], NucleotideType::C, 1);
        let pair_gc = WatsonCrickPairing::new(b1_gc, b2_gc).unwrap();

        // GC should have more negative (stronger) energy
        assert!(pair_gc.hbond_energy < pair_at.hbond_energy);
    }

    #[test]
    fn test_watson_crick_mismatch_returns_none() {
        let b1 = DnaBase::new([0.0, 0.0, 0.0], NucleotideType::A, 0);
        let b2 = DnaBase::new([10.0, 0.0, 0.0], NucleotideType::G, 1);
        assert!(WatsonCrickPairing::new(b1, b2).is_none());
    }

    // ---- BaseStacking tests ----

    #[test]
    fn test_base_stacking_rise() {
        let stack = BaseStacking::new(NucleotideType::A, NucleotideType::T);
        assert!(
            (stack.rise_angstrom - 3.4).abs() < 1e-10,
            "rise = {:.6}",
            stack.rise_angstrom
        );
    }

    #[test]
    fn test_base_stacking_twist() {
        let stack = BaseStacking::new(NucleotideType::G, NucleotideType::C);
        assert!(
            (stack.twist_degrees - 36.0).abs() < 1e-10,
            "twist = {:.6}",
            stack.twist_degrees
        );
    }

    #[test]
    fn test_base_stacking_energy_negative() {
        for (a, b) in [
            (NucleotideType::A, NucleotideType::A),
            (NucleotideType::G, NucleotideType::C),
            (NucleotideType::C, NucleotideType::G),
            (NucleotideType::A, NucleotideType::T),
        ] {
            let stack = BaseStacking::new(a, b);
            assert!(
                stack.stacking_energy < 0.0,
                "stacking energy should be negative for {:?}/{:?}",
                a,
                b
            );
        }
    }

    #[test]
    fn test_base_stacking_rise_nm() {
        let stack = BaseStacking::new(NucleotideType::A, NucleotideType::A);
        assert!((stack.rise_nm() - 0.34).abs() < 1e-10);
    }

    // ---- OxDnaPotential tests ----

    #[test]
    fn test_oxdna_excluded_volume_repulsive() {
        let pot = OxDnaPotential::new(300.0);
        // Very close: should be repulsive
        let e_close = pot.excluded_volume_energy(1.0);
        let e_far = pot.excluded_volume_energy(5.0);
        assert!(e_close > e_far);
        assert!(e_far == 0.0); // beyond cutoff
    }

    #[test]
    fn test_oxdna_backbone_energy_zero_at_r0() {
        let pot = OxDnaPotential::new(300.0);
        let e = pot.backbone_energy(pot.r0_backbone);
        assert!(e.abs() < 1e-12);
    }

    #[test]
    fn test_oxdna_hbond_energy_gc_stronger() {
        let pot = OxDnaPotential::new(300.0);
        let e_gc = pot.hbond_energy(5.86, NucleotideType::G, NucleotideType::C);
        let e_at = pot.hbond_energy(5.86, NucleotideType::A, NucleotideType::T);
        assert!(e_gc < e_at);
    }

    #[test]
    fn test_oxdna_hbond_mismatch_zero() {
        let pot = OxDnaPotential::new(300.0);
        let e = pot.hbond_energy(5.86, NucleotideType::A, NucleotideType::G);
        assert_eq!(e, 0.0);
    }

    // ---- KirchhoffRod tests ----

    #[test]
    fn test_kirchhoff_rod_persistence_length() {
        let rod = KirchhoffRod::new(100.0, 300.0);
        assert!((rod.persistence_length_nm - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_kirchhoff_rod_twist_persistence_length() {
        let rod = KirchhoffRod::new(100.0, 300.0);
        assert!((rod.twist_persistence_length_nm - 75.0).abs() < 1e-10);
    }

    #[test]
    fn test_kirchhoff_rod_bending_stiffness_positive() {
        let rod = KirchhoffRod::new(100.0, 300.0);
        assert!(rod.bending_stiffness() > 0.0);
    }

    #[test]
    fn test_kirchhoff_rod_bending_energy_zero_for_zero_curvature() {
        let rod = KirchhoffRod::new(100.0, 300.0);
        let e = rod.bending_energy(0.0);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_kirchhoff_rod_stretch_energy_zero_for_zero_delta() {
        let rod = KirchhoffRod::new(100.0, 300.0);
        let e = rod.stretch_energy(0.0);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_kirchhoff_rod_mean_square_end_to_end_positive() {
        let rod = KirchhoffRod::new(100.0, 300.0);
        assert!(rod.mean_square_end_to_end() > 0.0);
    }

    // ---- WormLikeChain tests ----

    #[test]
    fn test_wlc_force_increases_with_extension() {
        let wlc = WormLikeChain::new(1000.0, 300.0);
        let f_low = wlc.force_pn(100.0);
        let f_high = wlc.force_pn(900.0);
        assert!(f_high > f_low);
    }

    #[test]
    fn test_wlc_force_diverges_near_full_extension() {
        let wlc = WormLikeChain::new(1000.0, 300.0);
        // Extension close to contour length: force should be very large
        let f_near = wlc.force_pn(999.0);
        assert!(f_near > 100.0, "force near full extension = {:.6}", f_near);
    }

    #[test]
    fn test_wlc_extension_inverse_of_force() {
        let wlc = WormLikeChain::new(1000.0, 300.0);
        let ext_target = 500.0_f64;
        let f = wlc.force_pn(ext_target);
        let ext_recovered = wlc.extension_nm(f);
        assert!(
            (ext_recovered - ext_target).abs() < 0.5,
            "recovered {:.6} vs target {:.6}",
            ext_recovered,
            ext_target
        );
    }

    #[test]
    fn test_wlc_fractional_extension_bounded() {
        let wlc = WormLikeChain::new(1000.0, 300.0);
        for f in [1.0, 10.0, 100.0] {
            let frac = wlc.fractional_extension(f);
            assert!(
                frac > 0.0 && frac < 1.0,
                "fractional extension out of bounds: {:.6}",
                frac
            );
        }
    }

    #[test]
    fn test_wlc_spring_constant_increases_with_extension() {
        let wlc = WormLikeChain::new(1000.0, 300.0);
        let k_low = wlc.spring_constant_at_extension(100.0);
        let k_high = wlc.spring_constant_at_extension(800.0);
        assert!(k_high > k_low);
    }

    // ---- MeltingTemperature tests ----

    #[test]
    fn test_melting_temp_gc_rich_higher_than_at_rich() {
        let at_seq = vec![
            NucleotideType::A,
            NucleotideType::T,
            NucleotideType::A,
            NucleotideType::T,
            NucleotideType::A,
            NucleotideType::T,
            NucleotideType::A,
            NucleotideType::T,
        ];
        let gc_seq = vec![
            NucleotideType::G,
            NucleotideType::C,
            NucleotideType::G,
            NucleotideType::C,
            NucleotideType::G,
            NucleotideType::C,
            NucleotideType::G,
            NucleotideType::C,
        ];
        let tm_at = MeltingTemperature::new(at_seq, 1e-6).melting_temperature_k();
        let tm_gc = MeltingTemperature::new(gc_seq, 1e-6).melting_temperature_k();
        assert!(
            tm_gc > tm_at,
            "GC Tm = {:.6} should be > AT Tm = {:.6}",
            tm_gc,
            tm_at
        );
    }

    #[test]
    fn test_melting_temp_gc_content_correlation() {
        let at_seq = vec![
            NucleotideType::A,
            NucleotideType::T,
            NucleotideType::A,
            NucleotideType::T,
        ];
        let gc_seq = vec![
            NucleotideType::G,
            NucleotideType::C,
            NucleotideType::G,
            NucleotideType::C,
        ];
        let mt_at = MeltingTemperature::new(at_seq, 1e-6);
        let mt_gc = MeltingTemperature::new(gc_seq, 1e-6);
        assert!(mt_gc.gc_content() > mt_at.gc_content());
    }

    #[test]
    fn test_melting_temp_increases_with_length() {
        let short_seq: Vec<_> = (0..5)
            .map(|i| {
                if i % 2 == 0 {
                    NucleotideType::G
                } else {
                    NucleotideType::C
                }
            })
            .collect();
        let long_seq: Vec<_> = (0..20)
            .map(|i| {
                if i % 2 == 0 {
                    NucleotideType::G
                } else {
                    NucleotideType::C
                }
            })
            .collect();
        let tm_short = MeltingTemperature::new(short_seq, 1e-6).melting_temperature_k();
        let tm_long = MeltingTemperature::new(long_seq, 1e-6).melting_temperature_k();
        assert!(tm_long > tm_short);
    }

    // ---- TwistingForce tests ----

    #[test]
    fn test_twisting_force_torque_linear_in_delta_lk() {
        let tf = TwistingForce::new(100.0, 300.0);
        let tau1 = tf.torque_pn_nm(1.0);
        let tau2 = tf.torque_pn_nm(2.0);
        assert!((tau2 / tau1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_twisting_force_delta_lk_roundtrip() {
        let tf = TwistingForce::new(100.0, 300.0);
        let delta_lk = 3.0;
        let torque = tf.torque_pn_nm(delta_lk);
        let recovered = tf.delta_lk_from_torque(torque);
        assert!((recovered - delta_lk).abs() < 1e-10);
    }

    #[test]
    fn test_twisting_energy_positive() {
        let tf = TwistingForce::new(100.0, 300.0);
        assert!(tf.twist_energy_kcal(1.0) > 0.0);
    }

    // ---- DnaOrigami tests ----

    #[test]
    fn test_dna_origami_scaffold_length() {
        let scaffold: Vec<NucleotideType> = (0..100)
            .map(|i| {
                if i % 4 == 0 {
                    NucleotideType::A
                } else if i % 4 == 1 {
                    NucleotideType::T
                } else if i % 4 == 2 {
                    NucleotideType::G
                } else {
                    NucleotideType::C
                }
            })
            .collect();
        let origami = DnaOrigami::new(scaffold, "rectangle");
        assert_eq!(origami.scaffold_length(), 100);
    }

    #[test]
    fn test_dna_origami_staple_coverage() {
        let scaffold: Vec<NucleotideType> = (0..100).map(|_| NucleotideType::A).collect();
        let mut origami = DnaOrigami::new(scaffold, "test");
        let staple: Vec<NucleotideType> = (0..50).map(|_| NucleotideType::T).collect();
        origami.add_staple(staple);
        let coverage = origami.staple_coverage();
        assert!((coverage - 0.5).abs() < 1e-10);
    }

    // ---- Utility function tests ----

    #[test]
    fn test_gc_content_all_gc() {
        let seq = vec![
            NucleotideType::G,
            NucleotideType::C,
            NucleotideType::G,
            NucleotideType::C,
        ];
        assert!((gc_content(&seq) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_all_at() {
        let seq = vec![
            NucleotideType::A,
            NucleotideType::T,
            NucleotideType::A,
            NucleotideType::T,
        ];
        assert!((gc_content(&seq) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_gc_content_mixed() {
        let seq = vec![
            NucleotideType::G,
            NucleotideType::A,
            NucleotideType::T,
            NucleotideType::C,
        ];
        assert!((gc_content(&seq) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_generate_bform_dna() {
        let seq = vec![
            NucleotideType::A,
            NucleotideType::T,
            NucleotideType::G,
            NucleotideType::C,
        ];
        let helix = generate_bform_dna(&seq);
        assert_eq!(helix.len(), 4);
        // Check rise per base pair
        let rise_0_1 = (helix[1].0.position[2] - helix[0].0.position[2]).abs();
        assert!((rise_0_1 - 3.4).abs() < 1e-10, "rise = {:.6}", rise_0_1);
    }

    #[test]
    fn test_generate_bform_dna_complementarity() {
        let seq = vec![
            NucleotideType::A,
            NucleotideType::G,
            NucleotideType::C,
            NucleotideType::T,
        ];
        let helix = generate_bform_dna(&seq);
        for (sense, anti) in &helix {
            assert_eq!(sense.nucleotide_type.complement_dna(), anti.nucleotide_type);
        }
    }

    #[test]
    fn test_count_watson_crick_pairs() {
        let s1 = vec![NucleotideType::A, NucleotideType::G, NucleotideType::A];
        let s2 = vec![NucleotideType::T, NucleotideType::C, NucleotideType::A];
        // A-T and G-C are valid; A-A is not
        assert_eq!(count_watson_crick_pairs(&s1, &s2), 2);
    }

    #[test]
    fn test_wlc_zero_force_at_zero_extension() {
        let wlc = WormLikeChain::new(1000.0, 300.0);
        let f = wlc.force_pn(0.0);
        assert!(f >= 0.0);
    }

    #[test]
    fn test_kirchhoff_rod_stretch_modulus() {
        let rod = KirchhoffRod::new(100.0, 300.0);
        assert!((rod.stretch_modulus_pn - 1000.0).abs() < 1e-10);
    }
}
