#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tribochemistry molecular dynamics module.
//!
//! Covers:
//! - Tribochemical reactions under mechanochemical activation (Bell model)
//! - Bond breaking and formation under shear stress
//! - Tribopolymer formation from boundary lubricant molecules
//! - Lubricant additive decomposition under contact
//! - ZDDP-like tribofilm formation and growth kinetics
//! - Mechanically-activated reaction pathways
//! - Contact stress effects on reaction rate
//! - Shear-induced phase transformations
//! - Wear particle formation and ejection
//! - Chemical potential under mechanical stress

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
const K_B: f64 = 1.380_649e-23;

/// Planck constant (J s).
const H_PLANCK: f64 = 6.626_070_15e-34;

/// Avogadro's number (mol⁻¹).
const N_AV: f64 = 6.022_140_76e23;

/// Universal gas constant (J mol⁻¹ K⁻¹).
const R_GAS: f64 = 8.314_462_618;

/// Elementary charge (C).
const E_CHARGE: f64 = 1.602_176_634e-19;

/// Atomic mass unit (kg).
const AMU: f64 = 1.660_539_066_6e-27;

/// Angstrom in meters.
const ANGSTROM: f64 = 1.0e-10;

/// Typical covalent bond length for C-C (m).
const BOND_CC: f64 = 1.54e-10;

/// Typical thermal velocity prefactor.
const K_TST: f64 = K_B / H_PLANCK;

// ---------------------------------------------------------------------------
// Vector helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors (a − b).
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Normalize a 3-vector.
#[inline]
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1.0e-30 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}

/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Bell model — bond breaking under stress
// ---------------------------------------------------------------------------

/// Bell model parameters for a single bond type.
///
/// The Bell (1978) model gives the stress-dependent bond lifetime:
/// τ(F) = τ₀ · exp((E_a − F·x‡) / k_B T)
#[derive(Debug, Clone)]
pub struct BellModelParams {
    /// Zero-force activation energy (J).
    pub activation_energy: f64,
    /// Attempt frequency prefactor (Hz).
    pub attempt_freq: f64,
    /// Reaction distance / force sensitivity (m).
    pub reaction_distance: f64,
    /// Bond spring constant (N m⁻¹).
    pub spring_constant: f64,
}

impl BellModelParams {
    /// Default parameters for a C-C covalent bond.
    pub fn carbon_carbon() -> Self {
        BellModelParams {
            activation_energy: 3.6e-19, // ~2.25 eV
            attempt_freq: 1.0e13,
            reaction_distance: 0.5e-10,
            spring_constant: 400.0,
        }
    }

    /// Default parameters for a C-S bond (relevant for ZDDP).
    pub fn carbon_sulfur() -> Self {
        BellModelParams {
            activation_energy: 2.9e-19,
            attempt_freq: 1.0e13,
            reaction_distance: 0.6e-10,
            spring_constant: 280.0,
        }
    }

    /// Default parameters for a P=S bond.
    pub fn phosphorus_sulfur() -> Self {
        BellModelParams {
            activation_energy: 2.4e-19,
            attempt_freq: 8.0e12,
            reaction_distance: 0.7e-10,
            spring_constant: 210.0,
        }
    }

    /// Stress-modified activation energy (J).
    ///
    /// `force` is the applied tensile force (N).
    pub fn effective_activation(&self, force: f64) -> f64 {
        (self.activation_energy - force * self.reaction_distance).max(0.0)
    }

    /// Bond breaking rate constant (s⁻¹) under force `force` (N) at temperature `t` (K).
    pub fn breaking_rate(&self, force: f64, t: f64) -> f64 {
        let e_eff = self.effective_activation(force);
        self.attempt_freq * (-e_eff / (K_B * t)).exp()
    }

    /// Mean bond lifetime (s) under force `force` (N) at temperature `t` (K).
    pub fn mean_lifetime(&self, force: f64, t: f64) -> f64 {
        1.0 / self.breaking_rate(force, t)
    }

    /// Probability that bond has broken after time `dt` (s).
    pub fn breaking_probability(&self, force: f64, t: f64, dt: f64) -> f64 {
        let k = self.breaking_rate(force, t);
        1.0 - (-k * dt).exp()
    }
}

// ---------------------------------------------------------------------------
// Shear stress on reaction rate
// ---------------------------------------------------------------------------

/// Compute the contact normal stress (Pa) for a Hertzian contact.
///
/// `e_star` is the reduced Young's modulus (Pa), `r_star` the reduced radius (m),
/// `load` the normal load (N).
pub fn hertz_contact_pressure(e_star: f64, r_star: f64, load: f64) -> f64 {
    let a3 = 3.0 * load * r_star / (4.0 * e_star);
    let a = a3.cbrt();
    if a < 1.0e-30 {
        return 0.0;
    }
    3.0 * load / (2.0 * PI * a * a)
}

/// Contact half-width for Hertzian contact (m).
pub fn hertz_contact_radius(e_star: f64, r_star: f64, load: f64) -> f64 {
    let a3 = 3.0 * load * r_star / (4.0 * e_star);
    a3.cbrt()
}

/// Shear stress at a tribochemical contact (Pa).
///
/// Converts friction coefficient `mu_friction`, contact pressure `p0` (Pa)
/// to mean contact shear stress.
pub fn contact_shear_stress(mu_friction: f64, p0: f64) -> f64 {
    mu_friction * p0
}

/// Reaction rate enhancement factor from shear stress (dimensionless).
///
/// Uses the mechanochemical coupling factor `gamma` (m³), shear stress `tau` (Pa),
/// and temperature `t` (K).
pub fn shear_rate_enhancement(gamma: f64, tau: f64, t: f64) -> f64 {
    (gamma * tau / (K_B * t)).exp()
}

// ---------------------------------------------------------------------------
// Chemical potential under mechanical stress
// ---------------------------------------------------------------------------

/// Chemical potential modification under hydrostatic stress.
///
/// `mu_0` is the unstressed chemical potential (J mol⁻¹),
/// `sigma_h` is the hydrostatic stress (Pa, positive = compression),
/// `v_molar` is the molar volume (m³ mol⁻¹).
pub fn chemical_potential_stressed(mu_0: f64, sigma_h: f64, v_molar: f64) -> f64 {
    mu_0 - sigma_h * v_molar
}

/// Excess chemical potential due to shear strain energy (J mol⁻¹).
///
/// `gamma_shear` shear strain (dimensionless), `g_modulus` shear modulus (Pa),
/// `v_molar` molar volume (m³ mol⁻¹).
pub fn shear_strain_chemical_potential(gamma_shear: f64, g_modulus: f64, v_molar: f64) -> f64 {
    0.5 * g_modulus * gamma_shear * gamma_shear * v_molar * N_AV
}

// ---------------------------------------------------------------------------
// Atom types for tribochemistry
// ---------------------------------------------------------------------------

/// Atom species relevant to tribochemical systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TriboChem {
    /// Carbon.
    C,
    /// Hydrogen.
    H,
    /// Oxygen.
    O,
    /// Sulfur.
    S,
    /// Phosphorus.
    P,
    /// Zinc.
    Zn,
    /// Iron (steel substrate).
    Fe,
    /// Nitrogen.
    N,
}

impl TriboChem {
    /// Atomic mass (kg) of the species.
    pub fn mass(self) -> f64 {
        let amu_mass: f64 = match self {
            TriboChem::C => 12.011,
            TriboChem::H => 1.008,
            TriboChem::O => 15.999,
            TriboChem::S => 32.06,
            TriboChem::P => 30.974,
            TriboChem::Zn => 65.38,
            TriboChem::Fe => 55.845,
            TriboChem::N => 14.007,
        };
        amu_mass * AMU
    }

    /// Covalent radius (m) of the species.
    pub fn covalent_radius(self) -> f64 {
        match self {
            TriboChem::C => 0.77e-10,
            TriboChem::H => 0.31e-10,
            TriboChem::O => 0.66e-10,
            TriboChem::S => 1.05e-10,
            TriboChem::P => 1.07e-10,
            TriboChem::Zn => 1.22e-10,
            TriboChem::Fe => 1.26e-10,
            TriboChem::N => 0.71e-10,
        }
    }
}

// ---------------------------------------------------------------------------
// Bond
// ---------------------------------------------------------------------------

/// A covalent bond between two atoms in a tribochemical simulation.
#[derive(Debug, Clone)]
pub struct TriboBond {
    /// Index of atom i.
    pub atom_i: usize,
    /// Index of atom j.
    pub atom_j: usize,
    /// Equilibrium bond length (m).
    pub length_eq: f64,
    /// Bond spring constant (N m⁻¹).
    pub k_spring: f64,
    /// Dissociation energy (J).
    pub dissociation_energy: f64,
    /// Whether this bond is currently broken.
    pub broken: bool,
    /// Bell model parameters for this bond.
    pub bell: BellModelParams,
    /// Cumulative force applied to bond (N).
    pub cumulative_force: f64,
}

impl TriboBond {
    /// Create a new C-C bond between atoms `i` and `j`.
    pub fn new_cc(atom_i: usize, atom_j: usize) -> Self {
        TriboBond {
            atom_i,
            atom_j,
            length_eq: BOND_CC,
            k_spring: 400.0,
            dissociation_energy: 3.6e-19,
            broken: false,
            bell: BellModelParams::carbon_carbon(),
            cumulative_force: 0.0,
        }
    }

    /// Create a new C-S bond.
    pub fn new_cs(atom_i: usize, atom_j: usize) -> Self {
        TriboBond {
            atom_i,
            atom_j,
            length_eq: 1.82e-10,
            k_spring: 280.0,
            dissociation_energy: 2.9e-19,
            broken: false,
            bell: BellModelParams::carbon_sulfur(),
            cumulative_force: 0.0,
        }
    }

    /// Current stretching force (N) given actual bond vector.
    pub fn stretching_force(&self, pos_i: [f64; 3], pos_j: [f64; 3]) -> f64 {
        let r = norm3(sub3(pos_i, pos_j));
        self.k_spring * (r - self.length_eq)
    }

    /// Try to break this bond stochastically.
    ///
    /// `force` is the current tensile force (N), `t` is temperature (K),
    /// `dt` is the MD time step (s). Returns `true` if bond broke.
    pub fn attempt_break(&mut self, force: f64, t: f64, dt: f64, rng_val: f64) -> bool {
        if self.broken {
            return false;
        }
        let prob = self.bell.breaking_probability(force.max(0.0), t, dt);
        if rng_val < prob {
            self.broken = true;
            true
        } else {
            false
        }
    }
}

// ---------------------------------------------------------------------------
// MD atom
// ---------------------------------------------------------------------------

/// A single atom in a tribochemical MD simulation.
#[derive(Debug, Clone)]
pub struct TriboAtom {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m s⁻¹).
    pub vel: [f64; 3],
    /// Force accumulator (N).
    pub force: [f64; 3],
    /// Atom species.
    pub species: TriboChem,
    /// Charge (elementary units).
    pub charge: f64,
    /// Stress tensor contribution (Pa·m³), stored as 6 independent components
    /// \[xx, yy, zz, xy, xz, yz\].
    pub virial: [f64; 6],
}

impl TriboAtom {
    /// Create a new atom at rest.
    pub fn new(pos: [f64; 3], species: TriboChem) -> Self {
        TriboAtom {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            species,
            charge: 0.0,
            virial: [0.0; 6],
        }
    }

    /// Mass of this atom (kg).
    pub fn mass(&self) -> f64 {
        self.species.mass()
    }

    /// Kinetic energy (J).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass() * dot3(self.vel, self.vel)
    }

    /// Instantaneous temperature (K) from equipartition.
    pub fn instantaneous_temperature(&self) -> f64 {
        2.0 * self.kinetic_energy() / (3.0 * K_B)
    }

    /// Atomic stress (Pa) — hydrostatic component from virial.
    pub fn hydrostatic_stress(&self, volume: f64) -> f64 {
        -(self.virial[0] + self.virial[1] + self.virial[2]) / (3.0 * volume)
    }
}

// ---------------------------------------------------------------------------
// Lennard-Jones with reaction field (ReaxFF-like simplified)
// ---------------------------------------------------------------------------

/// Lennard-Jones parameters for a pair.
#[derive(Debug, Clone, Copy)]
pub struct LjParams {
    /// Well depth (J).
    pub epsilon: f64,
    /// Distance at zero potential (m).
    pub sigma: f64,
    /// Cutoff distance (m).
    pub cutoff: f64,
}

impl LjParams {
    /// C–C pair parameters.
    pub fn cc() -> Self {
        LjParams {
            epsilon: 3.6e-22,
            sigma: 3.4e-10,
            cutoff: 12.0e-10,
        }
    }

    /// C–Fe pair parameters.
    pub fn c_fe() -> Self {
        LjParams {
            epsilon: 5.2e-22,
            sigma: 2.9e-10,
            cutoff: 12.0e-10,
        }
    }

    /// Fe–Fe pair parameters.
    pub fn fe_fe() -> Self {
        LjParams {
            epsilon: 8.4e-22,
            sigma: 2.5e-10,
            cutoff: 12.0e-10,
        }
    }

    /// Lennard-Jones potential energy (J) at distance `r` (m).
    pub fn energy(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            return 0.0;
        }
        let sr = self.sigma / r;
        4.0 * self.epsilon * (sr.powi(12) - sr.powi(6))
    }

    /// Lennard-Jones force magnitude (N) — negative means attractive.
    ///
    /// Returns `dU/dr` (positive = repulsive).
    pub fn force_magnitude(&self, r: f64) -> f64 {
        if r >= self.cutoff {
            return 0.0;
        }
        let sr = self.sigma / r;
        4.0 * self.epsilon / r * (12.0 * sr.powi(12) - 6.0 * sr.powi(6))
    }

    /// Force vector on atom i from atom j.
    pub fn force_vector(&self, pos_i: [f64; 3], pos_j: [f64; 3]) -> [f64; 3] {
        let r_vec = sub3(pos_i, pos_j);
        let r = norm3(r_vec);
        if r < 1.0e-20 || r >= self.cutoff {
            return [0.0; 3];
        }
        let f_mag = self.force_magnitude(r);
        scale3(normalize3(r_vec), f_mag)
    }
}

// ---------------------------------------------------------------------------
// Morse potential (for reactive bonds)
// ---------------------------------------------------------------------------

/// Morse potential parameters for a reactive bond.
#[derive(Debug, Clone, Copy)]
pub struct MorseParams {
    /// Dissociation energy (J).
    pub de: f64,
    /// Width parameter (m⁻¹).
    pub alpha: f64,
    /// Equilibrium distance (m).
    pub r_eq: f64,
}

impl MorseParams {
    /// C–C Morse parameters.
    pub fn cc() -> Self {
        MorseParams {
            de: 6.03e-19,
            alpha: 1.8e10,
            r_eq: BOND_CC,
        }
    }

    /// C–O Morse parameters.
    pub fn co() -> Self {
        MorseParams {
            de: 5.73e-19,
            alpha: 2.3e10,
            r_eq: 1.43e-10,
        }
    }

    /// Morse potential energy (J).
    pub fn energy(&self, r: f64) -> f64 {
        let e = (-self.alpha * (r - self.r_eq)).exp();
        self.de * (1.0 - e) * (1.0 - e) - self.de
    }

    /// Morse force magnitude (N, positive = repulsive).
    pub fn force_magnitude(&self, r: f64) -> f64 {
        let e = (-self.alpha * (r - self.r_eq)).exp();
        2.0 * self.de * self.alpha * e * (1.0 - e)
    }

    /// Morse force vector on atom i from atom j.
    pub fn force_vector(&self, pos_i: [f64; 3], pos_j: [f64; 3]) -> [f64; 3] {
        let r_vec = sub3(pos_i, pos_j);
        let r = norm3(r_vec);
        if r < 1.0e-20 {
            return [0.0; 3];
        }
        scale3(normalize3(r_vec), self.force_magnitude(r))
    }
}

// ---------------------------------------------------------------------------
// ZDDP tribofilm model
// ---------------------------------------------------------------------------

/// Zinc dialkyldithiophosphate (ZDDP) tribofilm formation model.
///
/// Tracks the growth of a protective phosphate glass-like film on the substrate.
#[derive(Debug, Clone)]
pub struct ZddpTribofilm {
    /// Film thickness (m).
    pub thickness: f64,
    /// Film growth rate constant at reference conditions (m s⁻¹).
    pub k_growth: f64,
    /// Film removal rate constant (m s⁻¹ Pa⁻¹).
    pub k_removal: f64,
    /// Maximum film thickness (m).
    pub max_thickness: f64,
    /// Activation energy for film growth (J mol⁻¹).
    pub e_activation: f64,
    /// Film density (kg m⁻³).
    pub density: f64,
    /// Phosphate concentration in film (mass fraction).
    pub phosphate_fraction: f64,
    /// Zinc concentration in film (mass fraction).
    pub zinc_fraction: f64,
    /// Sulfur concentration in film (mass fraction).
    pub sulfur_fraction: f64,
}

impl ZddpTribofilm {
    /// Create a new ZDDP tribofilm tracker.
    pub fn new() -> Self {
        ZddpTribofilm {
            thickness: 0.0,
            k_growth: 1.0e-12,
            k_removal: 5.0e-19,
            max_thickness: 200.0e-9,
            e_activation: 60_000.0,
            density: 2800.0,
            phosphate_fraction: 0.35,
            zinc_fraction: 0.25,
            sulfur_fraction: 0.20,
        }
    }

    /// Film growth rate (m s⁻¹) at temperature `t` (K) and contact pressure `p` (Pa).
    pub fn growth_rate(&self, t: f64, p: f64) -> f64 {
        let arrhenius = (-self.e_activation / (R_GAS * t)).exp();
        let coverage = 1.0 - self.thickness / self.max_thickness;
        self.k_growth * arrhenius * p * coverage.max(0.0)
    }

    /// Film removal rate (m s⁻¹) under shear stress `tau` (Pa).
    pub fn removal_rate(&self, tau: f64) -> f64 {
        self.k_removal * tau * self.thickness
    }

    /// Advance film growth for time step `dt` (s).
    pub fn step(&mut self, t: f64, p: f64, tau: f64, dt: f64) {
        let grow = self.growth_rate(t, p);
        let remove = self.removal_rate(tau);
        let net = (grow - remove) * dt;
        self.thickness = (self.thickness + net).clamp(0.0, self.max_thickness);
    }

    /// Protective pressure load capacity (Pa) of the film.
    pub fn load_capacity(&self) -> f64 {
        let h_film = 5.0e9; // Film hardness ~5 GPa
        h_film * self.thickness / self.max_thickness
    }

    /// Film mass per unit area (kg m⁻²).
    pub fn areal_density(&self) -> f64 {
        self.density * self.thickness
    }
}

impl Default for ZddpTribofilm {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tribopolymer formation
// ---------------------------------------------------------------------------

/// Tribopolymer chain formed by shear-induced polymerization.
#[derive(Debug, Clone)]
pub struct TriboPolymer {
    /// Number of repeat units.
    pub degree_of_polymerization: usize,
    /// Monomer molecular weight (kg mol⁻¹).
    pub monomer_mw: f64,
    /// End-to-end distance (m).
    pub end_to_end: f64,
    /// Persistence length (m).
    pub persistence_length: f64,
    /// Is the chain crosslinked?
    pub crosslinked: bool,
    /// Shear viscosity contribution (Pa·s).
    pub viscosity_contribution: f64,
}

impl TriboPolymer {
    /// Create a new tribopolymer chain from a monomer.
    pub fn new(monomer_mw: f64) -> Self {
        TriboPolymer {
            degree_of_polymerization: 1,
            monomer_mw,
            end_to_end: 3.0e-10,
            persistence_length: 1.0e-9,
            crosslinked: false,
            viscosity_contribution: 0.0,
        }
    }

    /// Molecular weight of the chain (kg mol⁻¹).
    pub fn molecular_weight(&self) -> f64 {
        self.monomer_mw * self.degree_of_polymerization as f64
    }

    /// Add one repeat unit to the chain.
    pub fn propagate(&mut self, monomer_length: f64) {
        self.degree_of_polymerization += 1;
        let n = self.degree_of_polymerization as f64;
        // Freely-jointed chain: R ≈ l·√N
        self.end_to_end = monomer_length * n.sqrt();
    }

    /// Radius of gyration (m) for a Gaussian chain.
    pub fn radius_of_gyration(&self) -> f64 {
        let n = self.degree_of_polymerization as f64;
        let l = 1.54e-10; // backbone bond length
        l * n.sqrt() / (6.0f64).sqrt()
    }

    /// Reptation diffusion coefficient (m² s⁻¹) in a melt.
    pub fn reptation_diffusivity(&self, t: f64, friction: f64) -> f64 {
        let n = self.degree_of_polymerization as f64;
        K_B * t / (friction * n * n)
    }

    /// Viscosity contribution from this chain (Rouse model, Pa·s).
    pub fn rouse_viscosity(&self, t: f64, friction: f64, rho: f64) -> f64 {
        let n = self.degree_of_polymerization as f64;
        let mw = self.molecular_weight();
        rho * K_B * t * friction * n / (6.0 * mw * AMU)
    }
}

// ---------------------------------------------------------------------------
// Lubricant additive decomposition
// ---------------------------------------------------------------------------

/// Lubricant additive molecule undergoing stress-induced decomposition.
#[derive(Debug, Clone)]
pub struct LubricantAdditive {
    /// Additive type label.
    pub label: String,
    /// Molecular weight (kg mol⁻¹).
    pub mw: f64,
    /// Number of intact molecules.
    pub n_intact: usize,
    /// Number of decomposed molecules.
    pub n_decomposed: usize,
    /// Thermal decomposition rate (s⁻¹).
    pub k_thermal: f64,
    /// Mechanochemical rate sensitivity (m³).
    pub mechano_coupling: f64,
    /// Current shear stress (Pa).
    pub shear_stress: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl LubricantAdditive {
    /// Create a ZDDP-like additive.
    pub fn zddp(n: usize) -> Self {
        LubricantAdditive {
            label: "ZDDP".to_string(),
            mw: 0.302,
            n_intact: n,
            n_decomposed: 0,
            k_thermal: 1.0e-8,
            mechano_coupling: 5.0e-30,
            shear_stress: 0.0,
            temperature: 373.15,
        }
    }

    /// Create a friction modifier additive (glycerol monooleate proxy).
    pub fn friction_modifier(n: usize) -> Self {
        LubricantAdditive {
            label: "GMO".to_string(),
            mw: 0.356,
            n_intact: n,
            n_decomposed: 0,
            k_thermal: 1.0e-10,
            mechano_coupling: 2.0e-30,
            shear_stress: 0.0,
            temperature: 373.15,
        }
    }

    /// Total decomposition rate constant (s⁻¹) including mechanochemical term.
    pub fn total_rate(&self) -> f64 {
        let mechano =
            shear_rate_enhancement(self.mechano_coupling, self.shear_stress, self.temperature);
        self.k_thermal * mechano
    }

    /// Advance decomposition for time `dt` (s).
    pub fn step(&mut self, dt: f64) {
        let k = self.total_rate();
        let decompose = (self.n_intact as f64 * k * dt) as usize;
        let decompose = decompose.min(self.n_intact);
        self.n_intact -= decompose;
        self.n_decomposed += decompose;
    }

    /// Fractional conversion (0–1).
    pub fn conversion(&self) -> f64 {
        let total = self.n_intact + self.n_decomposed;
        if total == 0 {
            return 0.0;
        }
        self.n_decomposed as f64 / total as f64
    }
}

// ---------------------------------------------------------------------------
// Shear-induced phase transformation
// ---------------------------------------------------------------------------

/// Phase state of a tribological material.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriboPhase {
    /// Crystalline phase.
    Crystalline,
    /// Amorphous / glassy phase.
    Amorphous,
    /// Nanocrystalline.
    Nanocrystalline,
    /// Liquid-like shear band.
    LiquidLike,
}

/// Shear-induced phase transformation model.
///
/// Models the transformation from crystalline to amorphous under shear stress.
#[derive(Debug, Clone)]
pub struct ShearPhaseTransformer {
    /// Current phase.
    pub phase: TriboPhase,
    /// Volume fraction in amorphous state.
    pub amorphous_fraction: f64,
    /// Critical shear stress for transformation (Pa).
    pub tau_crit: f64,
    /// Transformation rate constant (Pa⁻¹ s⁻¹).
    pub k_transform: f64,
    /// Recovery rate constant (s⁻¹).
    pub k_recover: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl ShearPhaseTransformer {
    /// Create a new phase transformer for steel-like material.
    pub fn steel() -> Self {
        ShearPhaseTransformer {
            phase: TriboPhase::Crystalline,
            amorphous_fraction: 0.0,
            tau_crit: 5.0e9,
            k_transform: 1.0e-12,
            k_recover: 1.0e-5,
            temperature: 300.0,
        }
    }

    /// Create a model for an amorphous carbon coating (DLC).
    pub fn dlc() -> Self {
        ShearPhaseTransformer {
            phase: TriboPhase::Amorphous,
            amorphous_fraction: 1.0,
            tau_crit: 10.0e9,
            k_transform: 0.0,
            k_recover: 0.0,
            temperature: 300.0,
        }
    }

    /// Advance phase transformation for time step `dt` (s) under shear `tau` (Pa).
    pub fn step(&mut self, tau: f64, dt: f64) {
        let drive = (tau - self.tau_crit).max(0.0);
        let transform = self.k_transform * drive * (1.0 - self.amorphous_fraction) * dt;
        let recover = self.k_recover * self.amorphous_fraction * dt;
        self.amorphous_fraction = (self.amorphous_fraction + transform - recover).clamp(0.0, 1.0);
        self.phase = if self.amorphous_fraction > 0.9 {
            TriboPhase::Amorphous
        } else if self.amorphous_fraction > 0.5 {
            TriboPhase::Nanocrystalline
        } else if tau > 8.0e9 {
            TriboPhase::LiquidLike
        } else {
            TriboPhase::Crystalline
        };
    }

    /// Hardness modification factor from amorphous fraction.
    pub fn hardness_factor(&self) -> f64 {
        // Amorphous phases typically 20–30% softer
        1.0 - 0.25 * self.amorphous_fraction
    }
}

// ---------------------------------------------------------------------------
// Wear particle formation
// ---------------------------------------------------------------------------

/// A wear particle formed by tribochemical or mechanical wear.
#[derive(Debug, Clone)]
pub struct WearParticle {
    /// Position of detachment (m).
    pub origin: [f64; 3],
    /// Particle radius (m).
    pub radius: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Chemical composition tag.
    pub composition: WearComposition,
    /// Formation mechanism.
    pub mechanism: WearMechanism,
    /// Current velocity (m s⁻¹).
    pub velocity: [f64; 3],
}

/// Chemical composition of a wear particle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WearComposition {
    /// Metallic wear particle.
    Metal,
    /// Oxide wear particle (tribo-oxidation).
    Oxide,
    /// Carbide or graphitic particle.
    Carbon,
    /// Mixed organic-inorganic (from additive decomposition).
    Organometallic,
}

/// Mechanism by which a wear particle formed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WearMechanism {
    /// Adhesive wear (asperity junction).
    Adhesive,
    /// Abrasive wear (hard asperity plowing).
    Abrasive,
    /// Tribochemical wear (reaction-assisted).
    Tribochemical,
    /// Fatigue-induced delamination.
    Fatigue,
}

impl WearParticle {
    /// Create a metallic adhesive wear particle.
    pub fn adhesive(origin: [f64; 3], radius: f64) -> Self {
        let vol = (4.0 / 3.0) * PI * radius.powi(3);
        let rho_fe = 7874.0;
        WearParticle {
            origin,
            radius,
            mass: rho_fe * vol,
            composition: WearComposition::Metal,
            mechanism: WearMechanism::Adhesive,
            velocity: [0.0; 3],
        }
    }

    /// Create a tribochemically-formed oxide particle.
    pub fn oxide(origin: [f64; 3], radius: f64) -> Self {
        let vol = (4.0 / 3.0) * PI * radius.powi(3);
        let rho_fe2o3 = 5240.0;
        WearParticle {
            origin,
            radius,
            mass: rho_fe2o3 * vol,
            composition: WearComposition::Oxide,
            mechanism: WearMechanism::Tribochemical,
            velocity: [0.0; 3],
        }
    }

    /// Volume of the wear particle (m³).
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * PI * self.radius.powi(3)
    }
}

/// Archard wear model: wear volume rate (m³ s⁻¹).
///
/// `k_archard` wear coefficient (dimensionless), `load` (N),
/// `hardness` (Pa), `sliding_velocity` (m s⁻¹).
pub fn archard_wear_rate(k_archard: f64, load: f64, hardness: f64, sliding_velocity: f64) -> f64 {
    k_archard * load * sliding_velocity / hardness
}

/// Tribochemical wear rate enhanced by reaction.
///
/// Adds a chemical term `k_chem` (m s⁻¹ Pa⁻¹) times contact pressure.
pub fn tribochemical_wear_rate(
    k_archard: f64,
    load: f64,
    hardness: f64,
    sliding_velocity: f64,
    k_chem: f64,
    contact_pressure: f64,
) -> f64 {
    let mechanical = archard_wear_rate(k_archard, load, hardness, sliding_velocity);
    let chemical = k_chem * contact_pressure;
    mechanical + chemical
}

// ---------------------------------------------------------------------------
// Mechanically-activated reaction
// ---------------------------------------------------------------------------

/// Mechanically-activated reaction rate using Eyring model.
///
/// Extends the Bell model to include both activation volume and shear stress.
#[derive(Debug, Clone)]
pub struct EyringReaction {
    /// Pre-exponential factor (same units as rate, s⁻¹ typically).
    pub a0: f64,
    /// Thermal activation energy (J mol⁻¹).
    pub e_thermal: f64,
    /// Activation volume (m³ mol⁻¹) — mechanochemical coupling.
    pub v_activation: f64,
    /// Current temperature (K).
    pub temperature: f64,
    /// Current shear stress (Pa).
    pub shear_stress: f64,
    /// Current normal stress (Pa).
    pub normal_stress: f64,
}

impl EyringReaction {
    /// Create a new Eyring reaction model.
    pub fn new(a0: f64, e_thermal: f64, v_activation: f64, temperature: f64) -> Self {
        EyringReaction {
            a0,
            e_thermal,
            v_activation,
            temperature,
            shear_stress: 0.0,
            normal_stress: 0.0,
        }
    }

    /// Forward reaction rate (s⁻¹).
    pub fn forward_rate(&self) -> f64 {
        let e_eff = self.e_thermal - self.v_activation * self.shear_stress;
        self.a0 * (-e_eff / (R_GAS * self.temperature)).exp()
    }

    /// Reverse reaction rate (s⁻¹).
    pub fn reverse_rate(&self) -> f64 {
        let e_eff = self.e_thermal + self.v_activation * self.shear_stress;
        self.a0 * (-e_eff / (R_GAS * self.temperature)).exp()
    }

    /// Net reaction rate (s⁻¹), positive = forward.
    pub fn net_rate(&self) -> f64 {
        self.forward_rate() - self.reverse_rate()
    }

    /// Stress-induced rate enhancement over thermal rate.
    pub fn rate_enhancement(&self) -> f64 {
        let k_thermal = self.a0 * (-self.e_thermal / (R_GAS * self.temperature)).exp();
        if k_thermal < 1.0e-200 {
            return 1.0;
        }
        self.forward_rate() / k_thermal
    }
}

// ---------------------------------------------------------------------------
// Simulation box and shear
// ---------------------------------------------------------------------------

/// Periodic simulation box with optional Lees-Edwards shear boundary.
#[derive(Debug, Clone)]
pub struct TriboBox {
    /// Box lengths \[Lx, Ly, Lz\] (m).
    pub lengths: [f64; 3],
    /// Applied shear strain rate (s⁻¹) in x-direction of xz-plane.
    pub shear_rate: f64,
    /// Accumulated shear strain.
    pub strain: f64,
    /// Current box tilt factor (Lees-Edwards offset in x per Lz).
    pub tilt: f64,
}

impl TriboBox {
    /// Create a cubic simulation box.
    pub fn cubic(length: f64) -> Self {
        TriboBox {
            lengths: [length; 3],
            shear_rate: 0.0,
            strain: 0.0,
            tilt: 0.0,
        }
    }

    /// Apply Lees-Edwards boundary conditions: returns wrapped position.
    pub fn wrap_lees_edwards(&self, pos: [f64; 3], _vel: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let mut p = pos;
        let mut v = _vel;
        // Wrap in z first
        let lz = self.lengths[2];
        if p[2] < 0.0 {
            p[2] += lz;
            p[0] += self.tilt * lz;
            v[0] += self.shear_rate * lz;
        } else if p[2] >= lz {
            p[2] -= lz;
            p[0] -= self.tilt * lz;
            v[0] -= self.shear_rate * lz;
        }
        // Wrap x, y
        for k in 0..2 {
            while p[k] < 0.0 {
                p[k] += self.lengths[k];
            }
            while p[k] >= self.lengths[k] {
                p[k] -= self.lengths[k];
            }
        }
        (p, v)
    }

    /// Advance the shear strain by time step `dt` (s).
    pub fn shear_step(&mut self, dt: f64) {
        self.strain += self.shear_rate * dt;
        self.tilt += self.shear_rate * dt;
        // Fold tilt back to (-0.5, 0.5) range
        while self.tilt > 0.5 {
            self.tilt -= 1.0;
        }
        while self.tilt < -0.5 {
            self.tilt += 1.0;
        }
    }

    /// Volume of the simulation box (m³).
    pub fn volume(&self) -> f64 {
        self.lengths[0] * self.lengths[1] * self.lengths[2]
    }
}

// ---------------------------------------------------------------------------
// Tribochemical simulation
// ---------------------------------------------------------------------------

/// Configuration for a tribochemical MD simulation.
#[derive(Debug, Clone)]
pub struct TriboConfig {
    /// Time step (s).
    pub dt: f64,
    /// System temperature (K).
    pub temperature: f64,
    /// Applied normal stress (Pa).
    pub normal_stress: f64,
    /// Applied sliding velocity (m s⁻¹).
    pub sliding_velocity: f64,
    /// Thermostat relaxation time (s).
    pub tau_thermostat: f64,
    /// Enable bond breaking.
    pub reactive: bool,
    /// Enable ZDDP film growth.
    pub zddp_film: bool,
}

impl Default for TriboConfig {
    fn default() -> Self {
        TriboConfig {
            dt: 1.0e-15,
            temperature: 373.15,
            normal_stress: 1.0e9,
            sliding_velocity: 1.0,
            tau_thermostat: 0.1e-12,
            reactive: true,
            zddp_film: true,
        }
    }
}

/// Tribochemical MD simulation state.
#[derive(Debug, Clone)]
pub struct TriboSimulation {
    /// Simulation configuration.
    pub config: TriboConfig,
    /// All atoms.
    pub atoms: Vec<TriboAtom>,
    /// All reactive bonds.
    pub bonds: Vec<TriboBond>,
    /// ZDDP tribofilm.
    pub film: ZddpTribofilm,
    /// Lubricant additives.
    pub additives: Vec<LubricantAdditive>,
    /// Formed wear particles.
    pub wear_particles: Vec<WearParticle>,
    /// Phase transformer.
    pub phase: ShearPhaseTransformer,
    /// Simulation box.
    pub box_: TriboBox,
    /// Current simulation time (s).
    pub time: f64,
    /// Step counter.
    pub step: usize,
    /// Shear stress history (Pa).
    pub shear_stress_history: Vec<f64>,
}

impl TriboSimulation {
    /// Create a new tribochemical simulation.
    pub fn new(config: TriboConfig, box_length: f64) -> Self {
        TriboSimulation {
            config,
            atoms: Vec::new(),
            bonds: Vec::new(),
            film: ZddpTribofilm::new(),
            additives: Vec::new(),
            wear_particles: Vec::new(),
            phase: ShearPhaseTransformer::steel(),
            box_: TriboBox::cubic(box_length),
            time: 0.0,
            step: 0,
            shear_stress_history: Vec::new(),
        }
    }

    /// Add an atom to the simulation.
    pub fn add_atom(&mut self, atom: TriboAtom) {
        self.atoms.push(atom);
    }

    /// Add a reactive bond.
    pub fn add_bond(&mut self, bond: TriboBond) {
        self.bonds.push(bond);
    }

    /// Add a lubricant additive batch.
    pub fn add_additive(&mut self, add: LubricantAdditive) {
        self.additives.push(add);
    }

    /// Total kinetic energy (J).
    pub fn kinetic_energy(&self) -> f64 {
        self.atoms.iter().map(|a| a.kinetic_energy()).sum()
    }

    /// Instantaneous temperature (K).
    pub fn temperature(&self) -> f64 {
        let n = self.atoms.len();
        if n == 0 {
            return 0.0;
        }
        let ke: f64 = self.atoms.iter().map(|a| a.kinetic_energy()).sum();
        2.0 * ke / (3.0 * n as f64 * K_B)
    }

    /// Virial pressure (Pa).
    pub fn pressure(&self) -> f64 {
        let vol = self.box_.volume();
        let n = self.atoms.len();
        if vol < 1.0e-60 {
            return 0.0;
        }
        let ke: f64 = self.atoms.iter().map(|a| a.kinetic_energy()).sum();
        let virial: f64 = self
            .atoms
            .iter()
            .map(|a| a.virial[0] + a.virial[1] + a.virial[2])
            .sum();
        (2.0 * ke + virial) / (3.0 * vol * n.max(1) as f64)
    }

    /// Current shear stress estimate (Pa).
    pub fn shear_stress(&self) -> f64 {
        let vol = self.box_.volume();
        if vol < 1.0e-60 {
            return 0.0;
        }
        let virial_xz: f64 = self.atoms.iter().map(|a| a.virial[4]).sum();
        virial_xz.abs() / vol
    }

    /// Apply Berendsen velocity scaling thermostat.
    fn apply_thermostat(&mut self) {
        let t_inst = self.temperature();
        let t_target = self.config.temperature;
        if t_inst < 1.0e-10 {
            return;
        }
        let lambda = (1.0
            + (self.config.dt / self.config.tau_thermostat) * (t_target / t_inst - 1.0))
            .max(0.0)
            .sqrt();
        for atom in &mut self.atoms {
            atom.vel = scale3(atom.vel, lambda);
        }
    }

    /// Compute forces from all bond interactions.
    fn compute_bond_forces(&mut self) {
        // Reset forces
        for atom in &mut self.atoms {
            atom.force = [0.0; 3];
        }
        let positions: Vec<[f64; 3]> = self.atoms.iter().map(|a| a.pos).collect();
        for bond in &self.bonds {
            if bond.broken {
                continue;
            }
            let i = bond.atom_i;
            let j = bond.atom_j;
            if i >= positions.len() || j >= positions.len() {
                continue;
            }
            let r_vec = sub3(positions[i], positions[j]);
            let r = norm3(r_vec);
            if r < 1.0e-20 {
                continue;
            }
            let f_mag = bond.k_spring * (r - bond.length_eq);
            let f_vec = scale3(normalize3(r_vec), -f_mag);
            self.atoms[i].force = add3(self.atoms[i].force, f_vec);
            self.atoms[j].force = add3(self.atoms[j].force, scale3(f_vec, -1.0));
        }
    }

    /// Advance simulation by one MD step.
    pub fn step_forward(&mut self, rng_vals: &[f64]) {
        let dt = self.config.dt;
        let t = self.config.temperature;

        // 1. Compute forces
        self.compute_bond_forces();

        // 2. Velocity Verlet — half-step velocity
        for atom in &mut self.atoms {
            let m = atom.mass();
            atom.vel = add3(atom.vel, scale3(atom.force, 0.5 * dt / m));
        }

        // 3. Update positions
        let positions: Vec<[f64; 3]> = self.atoms.iter().map(|a| a.pos).collect();
        let velocities: Vec<[f64; 3]> = self.atoms.iter().map(|a| a.vel).collect();
        for (i, atom) in self.atoms.iter_mut().enumerate() {
            let new_pos = add3(positions[i], scale3(velocities[i], dt));
            let (wp, wv) = self.box_.wrap_lees_edwards(new_pos, velocities[i]);
            atom.pos = wp;
            atom.vel = wv;
        }

        // 4. Recompute forces
        self.compute_bond_forces();

        // 5. Complete velocity Verlet
        for atom in &mut self.atoms {
            let m = atom.mass();
            atom.vel = add3(atom.vel, scale3(atom.force, 0.5 * dt / m));
        }

        // 6. Thermostat
        self.apply_thermostat();

        // 7. Bond breaking (stochastic)
        if self.config.reactive {
            let positions2: Vec<[f64; 3]> = self.atoms.iter().map(|a| a.pos).collect();
            for (idx, bond) in self.bonds.iter_mut().enumerate() {
                if bond.broken {
                    continue;
                }
                let i = bond.atom_i;
                let j = bond.atom_j;
                if i >= positions2.len() || j >= positions2.len() {
                    continue;
                }
                let f = bond.stretching_force(positions2[i], positions2[j]);
                let rng = rng_vals.get(idx).copied().unwrap_or(0.5);
                bond.attempt_break(f, t, dt, rng);
            }
        }

        // 8. ZDDP film growth
        if self.config.zddp_film {
            let tau = self.shear_stress();
            self.film.step(t, self.config.normal_stress, tau, dt);
        }

        // 9. Additive decomposition
        for add in &mut self.additives {
            add.shear_stress = self.config.normal_stress * 0.1;
            add.temperature = t;
            add.step(dt);
        }

        // 10. Phase transformation
        let tau = self.config.normal_stress * 0.15;
        self.phase.step(tau, dt);

        // 11. Shear box advance
        let sr = self.config.sliding_velocity / self.box_.lengths[2];
        self.box_.shear_rate = sr;
        self.box_.shear_step(dt);

        self.shear_stress_history.push(self.shear_stress());
        self.time += dt;
        self.step += 1;
    }

    /// Number of broken bonds.
    pub fn n_broken_bonds(&self) -> usize {
        self.bonds.iter().filter(|b| b.broken).count()
    }

    /// Mean sliding distance (m) accumulated.
    pub fn sliding_distance(&self) -> f64 {
        self.config.sliding_velocity * self.time
    }
}

// ---------------------------------------------------------------------------
// Tribocorrosion model
// ---------------------------------------------------------------------------

/// Tribocorrosion model coupling mechanical wear and electrochemical dissolution.
#[derive(Debug, Clone)]
pub struct TribocorrosionModel {
    /// Passive film thickness (m).
    pub passive_film_thickness: f64,
    /// Film repassivation rate constant (m s⁻¹).
    pub k_repassivate: f64,
    /// Film removal rate per unit wear (dimensionless).
    pub wear_synergy: f64,
    /// Corrosion current density (A m⁻²).
    pub corrosion_current: f64,
    /// Faraday constant (C mol⁻¹).
    pub faraday: f64,
    /// Molar mass of Fe (kg mol⁻¹).
    pub molar_mass_fe: f64,
    /// Valence of Fe dissolution.
    pub valence: f64,
    /// Metal density (kg m⁻³).
    pub rho_metal: f64,
}

impl TribocorrosionModel {
    /// Create a default stainless steel tribocorrosion model.
    pub fn stainless_steel() -> Self {
        TribocorrosionModel {
            passive_film_thickness: 5.0e-9,
            k_repassivate: 1.0e-9,
            wear_synergy: 0.1,
            corrosion_current: 1.0e-4,
            faraday: 96_485.0,
            molar_mass_fe: 0.05585,
            valence: 2.0,
            rho_metal: 7900.0,
        }
    }

    /// Corrosion-enhanced wear rate (m³ s⁻¹ m⁻²).
    pub fn corrosion_wear_rate(&self) -> f64 {
        self.molar_mass_fe * self.corrosion_current / (self.valence * self.faraday * self.rho_metal)
    }

    /// Mechanical wear rate (m³ s⁻¹ m⁻²) from Archard model.
    pub fn mechanical_wear_rate(&self, k: f64, p: f64, v: f64, h: f64) -> f64 {
        k * p * v / h
    }

    /// Total synergistic wear rate (m³ s⁻¹ m⁻²).
    pub fn total_wear_rate(&self, k: f64, p: f64, v: f64, h: f64) -> f64 {
        let mech = self.mechanical_wear_rate(k, p, v, h);
        let chem = self.corrosion_wear_rate();
        let synergy = self.wear_synergy * (mech + chem);
        mech + chem + synergy
    }

    /// Repassivation time (s) after film removal.
    pub fn repassivation_time(&self) -> f64 {
        self.passive_film_thickness / self.k_repassivate
    }
}

// ---------------------------------------------------------------------------
// Stribeck curve
// ---------------------------------------------------------------------------

/// Stribeck curve parameters for a lubricated contact.
#[derive(Debug, Clone)]
pub struct StribeckCurve {
    /// Dynamic viscosity of lubricant (Pa·s).
    pub viscosity: f64,
    /// Sliding speed (m s⁻¹).
    pub speed: f64,
    /// Normal load (N).
    pub load: f64,
    /// Contact radius (m).
    pub contact_radius: f64,
}

impl StribeckCurve {
    /// Create a typical oil-lubricated steel-on-steel contact.
    pub fn steel_oil(load: f64, speed: f64) -> Self {
        StribeckCurve {
            viscosity: 0.04,
            speed,
            load,
            contact_radius: 1.0e-3,
        }
    }

    /// Stribeck number (dimensionless).
    pub fn stribeck_number(&self) -> f64 {
        self.viscosity * self.speed / (self.load / (PI * self.contact_radius * self.contact_radius))
    }

    /// Friction coefficient from Stribeck curve approximation.
    pub fn friction_coefficient(&self) -> f64 {
        let sn = self.stribeck_number();
        // Boundary: high friction, EHD: low, mixed: intermediate
        let mu_boundary = 0.12;
        let mu_ehd = 0.004;
        let sn_crit = 1.0e-8;
        if sn < sn_crit {
            mu_boundary
        } else {
            mu_ehd + (mu_boundary - mu_ehd) * (-(sn / sn_crit).ln().max(0.0) / 10.0).exp()
        }
    }

    /// Film thickness (m) from Hamrock-Dowson correlation (point contact).
    pub fn film_thickness(&self) -> f64 {
        let u_star = self.viscosity * self.speed; // entrainment parameter
        let r = self.contact_radius;
        let e_prime = 200.0e9; // reduced modulus
        let w = self.load;
        let u_bar = u_star / (e_prime * r);
        let w_bar = w / (e_prime * r * r);
        // Central film thickness: h_c = 2.69 R U^0.67 G^0.53 W^(-0.067)
        2.69 * r * u_bar.powf(0.67) * (1.0e4f64).powf(0.53) * w_bar.powf(-0.067)
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bell_model_zero_force() {
        let bell = BellModelParams::carbon_carbon();
        let k0 = bell.breaking_rate(0.0, 300.0);
        // Rate should be very small for C-C bond at room temperature
        assert!(k0 < 1.0e-5, "Zero-force rate should be tiny: {k0}");
    }

    #[test]
    fn test_bell_model_force_increases_rate() {
        let bell = BellModelParams::carbon_carbon();
        let k_low = bell.breaking_rate(0.0, 300.0);
        let k_high = bell.breaking_rate(1.0e-9, 300.0);
        assert!(k_high > k_low, "Higher force should increase rate");
    }

    #[test]
    fn test_bell_model_lifetime_positive() {
        let bell = BellModelParams::carbon_sulfur();
        let tau = bell.mean_lifetime(0.0, 300.0);
        assert!(tau > 0.0);
    }

    #[test]
    fn test_bell_probability_in_range() {
        let bell = BellModelParams::phosphorus_sulfur();
        let p = bell.breaking_probability(5.0e-10, 400.0, 1.0e-12);
        assert!((0.0..=1.0).contains(&p), "Probability out of range: {p}");
    }

    #[test]
    fn test_hertz_contact_pressure_positive() {
        let p0 = hertz_contact_pressure(100.0e9, 1.0e-3, 1.0);
        assert!(p0 > 0.0);
    }

    #[test]
    fn test_hertz_contact_radius_positive() {
        let a = hertz_contact_radius(100.0e9, 1.0e-3, 10.0);
        assert!(a > 0.0 && a < 1.0e-3, "Contact radius: {a}");
    }

    #[test]
    fn test_shear_rate_enhancement_unity_at_zero() {
        let f = shear_rate_enhancement(1.0e-29, 0.0, 300.0);
        assert!((f - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_chemical_potential_compression() {
        let mu = chemical_potential_stressed(0.0, 1.0e9, 1.0e-5);
        assert!(mu < 0.0, "Compression should lower chemical potential");
    }

    #[test]
    fn test_lj_params_energy_zero_outside_cutoff() {
        let lj = LjParams::cc();
        assert_eq!(lj.energy(lj.cutoff + 1.0e-10), 0.0);
    }

    #[test]
    fn test_lj_params_energy_minimum() {
        let lj = LjParams::cc();
        let r_min = 2.0f64.powf(1.0 / 6.0) * lj.sigma;
        let e_min = lj.energy(r_min);
        assert!(
            (e_min + lj.epsilon).abs() < 1.0e-24,
            "LJ minimum energy: {e_min}"
        );
    }

    #[test]
    fn test_morse_energy_at_equilibrium() {
        let morse = MorseParams::cc();
        let e_eq = morse.energy(morse.r_eq);
        assert!(
            (e_eq + morse.de).abs() < 1.0e-30,
            "Morse energy at r_eq: {e_eq}"
        );
    }

    #[test]
    fn test_morse_force_zero_at_equilibrium() {
        let morse = MorseParams::cc();
        let f = morse.force_magnitude(morse.r_eq);
        assert!(f.abs() < 1.0e-5, "Morse force at r_eq should be zero: {f}");
    }

    #[test]
    fn test_zddp_film_growth() {
        let mut film = ZddpTribofilm::new();
        let t0 = film.thickness;
        film.step(373.15, 1.0e9, 1.0e8, 1.0);
        assert!(film.thickness > t0, "Film should grow");
    }

    #[test]
    fn test_zddp_film_max_thickness() {
        let mut film = ZddpTribofilm::new();
        for _ in 0..10000 {
            film.step(373.15, 5.0e9, 0.0, 100.0);
        }
        assert!(film.thickness <= film.max_thickness);
    }

    #[test]
    fn test_zddp_areal_density_positive() {
        let mut film = ZddpTribofilm::new();
        film.step(373.15, 1.0e9, 1.0e8, 1.0);
        assert!(film.areal_density() > 0.0);
    }

    #[test]
    fn test_tribopolymer_propagate() {
        let mut poly = TriboPolymer::new(0.1);
        let n0 = poly.degree_of_polymerization;
        poly.propagate(1.54e-10);
        assert_eq!(poly.degree_of_polymerization, n0 + 1);
    }

    #[test]
    fn test_tribopolymer_mw() {
        let mut poly = TriboPolymer::new(0.1);
        poly.propagate(1.54e-10);
        assert!((poly.molecular_weight() - 0.2).abs() < 1.0e-10);
    }

    #[test]
    fn test_lubricant_additive_decomposition() {
        // Use high shear stress (21 GPa) and long time steps (1 ms) so the
        // mechanochemical rate is large enough to produce measurable decomposition.
        let mut add = LubricantAdditive::zddp(1000);
        for _ in 0..100 {
            add.shear_stress = 2.1e10; // 21 GPa
            add.step(1.0e-3); // 1 ms per step
        }
        assert!(add.n_decomposed > 0, "Some molecules should decompose");
    }

    #[test]
    fn test_lubricant_additive_conservation() {
        let n0 = 500;
        let mut add = LubricantAdditive::friction_modifier(n0);
        add.step(1.0e-9);
        assert_eq!(add.n_intact + add.n_decomposed, n0);
    }

    #[test]
    fn test_shear_phase_transformation() {
        let mut spt = ShearPhaseTransformer::steel();
        assert_eq!(spt.phase, TriboPhase::Crystalline);
        // Apply huge stress many times
        for _ in 0..100_000 {
            spt.step(20.0e9, 1.0e-12);
        }
        assert!(
            spt.amorphous_fraction > 0.0,
            "Should transform under extreme stress"
        );
    }

    #[test]
    fn test_wear_particle_volume() {
        let wp = WearParticle::adhesive([0.0; 3], 1.0e-7);
        let v_expected = (4.0 / 3.0) * PI * (1.0e-7f64).powi(3);
        assert!((wp.volume() - v_expected).abs() < 1.0e-30);
    }

    #[test]
    fn test_archard_wear_rate_positive() {
        let wr = archard_wear_rate(1.0e-4, 10.0, 1.0e9, 1.0);
        assert!(wr > 0.0);
    }

    #[test]
    fn test_eyring_reaction_forward_gt_reverse() {
        let mut rxn = EyringReaction::new(1.0e12, 60_000.0, 1.0e-5, 400.0);
        rxn.shear_stress = 1.0e9;
        assert!(
            rxn.forward_rate() > rxn.reverse_rate(),
            "Forward should dominate under shear"
        );
    }

    #[test]
    fn test_eyring_rate_enhancement_gt_one() {
        let mut rxn = EyringReaction::new(1.0e12, 50_000.0, 1.0e-5, 400.0);
        rxn.shear_stress = 1.0e9;
        let enh = rxn.rate_enhancement();
        assert!(enh >= 1.0, "Rate should be enhanced: {enh}");
    }

    #[test]
    fn test_tribo_box_volume() {
        let b = TriboBox::cubic(1.0e-9);
        assert!((b.volume() - 1.0e-27).abs() < 1.0e-40);
    }

    #[test]
    fn test_tribocorrosion_repassivation_time() {
        let tc = TribocorrosionModel::stainless_steel();
        let t = tc.repassivation_time();
        assert!(t > 0.0 && t < 100.0, "Repassivation time: {t} s");
    }

    #[test]
    fn test_tribocorrosion_total_wear() {
        let tc = TribocorrosionModel::stainless_steel();
        let total = tc.total_wear_rate(1.0e-4, 1.0e9, 1.0, 1.0e10);
        let mech = tc.mechanical_wear_rate(1.0e-4, 1.0e9, 1.0, 1.0e10);
        assert!(
            total > mech,
            "Total wear should exceed mechanical alone (synergy)"
        );
    }

    #[test]
    fn test_stribeck_number_positive() {
        let sc = StribeckCurve::steel_oil(100.0, 0.1);
        let sn = sc.stribeck_number();
        assert!(sn > 0.0);
    }

    #[test]
    fn test_stribeck_friction_in_range() {
        let sc = StribeckCurve::steel_oil(100.0, 0.01);
        let mu = sc.friction_coefficient();
        assert!(mu > 0.001 && mu < 0.2, "Friction coefficient: {mu}");
    }

    #[test]
    fn test_tribo_simulation_step() {
        let config = TriboConfig {
            dt: 1.0e-15,
            temperature: 300.0,
            ..Default::default()
        };
        let mut sim = TriboSimulation::new(config, 5.0e-9);
        let mut a1 = TriboAtom::new([0.0, 0.0, 0.0], TriboChem::C);
        let mut a2 = TriboAtom::new([1.54e-10, 0.0, 0.0], TriboChem::C);
        a1.vel = [0.0, 0.0, 0.0];
        a2.vel = [100.0, 0.0, 0.0];
        sim.add_atom(a1);
        sim.add_atom(a2);
        sim.add_bond(TriboBond::new_cc(0, 1));
        sim.step_forward(&[0.9999]);
        assert_eq!(sim.step, 1);
        assert!(sim.time > 0.0);
    }

    #[test]
    fn test_tribo_atom_kinetic_energy() {
        let mut atom = TriboAtom::new([0.0; 3], TriboChem::Fe);
        atom.vel = [100.0, 0.0, 0.0];
        let ke = atom.kinetic_energy();
        let expected = 0.5 * TriboChem::Fe.mass() * 10_000.0;
        assert!((ke - expected).abs() < 1.0e-35, "KE: {ke} vs {expected}");
    }
}
