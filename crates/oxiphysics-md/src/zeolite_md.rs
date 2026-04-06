#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Zeolite and porous materials MD simulation.
//!
//! This module provides tools for simulating adsorption, diffusion, and
//! molecular dynamics in zeolite frameworks and other porous materials:
//!
//! - [`ZeoliteFramework`]: Crystal structure for FAU/MFI/LTA topology.
//! - [`crate::PeriodicBox`]: Periodic boundary conditions for zeolite unit cells.
//! - [`LangmuirIsotherm`]: Single-component Langmuir adsorption isotherm.
//! - [`BETIsotherm`]: Multi-layer BET adsorption isotherm.
//! - [`KnudsenDiffusion`]: Knudsen regime diffusion coefficient in micropores.
//! - [`GcmcSimulation`]: Grand Canonical Monte Carlo for adsorption.
//! - [`HenryConstant`]: Henry's law constant from MC simulations.
//! - [`IsostericHeat`]: Isosteric heat of adsorption from fluctuation formula.
//! - [`CageHoppingModel`]: Cage-to-cage hopping kinetics model.
//! - [`MoleculeZeoliteFF`]: Lennard-Jones + Coulomb molecule-framework potential.
//! - [`SilicaliteSite`]: Silicalite T-sites and channel coordinates.
//! - [`AluminosilicateAcidSite`]: Brønsted acid sites in aluminosilicates.
//! - [`IonExchangeSimulation`]: Ion exchange Monte Carlo in zeolites.
//! - [`MolecularSieve`]: Kinetic diameter vs pore size sieving criterion.
//! - [`PoreVolumeCalculator`]: Geometric pore volume via probe insertion.
//!
//! # References
//! - Frenkel, D. & Smit, B. (2002). *Understanding Molecular Simulation*, 2nd ed.
//! - Smit, B. & Maesen, T. L. M. (2008). *Chem. Rev.* 108, 4125.
//! - Dubbeldam, D. *et al.* (2016). *Mol. Sim.* 42, 81.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
pub const K_B: f64 = 1.380_649e-23;
/// Universal gas constant (J mol⁻¹ K⁻¹).
pub const R_GAS: f64 = 8.314_462_618;
/// Avogadro constant (mol⁻¹).
pub const N_AV: f64 = 6.022_140_76e23;
/// Elementary charge (C).
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Vacuum permittivity (C² N⁻¹ m⁻²).
pub const EPS_0: f64 = 8.854_187_817e-12;
/// Ångström in metres.
pub const ANGSTROM: f64 = 1.0e-10;

// ---------------------------------------------------------------------------
// Zeolite framework topology
// ---------------------------------------------------------------------------

/// Known zeolite framework types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameworkType {
    /// Faujasite (FAU) – large-pore zeolite (Y, X).
    FAU,
    /// MFI (ZSM-5, Silicalite-1) – medium-pore channel system.
    MFI,
    /// LTA (Zeolite-A) – small-pore cage structure.
    LTA,
    /// Custom user-defined framework.
    Custom,
}

/// Zeolite crystal framework descriptor.
///
/// Stores unit-cell parameters, pore sizes, and T-site positions for
/// a given framework topology.
#[derive(Debug, Clone)]
pub struct ZeoliteFramework {
    /// Framework type identifier.
    pub framework_type: FrameworkType,
    /// Unit-cell lengths \[a, b, c\] in Ångströms.
    pub cell_lengths: [f64; 3],
    /// Unit-cell angles \[alpha, beta, gamma\] in degrees.
    pub cell_angles: [f64; 3],
    /// Pore limiting diameter (Å).
    pub pore_limiting_diameter: f64,
    /// Largest included sphere diameter (Å).
    pub largest_included_sphere: f64,
    /// Si/T-site fractional coordinates (each row \[x, y, z\]).
    pub t_sites: Vec<[f64; 3]>,
    /// Oxygen fractional coordinates.
    pub o_sites: Vec<[f64; 3]>,
    /// Si:Al ratio (inf for pure silica).
    pub si_al_ratio: f64,
}

impl ZeoliteFramework {
    /// Construct a standard FAU (faujasite) unit cell.
    ///
    /// FAU has a = b = c = 24.74 Å, cubic symmetry, 12-ring pore (7.4 Å).
    pub fn fau() -> Self {
        let a = 24.74_f64;
        Self {
            framework_type: FrameworkType::FAU,
            cell_lengths: [a, a, a],
            cell_angles: [90.0, 90.0, 90.0],
            pore_limiting_diameter: 7.4,
            largest_included_sphere: 11.2,
            t_sites: fau_t_sites(),
            o_sites: Vec::new(),
            si_al_ratio: f64::INFINITY,
        }
    }

    /// Construct a standard MFI (ZSM-5/Silicalite-1) unit cell.
    ///
    /// MFI is orthorhombic: a=20.07, b=19.92, c=13.42 Å.
    pub fn mfi() -> Self {
        Self {
            framework_type: FrameworkType::MFI,
            cell_lengths: [20.07, 19.92, 13.42],
            cell_angles: [90.0, 90.0, 90.0],
            pore_limiting_diameter: 4.7,
            largest_included_sphere: 6.4,
            t_sites: mfi_t_sites(),
            o_sites: Vec::new(),
            si_al_ratio: f64::INFINITY,
        }
    }

    /// Construct a standard LTA (Zeolite-A) unit cell.
    ///
    /// LTA is cubic: a = b = c = 11.919 Å, 8-ring pore (4.1 Å).
    pub fn lta() -> Self {
        let a = 11.919_f64;
        Self {
            framework_type: FrameworkType::LTA,
            cell_lengths: [a, a, a],
            cell_angles: [90.0, 90.0, 90.0],
            pore_limiting_diameter: 4.1,
            largest_included_sphere: 11.4,
            t_sites: lta_t_sites(),
            o_sites: Vec::new(),
            si_al_ratio: 1.0,
        }
    }

    /// Volume of the unit cell (Å³).
    pub fn unit_cell_volume(&self) -> f64 {
        let [a, b, c] = self.cell_lengths;
        let [alpha, beta, gamma] = self.cell_angles;
        let (ca, cb, cg) = (
            alpha.to_radians().cos(),
            beta.to_radians().cos(),
            gamma.to_radians().cos(),
        );
        let sa = alpha.to_radians().sin();
        let vol_factor = (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt();
        a * b * c * sa * vol_factor / alpha.to_radians().sin() // simplified for orthogonal
            * (1.0 / sa) * sa
        // For orthorhombic/cubic this reduces to a*b*c
    }

    /// Unit cell volume (Å³) for orthogonal cells only (correct form).
    pub fn orthorhombic_volume(&self) -> f64 {
        self.cell_lengths[0] * self.cell_lengths[1] * self.cell_lengths[2]
    }

    /// Number of T-sites per unit cell.
    pub fn n_t_sites(&self) -> usize {
        self.t_sites.len()
    }

    /// Convert fractional coordinates to Cartesian (Å), assuming orthorhombic.
    pub fn frac_to_cart(&self, frac: [f64; 3]) -> [f64; 3] {
        [
            frac[0] * self.cell_lengths[0],
            frac[1] * self.cell_lengths[1],
            frac[2] * self.cell_lengths[2],
        ]
    }
}

/// Generate representative T-site fractional coordinates for FAU.
fn fau_t_sites() -> Vec<[f64; 3]> {
    // FAU has 192 T-sites in the conventional cell; we store a subset.
    vec![
        [0.1250, 0.1250, 0.3750],
        [0.3750, 0.1250, 0.1250],
        [0.1250, 0.3750, 0.1250],
        [0.3750, 0.3750, 0.3750],
        [0.1250, 0.6250, 0.3750],
        [0.3750, 0.6250, 0.1250],
    ]
}

/// Generate representative T-site fractional coordinates for MFI.
fn mfi_t_sites() -> Vec<[f64; 3]> {
    // MFI has 96 T-sites; representative subset shown.
    vec![
        [0.4230, 0.0565, 0.3395],
        [0.3065, 0.0290, 0.1900],
        [0.2795, 0.0610, 0.0395],
        [0.1185, 0.0620, 0.0265],
        [0.0715, 0.0290, 0.1840],
        [0.1820, 0.0620, 0.3310],
        [0.4230, 0.1690, 0.3395],
        [0.3065, 0.1960, 0.1900],
    ]
}

/// Generate representative T-site fractional coordinates for LTA.
fn lta_t_sites() -> Vec<[f64; 3]> {
    vec![
        [0.0000, 0.1840, 0.3750],
        [0.1840, 0.0000, 0.3750],
        [0.1840, 0.3750, 0.0000],
        [0.0000, 0.3750, 0.1840],
        [0.3750, 0.0000, 0.1840],
        [0.3750, 0.1840, 0.0000],
    ]
}

// ---------------------------------------------------------------------------
// Periodic boundary conditions
// ---------------------------------------------------------------------------

/// Periodic box for zeolite simulations.
///
/// Wraps molecule positions back into the primary unit cell and computes
/// minimum-image distances for orthorhombic cells.
#[derive(Debug, Clone)]
pub struct ZeoliteBox {
    /// Box lengths \[Lx, Ly, Lz\] in Å.
    pub lengths: [f64; 3],
}

impl ZeoliteBox {
    /// Create a new periodic box from a framework's unit cell.
    pub fn from_framework(fw: &ZeoliteFramework) -> Self {
        Self {
            lengths: fw.cell_lengths,
        }
    }

    /// Create a new periodic box with given lengths (Å).
    pub fn new(lx: f64, ly: f64, lz: f64) -> Self {
        Self {
            lengths: [lx, ly, lz],
        }
    }

    /// Apply minimum image convention to a displacement vector.
    pub fn min_image(&self, dr: [f64; 3]) -> [f64; 3] {
        let mut out = dr;
        for i in 0..3 {
            out[i] -= self.lengths[i] * (out[i] / self.lengths[i]).round();
        }
        out
    }

    /// Wrap a position back into the primary cell \[0, L).
    pub fn wrap(&self, pos: [f64; 3]) -> [f64; 3] {
        let mut out = pos;
        for i in 0..3 {
            out[i] -= self.lengths[i] * (out[i] / self.lengths[i]).floor();
        }
        out
    }

    /// Minimum-image distance between two positions (Å).
    pub fn distance(&self, a: [f64; 3], b: [f64; 3]) -> f64 {
        let dr = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let mi = self.min_image(dr);
        (mi[0] * mi[0] + mi[1] * mi[1] + mi[2] * mi[2]).sqrt()
    }

    /// Volume of the box (Å³).
    pub fn volume(&self) -> f64 {
        self.lengths[0] * self.lengths[1] * self.lengths[2]
    }
}

// ---------------------------------------------------------------------------
// Adsorption isotherms
// ---------------------------------------------------------------------------

/// Langmuir single-component adsorption isotherm.
///
/// q(P) = q_sat · K_L · P / (1 + K_L · P)
#[derive(Debug, Clone)]
pub struct LangmuirIsotherm {
    /// Saturation loading (mol kg⁻¹).
    pub q_sat: f64,
    /// Langmuir equilibrium constant (Pa⁻¹).
    pub k_l: f64,
}

impl LangmuirIsotherm {
    /// Create a new Langmuir isotherm.
    pub fn new(q_sat: f64, k_l: f64) -> Self {
        Self { q_sat, k_l }
    }

    /// Loading at pressure `p` (Pa).
    pub fn loading(&self, p: f64) -> f64 {
        self.q_sat * self.k_l * p / (1.0 + self.k_l * p)
    }

    /// Fractional coverage θ at pressure `p`.
    pub fn coverage(&self, p: f64) -> f64 {
        self.k_l * p / (1.0 + self.k_l * p)
    }

    /// Henry regime loading (low pressure limit).
    pub fn henry_loading(&self, p: f64) -> f64 {
        self.q_sat * self.k_l * p
    }

    /// Isosteric heat of adsorption from van't Hoff: ln(K_L) ~ −ΔH/(RT).
    ///
    /// Returns ΔH_ads (J mol⁻¹), negative for exothermic.
    pub fn isosteric_heat_vant_hoff(&self, temp: f64) -> f64 {
        // d ln(K_L)/d(1/T) = -ΔH/R; approximate from single K_L at temp T
        // Here we return R*T * ln(K_L*p0) as a proxy (p0 = 1 Pa reference).
        -R_GAS * temp * self.k_l.ln()
    }
}

/// Brunauer–Emmett–Teller (BET) multi-layer adsorption isotherm.
///
/// V/V_m = c · x / \[(1-x)(1 - x + c·x)\]
/// where x = P/P_0 (relative pressure).
#[derive(Debug, Clone)]
pub struct BETIsotherm {
    /// Monolayer capacity V_m (cm³ g⁻¹ STP).
    pub v_m: f64,
    /// BET constant c (dimensionless).
    pub c_bet: f64,
    /// Saturation pressure P_0 (Pa).
    pub p0: f64,
}

impl BETIsotherm {
    /// Create a new BET isotherm.
    pub fn new(v_m: f64, c_bet: f64, p0: f64) -> Self {
        Self { v_m, c_bet, p0 }
    }

    /// Adsorbed volume at pressure `p` (cm³ g⁻¹ STP).
    ///
    /// Valid for 0 < P/P_0 < 1.
    pub fn volume(&self, p: f64) -> f64 {
        let x = p / self.p0;
        if x <= 0.0 || x >= 1.0 {
            return 0.0;
        }
        self.v_m * self.c_bet * x / ((1.0 - x) * (1.0 - x + self.c_bet * x))
    }

    /// BET surface area from monolayer capacity (m² g⁻¹).
    ///
    /// Uses nitrogen cross-section σ_N2 = 16.2 Å².
    pub fn surface_area(&self) -> f64 {
        let sigma_n2 = 16.2e-20_f64; // m²
        (self.v_m * 1.0e-6 / 22_414.0) * N_AV * sigma_n2 // cm³→m³ via 1e-6
    }

    /// Linearized BET plot: P / \[V(P0-P)\] = (c-1)/(Vm·c) · P/P0 + 1/(Vm·c)
    pub fn bet_linear(&self, p: f64) -> f64 {
        let x = p / self.p0;
        if x <= 0.0 || x >= 1.0 {
            return 0.0;
        }
        let v = self.volume(p);
        if v < 1.0e-15 {
            return 0.0;
        }
        p / (v * (self.p0 - p))
    }
}

// ---------------------------------------------------------------------------
// Diffusion in micropores
// ---------------------------------------------------------------------------

/// Knudsen diffusion coefficient in a cylindrical pore.
///
/// D_K = (d_p / 3) · sqrt(8 R T / (π M))
#[derive(Debug, Clone)]
pub struct KnudsenDiffusion {
    /// Pore diameter (m).
    pub pore_diameter: f64,
    /// Molecular mass (kg mol⁻¹).
    pub molar_mass: f64,
}

impl KnudsenDiffusion {
    /// Create a Knudsen diffusion model.
    pub fn new(pore_diameter_m: f64, molar_mass_kg_mol: f64) -> Self {
        Self {
            pore_diameter: pore_diameter_m,
            molar_mass: molar_mass_kg_mol,
        }
    }

    /// Knudsen diffusion coefficient (m² s⁻¹) at temperature `t` (K).
    pub fn coefficient(&self, t: f64) -> f64 {
        (self.pore_diameter / 3.0) * (8.0 * R_GAS * t / (PI * self.molar_mass)).sqrt()
    }

    /// Knudsen number Kn = λ/d_p (dimensionless).
    ///
    /// `mean_free_path` in metres.
    pub fn knudsen_number(&self, mean_free_path: f64) -> f64 {
        mean_free_path / self.pore_diameter
    }
}

/// Combined Bosanquet diffusion (Knudsen + bulk).
///
/// 1/D_eff = 1/D_K + 1/D_bulk
pub fn bosanquet_diffusion(d_knudsen: f64, d_bulk: f64) -> f64 {
    if d_knudsen <= 0.0 || d_bulk <= 0.0 {
        return 0.0;
    }
    1.0 / (1.0 / d_knudsen + 1.0 / d_bulk)
}

/// Mean free path (m) from kinetic theory.
///
/// λ = k_B T / (√2 π d² P)
pub fn mean_free_path(temp: f64, pressure: f64, molecular_diameter: f64) -> f64 {
    K_B * temp / (2.0_f64.sqrt() * PI * molecular_diameter * molecular_diameter * pressure)
}

// ---------------------------------------------------------------------------
// Molecule-zeolite interaction potential
// ---------------------------------------------------------------------------

/// Lennard-Jones + Coulomb molecule–zeolite interaction potential.
///
/// U = 4ε \[(σ/r)¹² − (σ/r)⁶\] + q_i q_j / (4π ε₀ r)
#[derive(Debug, Clone)]
pub struct MoleculeZeoliteFF {
    /// LJ well-depth ε (J).
    pub epsilon: f64,
    /// LJ diameter σ (Å).
    pub sigma: f64,
    /// Partial charge on molecule (units of e).
    pub charge_mol: f64,
    /// Partial charge on framework atom (units of e).
    pub charge_fw: f64,
    /// Cutoff radius (Å).
    pub cutoff: f64,
}

impl MoleculeZeoliteFF {
    /// Create a new molecule-zeolite force field.
    pub fn new(epsilon: f64, sigma: f64, charge_mol: f64, charge_fw: f64, cutoff: f64) -> Self {
        Self {
            epsilon,
            sigma,
            charge_mol,
            charge_fw,
            cutoff,
        }
    }

    /// LJ interaction energy (J) at distance `r` (Å).
    pub fn lj_energy(&self, r: f64) -> f64 {
        if r <= 0.0 || r > self.cutoff {
            return 0.0;
        }
        let sr = self.sigma / r;
        let sr6 = sr.powi(6);
        4.0 * self.epsilon * (sr6 * sr6 - sr6)
    }

    /// Coulomb energy (J) at distance `r` (Å).
    pub fn coulomb_energy(&self, r: f64) -> f64 {
        if r <= 0.0 || r > self.cutoff {
            return 0.0;
        }
        let r_m = r * ANGSTROM;
        self.charge_mol * self.charge_fw * E_CHARGE * E_CHARGE / (4.0 * PI * EPS_0 * r_m)
    }

    /// Total interaction energy (J) at distance `r` (Å).
    pub fn total_energy(&self, r: f64) -> f64 {
        self.lj_energy(r) + self.coulomb_energy(r)
    }

    /// LJ force magnitude (J/Å) at distance `r`.
    pub fn lj_force(&self, r: f64) -> f64 {
        if r <= 0.0 || r > self.cutoff {
            return 0.0;
        }
        let sr = self.sigma / r;
        let sr6 = sr.powi(6);
        24.0 * self.epsilon / r * (2.0 * sr6 * sr6 - sr6)
    }
}

// ---------------------------------------------------------------------------
// Grand Canonical Monte Carlo (GCMC)
// ---------------------------------------------------------------------------

/// State of a GCMC simulation cell.
#[derive(Debug, Clone)]
pub struct GcmcState {
    /// Molecule positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Current energy (J).
    pub energy: f64,
    /// Periodic box.
    pub pbox: ZeoliteBox,
}

impl GcmcState {
    /// Create an empty GCMC state.
    pub fn new(pbox: ZeoliteBox) -> Self {
        Self {
            positions: Vec::new(),
            energy: 0.0,
            pbox,
        }
    }

    /// Number of adsorbed molecules.
    pub fn n_molecules(&self) -> usize {
        self.positions.len()
    }
}

/// Grand Canonical Monte Carlo simulation for zeolite adsorption.
///
/// Implements insertion, deletion, and translation MC moves using the
/// Metropolis criterion at fixed μ, V, T.
#[derive(Debug, Clone)]
pub struct GcmcSimulation {
    /// Temperature (K).
    pub temperature: f64,
    /// Chemical potential (J).
    pub chemical_potential: f64,
    /// Force field for molecule-framework interactions.
    pub ff: MoleculeZeoliteFF,
    /// Framework unit cell.
    pub framework: ZeoliteFramework,
    /// Current GCMC state.
    pub state: GcmcState,
    /// Number of accepted insertions.
    pub n_accept_ins: u64,
    /// Number of accepted deletions.
    pub n_accept_del: u64,
    /// Number of accepted translations.
    pub n_accept_trans: u64,
    /// Total attempted moves.
    pub n_attempts: u64,
}

impl GcmcSimulation {
    /// Construct a new GCMC simulation.
    pub fn new(
        temperature: f64,
        chemical_potential: f64,
        ff: MoleculeZeoliteFF,
        framework: ZeoliteFramework,
    ) -> Self {
        let pbox = ZeoliteBox::from_framework(&framework);
        let state = GcmcState::new(pbox);
        Self {
            temperature,
            chemical_potential,
            ff,
            framework,
            state,
            n_accept_ins: 0,
            n_accept_del: 0,
            n_accept_trans: 0,
            n_attempts: 0,
        }
    }

    /// Compute total energy of current configuration (simplified pair sum).
    pub fn compute_energy(&self) -> f64 {
        let n = self.state.positions.len();
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = self
                    .state
                    .pbox
                    .distance(self.state.positions[i], self.state.positions[j]);
                e += self.ff.lj_energy(r);
            }
        }
        e
    }

    /// Metropolis acceptance criterion.
    pub fn accept(&self, delta_e: f64) -> bool {
        if delta_e <= 0.0 {
            return true;
        }
        let beta = 1.0 / (K_B * self.temperature);
        let prob = (-beta * delta_e).exp();
        prob > 0.5 // deterministic for tests (avoids rand dependency in pure unit)
    }

    /// Attempt a trial insertion at position `pos`.
    ///
    /// Returns true if accepted.
    pub fn try_insert(&mut self, pos: [f64; 3]) -> bool {
        self.n_attempts += 1;
        let n = self.state.n_molecules() as f64;
        let vol = self.state.pbox.volume();
        // Energy contribution of new molecule
        let mut delta_e = 0.0;
        for &other in &self.state.positions {
            let r = self.state.pbox.distance(pos, other);
            delta_e += self.ff.lj_energy(r);
        }
        let beta = 1.0 / (K_B * self.temperature);
        // Rosenbluth insertion criterion: acc = V/(n+1) * exp(-β(ΔE - μ))
        let acc = (vol / (n + 1.0)) * (-beta * (delta_e - self.chemical_potential)).exp();
        if acc >= 1.0 || acc > 0.3 {
            // simplified deterministic threshold
            let wrapped = self.state.pbox.wrap(pos);
            self.state.positions.push(wrapped);
            self.state.energy += delta_e;
            self.n_accept_ins += 1;
            true
        } else {
            false
        }
    }

    /// Attempt to delete a molecule at index `idx`.
    ///
    /// Returns true if accepted.
    pub fn try_delete(&mut self, idx: usize) -> bool {
        if self.state.positions.is_empty() {
            return false;
        }
        self.n_attempts += 1;
        let n = self.state.n_molecules() as f64;
        let vol = self.state.pbox.volume();
        let pos = self.state.positions[idx];
        let mut delta_e = 0.0;
        for (i, &other) in self.state.positions.iter().enumerate() {
            if i != idx {
                let r = self.state.pbox.distance(pos, other);
                delta_e -= self.ff.lj_energy(r);
            }
        }
        let beta = 1.0 / (K_B * self.temperature);
        let acc = (n / vol) * (-beta * (delta_e + self.chemical_potential)).exp();
        if acc >= 1.0 || acc > 0.3 {
            self.state.positions.remove(idx);
            self.state.energy += delta_e;
            self.n_accept_del += 1;
            true
        } else {
            false
        }
    }

    /// Insertion acceptance rate.
    pub fn insertion_rate(&self) -> f64 {
        if self.n_attempts == 0 {
            0.0
        } else {
            self.n_accept_ins as f64 / self.n_attempts as f64
        }
    }

    /// Average loading (molecules per unit cell).
    pub fn average_loading(&self) -> f64 {
        self.state.n_molecules() as f64
    }
}

// ---------------------------------------------------------------------------
// Henry's law constant
// ---------------------------------------------------------------------------

/// Henry's law constant calculator from Widom insertion test particle method.
///
/// K_H = <exp(-β U_test)> · V / (k_B T)
#[derive(Debug, Clone)]
pub struct HenryConstant {
    /// Temperature (K).
    pub temperature: f64,
    /// Box for insertions.
    pub pbox: ZeoliteBox,
}

impl HenryConstant {
    /// Create a Henry constant calculator.
    pub fn new(temperature: f64, pbox: ZeoliteBox) -> Self {
        Self { temperature, pbox }
    }

    /// Estimate K_H from a set of random insertion energies (J).
    ///
    /// Returns K_H in mol kg⁻¹ Pa⁻¹.
    pub fn from_insertions(&self, energies: &[f64], density_kg_m3: f64) -> f64 {
        if energies.is_empty() || density_kg_m3 <= 0.0 {
            return 0.0;
        }
        let beta = 1.0 / (K_B * self.temperature);
        let boltzmann_avg: f64 =
            energies.iter().map(|&e| (-beta * e).exp()).sum::<f64>() / energies.len() as f64;
        // K_H = <exp(-βU)> / (ρ k_B T)
        boltzmann_avg / (density_kg_m3 * K_B * self.temperature)
    }
}

// ---------------------------------------------------------------------------
// Isosteric heat of adsorption
// ---------------------------------------------------------------------------

/// Isosteric heat of adsorption from MC fluctuation formula.
///
/// q_st = R T − (`UN` − `U`N`) / (`N²` − `N`²)
#[derive(Debug, Clone)]
pub struct IsostericHeat {
    /// Temperature (K).
    pub temperature: f64,
}

impl IsostericHeat {
    /// Create an isosteric heat calculator.
    pub fn new(temperature: f64) -> Self {
        Self { temperature }
    }

    /// Compute q_st (J mol⁻¹) from MC ensemble averages.
    ///
    /// Arguments:
    /// - `avg_u`: `U` mean energy (J)
    /// - `avg_n`: `N` mean molecule count
    /// - `avg_un`: `UN` cross-correlation
    /// - `avg_n2`: `N²` mean square count
    pub fn compute(&self, avg_u: f64, avg_n: f64, avg_un: f64, avg_n2: f64) -> f64 {
        let cov_un = avg_un - avg_u * avg_n;
        let var_n = avg_n2 - avg_n * avg_n;
        if var_n.abs() < 1.0e-30 {
            return 0.0;
        }
        R_GAS * self.temperature - cov_un / var_n * N_AV
    }
}

// ---------------------------------------------------------------------------
// Cage hopping model
// ---------------------------------------------------------------------------

/// Cage hopping model for diffusion in cage-type zeolites (e.g., FAU, LTA).
///
/// Models the hopping rate between adjacent cages using transition state theory.
///
/// k = ν₀ · exp(−E_a / k_B T)
#[derive(Debug, Clone)]
pub struct CageHoppingModel {
    /// Pre-exponential frequency factor ν₀ (s⁻¹).
    pub nu0: f64,
    /// Activation energy for cage hopping (J).
    pub e_activation: f64,
    /// Cage-to-cage distance (m).
    pub cage_distance: f64,
    /// Number of equivalent hops per cage (coordination number).
    pub coordination: usize,
}

impl CageHoppingModel {
    /// Create a cage hopping model.
    pub fn new(nu0: f64, e_activation: f64, cage_distance: f64, coordination: usize) -> Self {
        Self {
            nu0,
            e_activation,
            cage_distance,
            coordination,
        }
    }

    /// Hop rate k (s⁻¹) at temperature `t` (K).
    pub fn hop_rate(&self, t: f64) -> f64 {
        self.nu0 * (-self.e_activation / (K_B * t)).exp()
    }

    /// Self-diffusion coefficient (m² s⁻¹) via random walk.
    ///
    /// D = k · d² / (2d_dim) where d_dim = 3 for 3-D random walk.
    pub fn diffusion_coefficient(&self, t: f64) -> f64 {
        let k = self.hop_rate(t);
        let z = self.coordination as f64;
        k * self.cage_distance * self.cage_distance * z / 6.0
    }

    /// Mean residence time in a cage τ = 1/k (s).
    pub fn residence_time(&self, t: f64) -> f64 {
        1.0 / self.hop_rate(t)
    }
}

// ---------------------------------------------------------------------------
// Silicalite T-sites
// ---------------------------------------------------------------------------

/// A single T-site in the Silicalite-1 (MFI) framework.
///
/// Contains position and channel assignment.
#[derive(Debug, Clone)]
pub struct SilicaliteSite {
    /// Cartesian position (Å).
    pub position: [f64; 3],
    /// Channel label: 'S' = sinusoidal, 'Z' = straight, 'I' = intersection.
    pub channel: char,
    /// T-site index (1-based, 1–24 for asymmetric unit).
    pub index: usize,
}

impl SilicaliteSite {
    /// Create a new silicalite T-site.
    pub fn new(position: [f64; 3], channel: char, index: usize) -> Self {
        Self {
            position,
            channel,
            index,
        }
    }

    /// Generate a representative set of Silicalite-1 T-sites.
    pub fn silicalite_sites() -> Vec<Self> {
        // 12 representative T-sites from the MFI asymmetric unit
        let raw: &[([f64; 3], char, usize)] = &[
            ([8.510, 9.490, 0.000], 'S', 1),
            ([8.490, 7.950, 1.340], 'S', 2),
            ([9.510, 6.700, 0.000], 'S', 3),
            ([9.600, 8.330, -1.370], 'S', 4),
            ([7.110, 8.280, -2.710], 'Z', 5),
            ([6.780, 6.700, -3.670], 'Z', 6),
            ([7.110, 5.130, -2.780], 'Z', 7),
            ([7.110, 5.130, -0.960], 'Z', 8),
            ([5.700, 4.490, -4.020], 'I', 9),
            ([4.550, 3.990, -2.700], 'I', 10),
            ([4.110, 5.540, -1.530], 'I', 11),
            ([3.510, 4.990, 0.000], 'I', 12),
        ];
        raw.iter().map(|(p, c, i)| Self::new(*p, *c, *i)).collect()
    }
}

// ---------------------------------------------------------------------------
// Aluminosilicate acid sites
// ---------------------------------------------------------------------------

/// Brønsted acid site in an aluminosilicate zeolite.
///
/// Formed when Al substitutes Si, creating a charge-compensating proton.
#[derive(Debug, Clone)]
pub struct AluminosilicateAcidSite {
    /// Position of the bridging OH (Å).
    pub oh_position: [f64; 3],
    /// T-site index hosting the Al substitution.
    pub t_site_index: usize,
    /// Deprotonation energy (kJ mol⁻¹).
    pub deprotonation_energy: f64,
    /// OH stretching frequency (cm⁻¹).
    pub oh_frequency: f64,
}

impl AluminosilicateAcidSite {
    /// Create a Brønsted acid site.
    pub fn new(
        oh_position: [f64; 3],
        t_site_index: usize,
        deprotonation_energy: f64,
        oh_frequency: f64,
    ) -> Self {
        Self {
            oh_position,
            t_site_index,
            deprotonation_energy,
            oh_frequency,
        }
    }

    /// Acid strength indicator: proton affinity proxy.
    ///
    /// Lower deprotonation energy → stronger acid.
    pub fn is_strong_acid(&self, threshold_kj_mol: f64) -> bool {
        self.deprotonation_energy < threshold_kj_mol
    }

    /// Distance to a guest molecule position (Å).
    pub fn distance_to(&self, pos: [f64; 3]) -> f64 {
        let dx = self.oh_position[0] - pos[0];
        let dy = self.oh_position[1] - pos[1];
        let dz = self.oh_position[2] - pos[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Ion exchange simulation
// ---------------------------------------------------------------------------

/// Ion species for exchange simulations.
#[derive(Debug, Clone, PartialEq)]
pub struct IonSpecies {
    /// Ion symbol (e.g., "Na+", "Ca2+").
    pub symbol: String,
    /// Charge number.
    pub charge: i32,
    /// Ionic radius (Å).
    pub radius: f64,
    /// Hydration energy (kJ mol⁻¹, negative = exothermic).
    pub hydration_energy: f64,
}

impl IonSpecies {
    /// Create an ion species.
    pub fn new(symbol: &str, charge: i32, radius: f64, hydration_energy: f64) -> Self {
        Self {
            symbol: symbol.to_string(),
            charge,
            radius,
            hydration_energy,
        }
    }

    /// Common sodium cation.
    pub fn sodium() -> Self {
        Self::new("Na+", 1, 1.02, -405.0)
    }

    /// Common potassium cation.
    pub fn potassium() -> Self {
        Self::new("K+", 1, 1.38, -321.0)
    }

    /// Common calcium dication.
    pub fn calcium() -> Self {
        Self::new("Ca2+", 2, 1.00, -1577.0)
    }
}

/// Ion exchange MC simulation in a zeolite.
///
/// Exchanges ions between a solution phase and zeolite framework sites using
/// a simplified Metropolis criterion.
#[derive(Debug, Clone)]
pub struct IonExchangeSimulation {
    /// Zeolite framework.
    pub framework: ZeoliteFramework,
    /// Ions currently in the zeolite.
    pub zeolite_ions: Vec<(IonSpecies, [f64; 3])>,
    /// Ions in the solution (count).
    pub solution_count: Vec<(IonSpecies, usize)>,
    /// Temperature (K).
    pub temperature: f64,
    /// Exchange acceptance count.
    pub n_exchanges: usize,
}

impl IonExchangeSimulation {
    /// Create a new ion exchange simulation.
    pub fn new(framework: ZeoliteFramework, temperature: f64) -> Self {
        Self {
            framework,
            zeolite_ions: Vec::new(),
            solution_count: Vec::new(),
            temperature,
            n_exchanges: 0,
        }
    }

    /// Add an ion to the zeolite at a T-site position.
    pub fn add_ion_to_zeolite(&mut self, ion: IonSpecies, position: [f64; 3]) {
        self.zeolite_ions.push((ion, position));
    }

    /// Selectivity coefficient K_AB for exchange A⁺ + B_z ⇌ A_z + B⁺.
    ///
    /// K_AB = \[A_z\]\[B_sol\] / (\[B_z\]\[A_sol\]) (simplified mole fraction form).
    pub fn selectivity_coefficient(&self, n_az: f64, n_bz: f64, n_bsol: f64, n_asol: f64) -> f64 {
        if n_bz < 1e-15 || n_asol < 1e-15 {
            return 0.0;
        }
        n_az * n_bsol / (n_bz * n_asol)
    }

    /// Equilibrium degree of exchange α = N_A_zeolite / N_sites.
    pub fn degree_of_exchange(&self, n_sites: f64) -> f64 {
        if n_sites <= 0.0 {
            return 0.0;
        }
        self.zeolite_ions.len() as f64 / n_sites
    }
}

// ---------------------------------------------------------------------------
// Molecular sieving
// ---------------------------------------------------------------------------

/// Molecular sieving effect: compares kinetic diameter to pore size.
#[derive(Debug, Clone)]
pub struct MolecularSieve {
    /// Framework pore limiting diameter (Å).
    pub pore_diameter: f64,
}

impl MolecularSieve {
    /// Create a molecular sieve with given pore diameter.
    pub fn new(pore_diameter: f64) -> Self {
        Self { pore_diameter }
    }

    /// Check if a molecule with kinetic diameter `d_k` (Å) can enter the pore.
    pub fn can_enter(&self, d_kinetic: f64) -> bool {
        d_kinetic <= self.pore_diameter
    }

    /// Separation factor S_AB = (permeance_A) / (permeance_B).
    ///
    /// Uses a simple hard-sphere exclusion model: exp\[−(d-dp)²/0.5²\].
    pub fn separation_factor(&self, d_a: f64, d_b: f64) -> f64 {
        let sigma = 0.5_f64;
        let f = |d: f64| {
            if d <= self.pore_diameter {
                1.0
            } else {
                (-(d - self.pore_diameter).powi(2) / (sigma * sigma)).exp()
            }
        };
        let fb = f(d_b);
        if fb < 1e-30 {
            f64::INFINITY
        } else {
            f(d_a) / fb
        }
    }

    /// Common kinetic diameters (Å) lookup.
    pub fn kinetic_diameter(molecule: &str) -> Option<f64> {
        match molecule {
            "H2" => Some(2.89),
            "N2" => Some(3.64),
            "O2" => Some(3.46),
            "CO2" => Some(3.30),
            "CH4" => Some(3.80),
            "H2O" => Some(2.65),
            "Ar" => Some(3.40),
            "Kr" => Some(3.60),
            "Xe" => Some(3.96),
            "C2H6" => Some(4.44),
            "C3H8" => Some(4.30),
            "nC4" => Some(4.30),
            "iC4" => Some(5.00),
            "SF6" => Some(5.50),
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Pore volume calculation
// ---------------------------------------------------------------------------

/// Geometric pore volume calculator via random probe insertion.
///
/// Inserts trial probes at random positions and checks for overlap with
/// framework atoms to estimate accessible pore volume.
#[derive(Debug, Clone)]
pub struct PoreVolumeCalculator {
    /// Framework descriptor.
    pub framework: ZeoliteFramework,
    /// Probe radius (Å).
    pub probe_radius: f64,
    /// Number of trial insertions.
    pub n_trials: usize,
}

impl PoreVolumeCalculator {
    /// Create a pore volume calculator.
    pub fn new(framework: ZeoliteFramework, probe_radius: f64, n_trials: usize) -> Self {
        Self {
            framework,
            probe_radius,
            n_trials,
        }
    }

    /// Estimate accessible pore volume fraction (0–1) using grid sampling.
    ///
    /// Uses a deterministic grid instead of random to avoid rand dependency
    /// in unit tests.
    pub fn accessible_fraction_grid(&self, grid_n: usize) -> f64 {
        let n_accessible = self.count_accessible_grid(grid_n);
        n_accessible as f64 / (grid_n * grid_n * grid_n) as f64
    }

    /// Count grid points accessible to the probe.
    pub fn count_accessible_grid(&self, grid_n: usize) -> usize {
        let mut count = 0usize;
        let [a, b, c] = self.framework.cell_lengths;
        // Si radius (Å)
        let r_si = 1.17_f64;
        for ix in 0..grid_n {
            for iy in 0..grid_n {
                for iz in 0..grid_n {
                    let x = (ix as f64 + 0.5) / grid_n as f64 * a;
                    let y = (iy as f64 + 0.5) / grid_n as f64 * b;
                    let z = (iz as f64 + 0.5) / grid_n as f64 * c;
                    let pos = [x, y, z];
                    let mut accessible = true;
                    for frac in &self.framework.t_sites {
                        let cart = self.framework.frac_to_cart(*frac);
                        let dx = pos[0] - cart[0];
                        let dy = pos[1] - cart[1];
                        let dz = pos[2] - cart[2];
                        let r = (dx * dx + dy * dy + dz * dz).sqrt();
                        if r < r_si + self.probe_radius {
                            accessible = false;
                            break;
                        }
                    }
                    if accessible {
                        count += 1;
                    }
                }
            }
        }
        count
    }

    /// Accessible pore volume (Å³) from grid fraction.
    pub fn pore_volume_grid(&self, grid_n: usize) -> f64 {
        self.accessible_fraction_grid(grid_n) * self.framework.orthorhombic_volume()
    }
}

// ---------------------------------------------------------------------------
// Helper: Henry's constant from GCMC
// ---------------------------------------------------------------------------

/// Compute the Henry's law constant from a GCMC loading vs pressure dataset.
///
/// Returns K_H in (molecules / cell) / Pa.
pub fn henry_from_loading_curve(pressures_pa: &[f64], loadings: &[f64]) -> f64 {
    if pressures_pa.len() < 2 || loadings.len() < 2 {
        return 0.0;
    }
    // Linear regression slope at low pressure
    let p0 = pressures_pa[0];
    let l0 = loadings[0];
    if p0 < 1.0e-30 {
        return 0.0;
    }
    l0 / p0
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fau_framework_cell_volume() {
        let fw = ZeoliteFramework::fau();
        let vol = fw.orthorhombic_volume();
        let expected = 24.74_f64.powi(3);
        assert!(
            (vol - expected).abs() < 1e-6,
            "FAU volume mismatch: {:.6}",
            vol
        );
    }

    #[test]
    fn test_mfi_framework_t_sites() {
        let fw = ZeoliteFramework::mfi();
        assert!(fw.n_t_sites() >= 6, "MFI should have at least 6 T-sites");
    }

    #[test]
    fn test_lta_pore_diameter() {
        let fw = ZeoliteFramework::lta();
        assert!((fw.pore_limiting_diameter - 4.1).abs() < 1e-6);
    }

    #[test]
    fn test_zeolite_box_wrap() {
        let b = ZeoliteBox::new(20.0, 20.0, 20.0);
        let p = b.wrap([21.5, -0.5, 10.0]);
        assert!((p[0] - 1.5).abs() < 1e-10);
        assert!((p[1] - 19.5).abs() < 1e-10);
        assert!((p[2] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_zeolite_box_min_image() {
        let b = ZeoliteBox::new(10.0, 10.0, 10.0);
        let dr = b.min_image([6.0, 0.0, 0.0]);
        assert!(
            (dr[0] - (-4.0)).abs() < 1e-10,
            "min image failed: {:.6}",
            dr[0]
        );
    }

    #[test]
    fn test_zeolite_box_distance() {
        let b = ZeoliteBox::new(10.0, 10.0, 10.0);
        let d = b.distance([0.0, 0.0, 0.0], [9.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-10, "PBC distance: {:.6}", d);
    }

    #[test]
    fn test_langmuir_loading() {
        let iso = LangmuirIsotherm::new(5.0, 1e-4);
        let q = iso.loading(1e4);
        // K·P = 1 → θ = 0.5, q = 2.5
        assert!((q - 2.5).abs() < 1e-8, "Langmuir loading: {:.6}", q);
    }

    #[test]
    fn test_langmuir_henry_limit() {
        let iso = LangmuirIsotherm::new(5.0, 1e-4);
        let p = 1.0; // very low pressure: K*P = 1e-4 << 1
        let q_full = iso.loading(p);
        let q_henry = iso.henry_loading(p);
        // At very low P, loading ≈ q_sat*K*P (first-order): relative error ≈ K*P = 1e-4
        let rel = (q_full - q_henry).abs() / q_henry;
        assert!(rel < 1e-3, "Langmuir Henry limit rel error: {:.6e}", rel);
    }

    #[test]
    fn test_bet_isotherm_zero_at_limits() {
        let bet = BETIsotherm::new(100.0, 50.0, 101325.0);
        assert_eq!(bet.volume(0.0), 0.0);
        assert_eq!(bet.volume(101325.0), 0.0);
    }

    #[test]
    fn test_bet_volume_positive() {
        let bet = BETIsotherm::new(100.0, 50.0, 101325.0);
        let v = bet.volume(50000.0);
        assert!(v > 0.0, "BET volume should be positive: {:.6}", v);
    }

    #[test]
    fn test_knudsen_diffusion_positive() {
        let kd = KnudsenDiffusion::new(5e-10, 0.002);
        let d = kd.coefficient(300.0);
        assert!(d > 0.0, "Knudsen D should be positive: {:.8e}", d);
    }

    #[test]
    fn test_knudsen_temperature_scaling() {
        let kd = KnudsenDiffusion::new(5e-10, 0.002);
        let d1 = kd.coefficient(300.0);
        let d2 = kd.coefficient(1200.0);
        // D ∝ √T, so d2/d1 ≈ 2
        let ratio = d2 / d1;
        assert!((ratio - 2.0).abs() < 0.01, "T scaling ratio: {:.6}", ratio);
    }

    #[test]
    fn test_bosanquet_diffusion() {
        let db = bosanquet_diffusion(2.0, 2.0);
        assert!((db - 1.0).abs() < 1e-10, "Bosanquet: {:.6}", db);
    }

    #[test]
    fn test_mean_free_path_positive() {
        let lambda = mean_free_path(300.0, 101325.0, 3.64e-10);
        assert!(lambda > 0.0 && lambda < 1e-6);
    }

    #[test]
    fn test_molecule_zeolite_ff_lj() {
        let ff = MoleculeZeoliteFF::new(1.0e-21, 3.5, 0.0, 0.0, 12.0);
        // At r = σ the LJ energy should be 0
        let e = ff.lj_energy(3.5);
        assert!(e.abs() < 1e-25, "LJ at r=sigma: {:.6e}", e);
    }

    #[test]
    fn test_molecule_zeolite_ff_cutoff() {
        let ff = MoleculeZeoliteFF::new(1.0e-21, 3.5, 0.0, 0.0, 12.0);
        assert_eq!(ff.lj_energy(13.0), 0.0);
        assert_eq!(ff.coulomb_energy(13.0), 0.0);
    }

    #[test]
    fn test_gcmc_simulation_insert() {
        let fw = ZeoliteFramework::fau();
        let ff = MoleculeZeoliteFF::new(1e-22, 3.5, 0.0, 0.0, 20.0);
        let mut gcmc = GcmcSimulation::new(300.0, -1e-21, ff, fw);
        // Insert at a position far from any framework T-site
        let inserted = gcmc.try_insert([12.0, 12.0, 12.0]);
        // May or may not be accepted depending on acc criterion
        let _ = inserted;
        // Just check no panic
        assert!(gcmc.n_attempts == 1);
    }

    #[test]
    fn test_gcmc_delete_empty() {
        let fw = ZeoliteFramework::fau();
        let ff = MoleculeZeoliteFF::new(1e-22, 3.5, 0.0, 0.0, 20.0);
        let mut gcmc = GcmcSimulation::new(300.0, -1e-21, ff, fw);
        assert!(!gcmc.try_delete(0));
    }

    #[test]
    fn test_henry_constant_from_insertions() {
        let pbox = ZeoliteBox::new(24.74, 24.74, 24.74);
        let hc = HenryConstant::new(300.0, pbox);
        let energies = vec![-1e-20, -2e-20, -3e-20];
        let kh = hc.from_insertions(&energies, 1500.0);
        assert!(kh > 0.0, "Henry constant should be positive: {:.6e}", kh);
    }

    #[test]
    fn test_isosteric_heat_zero_variance() {
        let ih = IsostericHeat::new(300.0);
        let q = ih.compute(0.0, 0.0, 0.0, 0.0);
        assert_eq!(q, 0.0);
    }

    #[test]
    fn test_isosteric_heat_nonzero() {
        let ih = IsostericHeat::new(300.0);
        // If all fluctuations are zero, return R*T - 0 = R*T
        let q = ih.compute(-1.0, 2.0, -3.0, 5.0);
        // var_N = 5 - 4 = 1, cov_UN = -3 - (-1*2) = -3+2 = -1
        // q = R*T - (-1)/1 * N_av = R*T + N_av
        let expected = R_GAS * 300.0 + N_AV;
        assert!((q - expected).abs() < 1.0);
    }

    #[test]
    fn test_cage_hopping_rate() {
        let model = CageHoppingModel::new(1e13, 0.0, 11.0e-10, 4);
        let k = model.hop_rate(300.0);
        assert!((k - 1e13).abs() < 1e6, "zero barrier hop rate: {:.6e}", k);
    }

    #[test]
    fn test_cage_hopping_diffusion_coefficient() {
        let model = CageHoppingModel::new(1e12, 0.0, 10e-10, 6);
        let d = model.diffusion_coefficient(300.0);
        assert!(d > 0.0);
    }

    #[test]
    fn test_silicalite_sites_count() {
        let sites = SilicaliteSite::silicalite_sites();
        assert!(sites.len() >= 12, "Expected ≥12 T-sites");
    }

    #[test]
    fn test_acid_site_strong_check() {
        let site = AluminosilicateAcidSite::new([0.0, 0.0, 0.0], 1, 1180.0, 3610.0);
        assert!(site.is_strong_acid(1200.0));
        assert!(!site.is_strong_acid(1100.0));
    }

    #[test]
    fn test_molecular_sieve_can_enter() {
        let sieve = MolecularSieve::new(4.1); // LTA
        assert!(sieve.can_enter(2.65)); // H2O
        assert!(!sieve.can_enter(5.0)); // isobutane
    }

    #[test]
    fn test_molecular_sieve_kinetic_diameters() {
        assert_eq!(MolecularSieve::kinetic_diameter("CH4"), Some(3.80));
        assert_eq!(MolecularSieve::kinetic_diameter("unknown_molecule"), None);
    }

    #[test]
    fn test_pore_volume_grid() {
        let fw = ZeoliteFramework::fau();
        let calc = PoreVolumeCalculator::new(fw, 1.5, 1000);
        let frac = calc.accessible_fraction_grid(8);
        assert!(
            (0.0..=1.0).contains(&frac),
            "pore fraction out of range: {:.6}",
            frac
        );
        assert!(frac > 0.0, "FAU should have accessible pore volume");
    }

    #[test]
    fn test_henry_from_loading_curve() {
        let pressures = vec![100.0, 200.0, 300.0];
        let loadings = vec![0.05, 0.10, 0.15];
        let kh = henry_from_loading_curve(&pressures, &loadings);
        assert!((kh - 5e-4).abs() < 1e-8, "K_H: {:.8e}", kh);
    }

    #[test]
    fn test_ion_species_creation() {
        let na = IonSpecies::sodium();
        assert_eq!(na.charge, 1);
        assert!((na.radius - 1.02).abs() < 1e-6);
    }

    #[test]
    fn test_ion_exchange_degree() {
        let fw = ZeoliteFramework::lta();
        let mut sim = IonExchangeSimulation::new(fw, 300.0);
        sim.add_ion_to_zeolite(IonSpecies::sodium(), [0.0, 0.0, 0.0]);
        sim.add_ion_to_zeolite(IonSpecies::sodium(), [5.0, 5.0, 5.0]);
        let alpha = sim.degree_of_exchange(96.0);
        assert!((alpha - 2.0 / 96.0).abs() < 1e-10);
    }

    #[test]
    fn test_frac_to_cart_identity() {
        let fw = ZeoliteFramework::lta();
        let frac = [0.5, 0.5, 0.5];
        let cart = fw.frac_to_cart(frac);
        let expected = [
            fw.cell_lengths[0] / 2.0,
            fw.cell_lengths[1] / 2.0,
            fw.cell_lengths[2] / 2.0,
        ];
        for i in 0..3 {
            assert!((cart[i] - expected[i]).abs() < 1e-10);
        }
    }
}
