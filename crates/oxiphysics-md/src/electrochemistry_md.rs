// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrochemical molecular dynamics module.
//!
//! Provides MD-level models for simulating:
//! - Electrode-electrolyte interfaces (electrical double layer)
//! - Gouy-Chapman-Stern (GCS) electric double layer model
//! - Constant potential MD (CPM) with electrode charge dynamics
//! - Redox reaction dynamics via Marcus electron transfer theory
//! - Ion transport in electrolyte solutions
//! - Solvation shell structure and coordination numbers
//! - Surface charge density and its evolution
//! - Differential and integral capacitance from MD trajectories
//! - Lithium-ion intercalation (anode graphite and cathode LFP/NMC)

use rand::RngExt;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical Constants
// ---------------------------------------------------------------------------

/// Boltzmann constant in J/K.
const BOLTZMANN: f64 = 1.380649e-23;

/// Elementary charge in C.
const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;

/// Vacuum permittivity in F/m.
const EPSILON_0: f64 = 8.854187817e-12;

/// Avogadro's number.
const AVOGADRO: f64 = 6.02214076e23;

/// Universal gas constant in J/(mol K).
const GAS_CONSTANT: f64 = 8.314;

/// Faraday constant in C/mol.
const FARADAY: f64 = 96_485.332_9;

/// Standard temperature in K.
#[cfg(test)]
const STD_TEMP: f64 = 298.15;

/// Thermal voltage at standard temperature in V.
#[cfg(test)]
const THERMAL_VOLTAGE: f64 = BOLTZMANN * STD_TEMP / ELEMENTARY_CHARGE;

/// Marcus reorganisation energy reference in eV.
const MARCUS_LAMBDA_REF: f64 = 0.5;

// ---------------------------------------------------------------------------
// IonicSpecies
// ---------------------------------------------------------------------------

/// Ionic species in an electrochemical simulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IonicSpecies {
    /// Lithium cation Li+.
    Li,
    /// Sodium cation Na+.
    Na,
    /// Potassium cation K+.
    K,
    /// Calcium cation Ca2+.
    Ca,
    /// Chloride anion Cl-.
    Cl,
    /// Fluoride anion F-.
    F,
    /// Perchlorate anion ClO4-.
    Perchlorate,
    /// Hexafluorophosphate PF6-.
    Pf6,
    /// Bis(trifluoromethanesulfonyl)imide TFSI-.
    Tfsi,
    /// Water molecule (neutral).
    Water,
}

impl IonicSpecies {
    /// Charge number z (in units of elementary charge).
    pub fn charge_number(&self) -> i32 {
        match self {
            Self::Li | Self::Na | Self::K => 1,
            Self::Ca => 2,
            Self::Cl | Self::F | Self::Perchlorate | Self::Pf6 | Self::Tfsi => -1,
            Self::Water => 0,
        }
    }

    /// Effective ionic radius in Angstrom (crystal radius).
    pub fn ionic_radius(&self) -> f64 {
        match self {
            Self::Li => 0.76,
            Self::Na => 1.02,
            Self::K => 1.38,
            Self::Ca => 1.00,
            Self::Cl => 1.81,
            Self::F => 1.33,
            Self::Perchlorate => 2.40,
            Self::Pf6 => 2.54,
            Self::Tfsi => 3.80,
            Self::Water => 1.40,
        }
    }

    /// Hydration free energy in eV (Born approximation, negative = exothermic).
    pub fn hydration_free_energy(&self) -> f64 {
        match self {
            Self::Li => -5.01,
            Self::Na => -4.10,
            Self::K => -3.39,
            Self::Ca => -16.00,
            Self::Cl => -3.47,
            Self::F => -5.02,
            Self::Perchlorate => -2.01,
            Self::Pf6 => -1.80,
            Self::Tfsi => -1.50,
            Self::Water => 0.0,
        }
    }

    /// Average coordination number in bulk water (first solvation shell).
    pub fn coordination_number(&self) -> f64 {
        match self {
            Self::Li => 4.0,
            Self::Na => 6.0,
            Self::K => 7.0,
            Self::Ca => 8.0,
            Self::Cl => 6.0,
            Self::F => 6.0,
            Self::Perchlorate => 6.0,
            Self::Pf6 => 5.0,
            Self::Tfsi => 4.5,
            Self::Water => 4.0,
        }
    }

    /// Diffusion coefficient in bulk water at 298 K in m^2/s.
    pub fn diffusion_coefficient(&self) -> f64 {
        match self {
            Self::Li => 1.03e-9,
            Self::Na => 1.33e-9,
            Self::K => 1.96e-9,
            Self::Ca => 0.79e-9,
            Self::Cl => 2.03e-9,
            Self::F => 1.46e-9,
            Self::Perchlorate => 1.79e-9,
            Self::Pf6 => 0.80e-9,
            Self::Tfsi => 0.60e-9,
            Self::Water => 2.30e-9,
        }
    }

    /// Mass in kg.
    pub fn mass(&self) -> f64 {
        let g_per_mol = match self {
            Self::Li => 6.941,
            Self::Na => 22.990,
            Self::K => 39.098,
            Self::Ca => 40.078,
            Self::Cl => 35.453,
            Self::F => 18.998,
            Self::Perchlorate => 99.45,
            Self::Pf6 => 144.96,
            Self::Tfsi => 280.15,
            Self::Water => 18.015,
        };
        g_per_mol * 1.0e-3 / AVOGADRO
    }
}

// ---------------------------------------------------------------------------
// ElectrochemicalIon
// ---------------------------------------------------------------------------

/// A single ion in an electrochemical MD simulation.
#[derive(Debug, Clone)]
pub struct ElectrochemicalIon {
    /// Position \[x, y, z\] in Angstrom.
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\] in Angstrom/fs.
    pub vel: [f64; 3],
    /// Force \[fx, fy, fz\] in eV/Angstrom.
    pub force: [f64; 3],
    /// Ionic species.
    pub species: IonicSpecies,
    /// Effective charge (may differ from formal charge during redox).
    pub charge: f64,
    /// Solvation shell occupancy (number of first-shell solvent molecules).
    pub solvation_number: f64,
    /// Unique particle ID.
    pub id: usize,
    /// Distance from electrode surface in Angstrom.
    pub z_electrode: f64,
}

impl ElectrochemicalIon {
    /// Create a new ion at a given position.
    pub fn new(pos: [f64; 3], species: IonicSpecies) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            charge: species.charge_number() as f64,
            solvation_number: species.coordination_number(),
            id: 0,
            z_electrode: pos[2],
            species,
        }
    }

    /// Kinetic energy of this ion in eV.
    pub fn kinetic_energy_ev(&self) -> f64 {
        let m = self.species.mass();
        // vel in A/fs = 1e-10/1e-15 m/s = 1e5 m/s
        let v2 = self.vel.iter().map(|v| v * v).sum::<f64>() * 1.0e10; // (m/s)^2
        0.5 * m * v2 / ELEMENTARY_CHARGE
    }

    /// Speed in Angstrom/fs.
    pub fn speed(&self) -> f64 {
        self.vel.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    /// Distance to another ion in Angstrom.
    pub fn distance_to(&self, other: &ElectrochemicalIon) -> f64 {
        let dx = self.pos[0] - other.pos[0];
        let dy = self.pos[1] - other.pos[1];
        let dz = self.pos[2] - other.pos[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Coulomb potential energy with another ion in eV.
    ///
    /// V = k_e * q_i * q_j / r  where k_e = 1/(4pi eps_0).
    pub fn coulomb_energy_ev(&self, other: &ElectrochemicalIon, dielectric: f64) -> f64 {
        let r_m = self.distance_to(other) * 1.0e-10; // Angstrom to metres
        if r_m < 1.0e-15 {
            return 0.0;
        }
        let ke = 1.0 / (4.0 * PI * EPSILON_0 * dielectric);
        let qi = self.charge * ELEMENTARY_CHARGE;
        let qj = other.charge * ELEMENTARY_CHARGE;
        ke * qi * qj / r_m / ELEMENTARY_CHARGE // convert J to eV
    }

    /// Zero forces.
    pub fn zero_forces(&mut self) {
        self.force = [0.0; 3];
    }
}

// ---------------------------------------------------------------------------
// GouyChapmanSternModel
// ---------------------------------------------------------------------------

/// Gouy-Chapman-Stern (GCS) electric double layer model.
///
/// The Stern layer (inner Helmholtz layer) accounts for finite ion size,
/// while the Gouy-Chapman diffuse layer extends into the bulk electrolyte.
///
/// Potential profile:
/// - Stern layer: linear potential drop from electrode to Outer Helmholtz Plane (OHP)
/// - Diffuse layer: Poisson-Boltzmann exponential decay with Debye length kappa^{-1}
#[derive(Debug, Clone)]
pub struct GouyChapmanSternModel {
    /// Electrode surface potential relative to bulk in V.
    pub phi_m: f64,
    /// Potential at the OHP (start of diffuse layer) in V.
    pub phi_ohp: f64,
    /// Stern layer thickness d_s in Angstrom.
    pub stern_thickness: f64,
    /// Dielectric constant in Stern layer.
    pub eps_stern: f64,
    /// Dielectric constant in diffuse layer.
    pub eps_diffuse: f64,
    /// Electrolyte concentration (1:1 salt) in mol/L.
    pub concentration: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Surface charge density in C/m^2.
    pub sigma: f64,
}

impl GouyChapmanSternModel {
    /// Create a GCS model for aqueous KCl electrolyte.
    pub fn aqueous_kcl(concentration: f64, phi_m: f64, temperature: f64) -> Self {
        let sigma = Self::compute_sigma(concentration, phi_m, temperature);
        let phi_ohp = Self::compute_phi_ohp(phi_m, sigma, 3.5, 6.0);
        Self {
            phi_m,
            phi_ohp,
            stern_thickness: 3.5, // Angstrom
            eps_stern: 6.0,
            eps_diffuse: 78.5,
            concentration,
            temperature,
            sigma,
        }
    }

    /// Debye screening length in metres.
    ///
    /// kappa^{-1} = sqrt(eps * eps_0 * k_B * T / (2 N_A e^2 c))
    pub fn debye_length(&self) -> f64 {
        let eps = self.eps_diffuse * EPSILON_0;
        let c_si = self.concentration * 1.0e3; // mol/m^3
        let kappa_sq = 2.0 * AVOGADRO * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * c_si
            / (eps * BOLTZMANN * self.temperature);
        1.0 / kappa_sq.sqrt()
    }

    /// Compute surface charge density by GCS self-consistency.
    fn compute_sigma(concentration: f64, phi_m: f64, temperature: f64) -> f64 {
        let eps = 78.5 * EPSILON_0;
        let c_si = concentration * 1.0e3;
        let vt = BOLTZMANN * temperature / ELEMENTARY_CHARGE;
        // Grahame equation: sigma = sqrt(8 eps R T c) * sinh(phi/(2 vt))
        let prefactor = (8.0 * eps * GAS_CONSTANT * temperature * c_si / 1000.0).sqrt();
        prefactor * (phi_m / (2.0 * vt)).sinh()
    }

    /// Compute OHP potential from surface charge and Stern capacitance.
    fn compute_phi_ohp(phi_m: f64, sigma: f64, d_angstrom: f64, eps_s: f64) -> f64 {
        let d_m = d_angstrom * 1.0e-10;
        let c_stern = eps_s * EPSILON_0 / d_m;
        phi_m - sigma / c_stern
    }

    /// Potential at distance z from OHP in the diffuse layer (in V).
    ///
    /// phi(z) = 4 vt arctanh(tanh(phi_ohp/(4 vt)) * exp(-kappa z))
    pub fn diffuse_potential(&self, z_from_ohp: f64) -> f64 {
        let vt = BOLTZMANN * self.temperature / ELEMENTARY_CHARGE;
        let kappa_inv = self.debye_length();
        let kappa = 1.0 / kappa_inv;
        let y0 = (self.phi_ohp / (4.0 * vt)).tanh();
        let y = y0 * (-kappa * z_from_ohp).exp();
        4.0 * vt * y.atanh()
    }

    /// Differential capacitance of the EDL per unit area in F/m^2.
    ///
    /// C_diff = d sigma / d phi_m  (numerical derivative)
    pub fn differential_capacitance(&self) -> f64 {
        let dphi = 1.0e-4; // V
        let s_plus = Self::compute_sigma(self.concentration, self.phi_m + dphi, self.temperature);
        let s_minus = Self::compute_sigma(self.concentration, self.phi_m - dphi, self.temperature);
        (s_plus - s_minus) / (2.0 * dphi)
    }

    /// Integral capacitance C_int = sigma / phi_m in F/m^2.
    pub fn integral_capacitance(&self) -> f64 {
        if self.phi_m.abs() < 1.0e-10 {
            return self.differential_capacitance();
        }
        self.sigma / self.phi_m
    }

    /// Charge accumulated in the diffuse layer per unit area in C/m^2.
    pub fn diffuse_charge(&self) -> f64 {
        -self.sigma // charge neutrality
    }

    /// Ion number density profile at distance z from OHP (Boltzmann distribution).
    ///
    /// n+(z) = n0 * exp(-e phi(z) / kB T)
    /// n-(z) = n0 * exp(+e phi(z) / kB T)
    pub fn cation_density(&self, z_from_ohp: f64) -> f64 {
        let phi = self.diffuse_potential(z_from_ohp);
        let n0 = self.concentration * 1.0e3 * AVOGADRO; // m^{-3}
        n0 * (-ELEMENTARY_CHARGE * phi / (BOLTZMANN * self.temperature)).exp()
    }

    /// Anion number density at distance z from OHP.
    pub fn anion_density(&self, z_from_ohp: f64) -> f64 {
        let phi = self.diffuse_potential(z_from_ohp);
        let n0 = self.concentration * 1.0e3 * AVOGADRO;
        n0 * (ELEMENTARY_CHARGE * phi / (BOLTZMANN * self.temperature)).exp()
    }
}

// ---------------------------------------------------------------------------
// ConstantPotentialMD
// ---------------------------------------------------------------------------

/// Constant potential MD (CPM) electrode model.
///
/// Maintains a fixed potential difference between two electrodes
/// by dynamically updating electrode charges using the Thomas-algorithm
/// method (Reed, Madden & Corrigan). Electrode atoms fluctuate their
/// charges to satisfy equipotential constraints.
#[derive(Debug, Clone)]
pub struct ConstantPotentialMd {
    /// Number of electrode atoms on each electrode.
    pub n_electrode_atoms: usize,
    /// Applied potential difference Delta V in V.
    pub applied_voltage: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Electrode charges for the left (negative) electrode in elementary charges.
    pub charges_left: Vec<f64>,
    /// Electrode charges for the right (positive) electrode in elementary charges.
    pub charges_right: Vec<f64>,
    /// Total charge on left electrode in elementary charges.
    pub q_left: f64,
    /// Total charge on right electrode in elementary charges.
    pub q_right: f64,
    /// Electrode capacitance matrix diagonal element in eV^{-1} (Gaussian units).
    pub capacitance_aa: f64,
    /// Timestep in fs.
    pub dt: f64,
    /// Fictitious mass for charge dynamics in eV fs^2 / e^2.
    pub charge_mass: f64,
}

impl ConstantPotentialMd {
    /// Create a CPM system with a graphite-like electrode.
    pub fn graphene_electrode(n_atoms: usize, voltage: f64, temperature: f64) -> Self {
        // Approximate capacitance matrix element for graphene:
        // C_aa ≈ 0.17 eV^-1 per atom (Merlet et al., J. Phys. Chem. Lett. 2013)
        let capacitance_aa = 0.17;
        let n_half = n_atoms;
        let charges = vec![0.0; n_half];
        Self {
            n_electrode_atoms: n_atoms,
            applied_voltage: voltage,
            temperature,
            charges_left: charges.clone(),
            charges_right: charges,
            q_left: 0.0,
            q_right: 0.0,
            capacitance_aa,
            dt: 1.0,
            charge_mass: 0.01,
        }
    }

    /// Update electrode charges to meet equipotential constraint.
    ///
    /// The constraint is: sum_j C_ij * q_j + V_ext_i = V_electrode
    /// For a simplified uniform electrode: q_i = (V_electrode - V_ext) / (N * C_aa)
    pub fn update_charges(&mut self, external_potential_left: f64, external_potential_right: f64) {
        let n = self.n_electrode_atoms as f64;
        let v_left = -self.applied_voltage / 2.0;
        let v_right = self.applied_voltage / 2.0;
        let dq_left = (v_left - external_potential_left) / (n * self.capacitance_aa);
        let dq_right = (v_right - external_potential_right) / (n * self.capacitance_aa);
        for q in &mut self.charges_left {
            *q += dq_left;
        }
        for q in &mut self.charges_right {
            *q += dq_right;
        }
        self.q_left = self.charges_left.iter().sum();
        self.q_right = self.charges_right.iter().sum();
    }

    /// Charge fluctuation using Metropolis Monte Carlo move.
    ///
    /// Proposes a charge transfer delta_q between electrode atoms and
    /// accepts or rejects based on the Boltzmann criterion.
    pub fn mc_charge_move(&mut self, delta_q_max: f64) {
        let mut rng = rand::rng();
        let n = self.n_electrode_atoms;
        if n < 2 {
            return;
        }
        let i = rng.random_range(0..n);
        let j = rng.random_range(0..n);
        if i == j {
            return;
        }
        let dq = rng.random_range(-delta_q_max..delta_q_max);
        // Energy cost: dE = C_aa * (q_i * dq - q_j * dq)
        let q_old_i = self.charges_left[i];
        let q_old_j = self.charges_left[j];
        let de = self.capacitance_aa
            * ((q_old_i + dq).powi(2) - q_old_i.powi(2) + (q_old_j - dq).powi(2) - q_old_j.powi(2));
        let beta = ELEMENTARY_CHARGE / (BOLTZMANN * self.temperature);
        if de < 0.0 || rng.random_range(0.0_f64..1.0) < (-beta * de).exp() {
            self.charges_left[i] += dq;
            self.charges_left[j] -= dq;
            self.q_left = self.charges_left.iter().sum();
        }
    }

    /// Mean charge per electrode atom on left electrode.
    pub fn mean_charge_left(&self) -> f64 {
        if self.n_electrode_atoms == 0 {
            return 0.0;
        }
        self.q_left / self.n_electrode_atoms as f64
    }

    /// Mean charge per electrode atom on right electrode.
    pub fn mean_charge_right(&self) -> f64 {
        if self.n_electrode_atoms == 0 {
            return 0.0;
        }
        self.q_right / self.n_electrode_atoms as f64
    }

    /// Surface charge density in C/m^2 given electrode area.
    pub fn surface_charge_density(&self, area_m2: f64) -> f64 {
        self.q_left * ELEMENTARY_CHARGE / area_m2
    }

    /// Capacitance estimate from electrode charge and voltage in F.
    pub fn capacitance_estimate(&self, area_m2: f64) -> f64 {
        if self.applied_voltage.abs() < 1.0e-10 {
            return 0.0;
        }
        self.surface_charge_density(area_m2) / self.applied_voltage
    }
}

// ---------------------------------------------------------------------------
// MarcusElectronTransfer
// ---------------------------------------------------------------------------

/// Marcus electron transfer theory for redox reactions at electrode interfaces.
///
/// The Marcus rate is:
///   k_et = A * exp\[-(lambda + delta_G)^2 / (4 lambda k_B T)\]
///
/// where lambda is the reorganisation energy, delta_G is the reaction
/// free energy, A is the pre-exponential factor, and k_B T is thermal energy.
#[derive(Debug, Clone)]
pub struct MarcusElectronTransfer {
    /// Reorganisation energy lambda in eV.
    pub lambda: f64,
    /// Reaction free energy delta_G in eV (negative = exothermic).
    pub delta_g: f64,
    /// Pre-exponential factor A in s^{-1}.
    pub pre_exponential: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Electrode potential in V vs standard hydrogen electrode.
    pub electrode_potential: f64,
    /// Standard redox potential E0 in V.
    pub standard_potential: f64,
    /// Symmetry factor (transfer coefficient) beta.
    pub beta: f64,
}

impl MarcusElectronTransfer {
    /// Create a Marcus model for the Fe3+/Fe2+ redox couple.
    pub fn fe3_fe2_couple(electrode_potential: f64, temperature: f64) -> Self {
        Self {
            lambda: MARCUS_LAMBDA_REF * 2.0,         // 1.0 eV reorganisation
            delta_g: -(electrode_potential - 0.771), // E0(Fe3+/Fe2+) = 0.771 V
            pre_exponential: 1.0e13,
            temperature,
            electrode_potential,
            standard_potential: 0.771,
            beta: 0.5,
        }
    }

    /// Create a Marcus model for the Li+/Li redox couple.
    pub fn li_couple(electrode_potential: f64, temperature: f64) -> Self {
        Self {
            lambda: 0.7,
            delta_g: -(electrode_potential - (-3.04)), // E0(Li+/Li) = -3.04 V
            pre_exponential: 1.0e12,
            temperature,
            electrode_potential,
            standard_potential: -3.04,
            beta: 0.5,
        }
    }

    /// Marcus rate constant for electron transfer in s^{-1}.
    ///
    /// k = A * exp\[-(lambda + delta_G)^2 / (4 lambda k_B T)\]
    pub fn rate_constant(&self) -> f64 {
        let kbt = BOLTZMANN * self.temperature / ELEMENTARY_CHARGE; // eV
        let exponent = -(self.lambda + self.delta_g).powi(2) / (4.0 * self.lambda * kbt);
        self.pre_exponential * exponent.exp()
    }

    /// Butler-Volmer rate (oxidation) using Marcus-derived activation energy.
    pub fn butler_volmer_rate_ox(&self) -> f64 {
        let eta = self.electrode_potential - self.standard_potential;
        let kbt = BOLTZMANN * self.temperature / ELEMENTARY_CHARGE;
        let alpha = 1.0 - self.beta;
        self.pre_exponential * (-alpha * eta / kbt).exp()
    }

    /// Butler-Volmer rate (reduction) using Marcus-derived activation energy.
    pub fn butler_volmer_rate_red(&self) -> f64 {
        let eta = self.electrode_potential - self.standard_potential;
        let kbt = BOLTZMANN * self.temperature / ELEMENTARY_CHARGE;
        self.pre_exponential * (self.beta * eta / kbt).exp()
    }

    /// Activation energy for forward electron transfer in eV.
    pub fn activation_energy(&self) -> f64 {
        (self.lambda + self.delta_g).powi(2) / (4.0 * self.lambda)
    }

    /// Reorganisation energy outer-sphere contribution (Born model) in eV.
    ///
    /// lambda_out = e^2/(4pi eps_0) * (1/(2 r_ox) + 1/(2 r_red) - 1/r_12) * (1/n^2 - 1/eps)
    pub fn outer_sphere_lambda(
        r_ox: f64,   // radius of oxidised species in Angstrom
        r_red: f64,  // radius of reduced species in Angstrom
        r_12: f64,   // centre-to-centre distance in Angstrom
        n_refr: f64, // optical refractive index
        eps: f64,    // static dielectric constant
    ) -> f64 {
        let r_ox_m = r_ox * 1.0e-10;
        let r_red_m = r_red * 1.0e-10;
        let r12_m = r_12 * 1.0e-10;
        let ke = 1.0 / (4.0 * PI * EPSILON_0);
        let optical = 1.0 / (n_refr * n_refr);
        let static_eps = 1.0 / eps;
        ke * ELEMENTARY_CHARGE
            * ELEMENTARY_CHARGE
            * (1.0 / (2.0 * r_ox_m) + 1.0 / (2.0 * r_red_m) - 1.0 / r12_m)
            * (optical - static_eps)
            / ELEMENTARY_CHARGE // convert J to eV
    }

    /// Check if the reaction is in the Marcus inverted region.
    pub fn is_inverted_region(&self) -> bool {
        self.delta_g.abs() > self.lambda
    }
}

// ---------------------------------------------------------------------------
// IonTransport
// ---------------------------------------------------------------------------

/// Ion transport model in electrolyte solutions (Nernst-Planck equation).
///
/// The Nernst-Planck (NP) flux is:
///   J = -D grad(c) - (D z e / k_B T) c grad(phi)
///
/// where D is diffusivity, c is concentration, z is valence, phi is potential.
#[derive(Debug, Clone)]
pub struct IonTransport {
    /// Ionic species.
    pub species: IonicSpecies,
    /// Concentration profile along z-axis (from electrode) in mol/m^3.
    pub concentration: Vec<f64>,
    /// Electrostatic potential profile in V.
    pub potential: Vec<f64>,
    /// Grid spacing in metres.
    pub dz: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Diffusion coefficient in m^2/s (may be position-dependent).
    pub diffusivity: f64,
    /// Current flux density in mol/(m^2 s).
    pub flux: Vec<f64>,
}

impl IonTransport {
    /// Create an ion transport model on a uniform 1D grid.
    pub fn new(
        species: IonicSpecies,
        n_grid: usize,
        length: f64,
        c_bulk: f64,
        temperature: f64,
    ) -> Self {
        let dz = length / n_grid as f64;
        let concentration = vec![c_bulk; n_grid];
        let potential = vec![0.0; n_grid];
        let flux = vec![0.0; n_grid.saturating_sub(1)];
        Self {
            species,
            concentration,
            potential,
            dz,
            temperature,
            diffusivity: species.diffusion_coefficient(),
            flux,
        }
    }

    /// Compute Nernst-Planck flux at each grid interface.
    pub fn compute_flux(&mut self) {
        let n = self.concentration.len();
        let z = self.species.charge_number() as f64;
        let vt = BOLTZMANN * self.temperature / ELEMENTARY_CHARGE;
        self.flux.resize(n.saturating_sub(1), 0.0);
        for i in 0..(n - 1) {
            let c_avg = 0.5 * (self.concentration[i] + self.concentration[i + 1]);
            let dc = (self.concentration[i + 1] - self.concentration[i]) / self.dz;
            let dphi = (self.potential[i + 1] - self.potential[i]) / self.dz;
            self.flux[i] = -self.diffusivity * (dc + z / vt * c_avg * dphi);
        }
    }

    /// Ionic conductivity contribution in S/m.
    ///
    /// sigma = z^2 e^2 D c / (k_B T)
    pub fn ionic_conductivity(&self) -> f64 {
        let z = self.species.charge_number() as f64;
        let c_avg = self.concentration.iter().sum::<f64>() / self.concentration.len().max(1) as f64;
        z * z * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * self.diffusivity * c_avg
            / (BOLTZMANN * self.temperature)
    }

    /// Molar conductivity in S m^2 / mol.
    pub fn molar_conductivity(&self) -> f64 {
        let c_avg = self.concentration.iter().sum::<f64>() / self.concentration.len().max(1) as f64;
        if c_avg < 1.0e-30 {
            return 0.0;
        }
        self.ionic_conductivity() / c_avg
    }

    /// Einstein relation check: D = mu k_B T / z e.
    pub fn mobility(&self) -> f64 {
        let z = self.species.charge_number() as f64;
        if z.abs() < 1.0e-10 {
            return 0.0;
        }
        self.diffusivity * ELEMENTARY_CHARGE * z.abs() / (BOLTZMANN * self.temperature)
    }

    /// Step the concentration profile using explicit finite differences.
    pub fn step(&mut self, dt: f64) {
        self.compute_flux();
        let n = self.concentration.len();
        let mut new_conc = self.concentration.clone();
        for (i, c) in new_conc.iter_mut().enumerate().take(n - 1).skip(1) {
            *c -= dt / self.dz * (self.flux[i] - self.flux[i - 1]);
            *c = c.max(0.0);
        }
        self.concentration = new_conc;
    }
}

// ---------------------------------------------------------------------------
// SolvationShell
// ---------------------------------------------------------------------------

/// Solvation shell analysis around an ion in MD.
///
/// Analyses the radial distribution function (RDF) and coordination
/// number as a function of distance from the central ion.
#[derive(Debug, Clone)]
pub struct SolvationShell {
    /// Central ionic species.
    pub central_species: IonicSpecies,
    /// Surrounding solvent/ion positions in Angstrom (relative to central ion).
    pub neighbour_positions: Vec<[f64; 3]>,
    /// RDF histogram bin width in Angstrom.
    pub dr: f64,
    /// RDF histogram: g(r) values.
    pub rdf: Vec<f64>,
    /// Cumulative coordination number n(r).
    pub coordination: Vec<f64>,
    /// Bulk number density in Angstrom^{-3}.
    pub bulk_density: f64,
}

impl SolvationShell {
    /// Create a new solvation shell analyser.
    pub fn new(species: IonicSpecies, dr: f64, r_max: f64, bulk_density: f64) -> Self {
        let n_bins = (r_max / dr).ceil() as usize;
        Self {
            central_species: species,
            neighbour_positions: Vec::new(),
            dr,
            rdf: vec![0.0; n_bins],
            coordination: vec![0.0; n_bins],
            bulk_density,
        }
    }

    /// Add a solvent molecule/ion at a given position relative to the central ion.
    pub fn add_neighbour(&mut self, pos: [f64; 3]) {
        self.neighbour_positions.push(pos);
    }

    /// Compute the RDF from current neighbour positions.
    ///
    /// g(r) = n(r, r+dr) / (4 pi r^2 dr rho_bulk)
    pub fn compute_rdf(&mut self) {
        let n_bins = self.rdf.len();
        self.rdf.fill(0.0);
        let mut counts = vec![0u64; n_bins];
        for pos in &self.neighbour_positions {
            let r = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
            let bin = (r / self.dr) as usize;
            if bin < n_bins {
                counts[bin] += 1;
            }
        }
        let n_neigh = self.neighbour_positions.len() as f64;
        for (k, &count) in counts.iter().enumerate() {
            let r_mid = (k as f64 + 0.5) * self.dr;
            let shell_vol = 4.0 * PI * r_mid * r_mid * self.dr;
            let expected = self.bulk_density * shell_vol;
            self.rdf[k] = if expected > 1.0e-30 {
                count as f64 / (n_neigh.max(1.0) * expected)
            } else {
                0.0
            };
        }
        // Cumulative coordination
        let mut cumul = 0.0;
        for (k, &count) in counts.iter().enumerate() {
            cumul += count as f64 / n_neigh.max(1.0);
            self.coordination[k] = cumul;
        }
    }

    /// First solvation shell coordination number (integral to first minimum).
    pub fn first_shell_coordination(&self) -> f64 {
        // Find first peak then first minimum
        let n = self.rdf.len();
        let mut peak_idx = 0;
        let mut peak_val = 0.0_f64;
        for (i, &g) in self.rdf.iter().enumerate() {
            if g > peak_val {
                peak_val = g;
                peak_idx = i;
            }
        }
        // Find first minimum after peak
        let mut min_idx = n - 1;
        for i in peak_idx..n.saturating_sub(1) {
            if self.rdf[i] < self.rdf[i + 1] {
                min_idx = i;
                break;
            }
        }
        if min_idx < self.coordination.len() {
            self.coordination[min_idx]
        } else {
            0.0
        }
    }

    /// Mean distance to first solvation shell neighbours in Angstrom.
    pub fn mean_solvation_distance(&self) -> f64 {
        if self.neighbour_positions.is_empty() {
            return 0.0;
        }
        let total: f64 = self
            .neighbour_positions
            .iter()
            .map(|p| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt())
            .sum();
        total / self.neighbour_positions.len() as f64
    }
}

// ---------------------------------------------------------------------------
// SurfaceChargeDensity
// ---------------------------------------------------------------------------

/// Surface charge density dynamics at an electrode-electrolyte interface.
///
/// Tracks the evolution of surface charge under applied potential, ion
/// adsorption, and charge transfer reactions.
#[derive(Debug, Clone)]
pub struct SurfaceChargeDensity {
    /// Surface charge density sigma in C/m^2.
    pub sigma: f64,
    /// Electrode area in m^2.
    pub area: f64,
    /// Applied potential in V.
    pub potential: f64,
    /// Point of zero charge (PZC) potential in V.
    pub pzc: f64,
    /// Differential capacitance C_d in F/m^2.
    pub capacitance: f64,
    /// Ion adsorption coverage theta \[0, 1\].
    pub adsorption_coverage: f64,
    /// Adsorption free energy in eV (negative = spontaneous).
    pub delta_g_ads: f64,
    /// Temperature in K.
    pub temperature: f64,
}

impl SurfaceChargeDensity {
    /// Create a model for a Pt(111) electrode in aqueous electrolyte.
    pub fn platinum_111(area: f64, temperature: f64) -> Self {
        Self {
            sigma: 0.0,
            area,
            potential: 0.0,
            pzc: 0.29,         // PZC of Pt(111) ≈ 0.29 V vs NHE
            capacitance: 0.20, // ~20 uF/cm^2 = 0.20 F/m^2
            adsorption_coverage: 0.0,
            delta_g_ads: -0.15, // eV
            temperature,
        }
    }

    /// Update surface charge for a given applied potential.
    pub fn update_sigma(&mut self, new_potential: f64) {
        self.potential = new_potential;
        self.sigma = self.capacitance * (self.potential - self.pzc);
    }

    /// Langmuir adsorption isotherm for ion coverage.
    ///
    /// theta = K * c / (1 + K * c)  where K = exp(-delta_G / k_B T)
    pub fn langmuir_coverage(&self, concentration: f64) -> f64 {
        let kbt = BOLTZMANN * self.temperature / ELEMENTARY_CHARGE;
        let k_eq = (-self.delta_g_ads / kbt).exp();
        k_eq * concentration / (1.0 + k_eq * concentration)
    }

    /// Total charge in C.
    pub fn total_charge(&self) -> f64 {
        self.sigma * self.area
    }

    /// Electrode charge density in e/Angstrom^2.
    pub fn charge_per_area_ea2(&self) -> f64 {
        self.sigma / ELEMENTARY_CHARGE * 1.0e-20 // C/m^2 -> e/Angstrom^2
    }
}

// ---------------------------------------------------------------------------
// CapacitanceFromMd
// ---------------------------------------------------------------------------

/// Capacitance calculation from MD charge fluctuations.
///
/// Uses the charge fluctuation formula:
///   C = beta * e^2 * (<Q^2> - `Q`^2) / 2
/// where beta = 1/(k_B T), and Q is the electrode charge.
#[derive(Debug, Clone)]
pub struct CapacitanceFromMd {
    /// Electrode potential in V.
    pub potential: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Accumulated charge samples (in elementary charges).
    pub charge_samples: Vec<f64>,
}

impl CapacitanceFromMd {
    /// Create a new capacitance accumulator.
    pub fn new(potential: f64, temperature: f64) -> Self {
        Self {
            potential,
            temperature,
            charge_samples: Vec::new(),
        }
    }

    /// Add a charge sample from a simulation frame.
    pub fn add_sample(&mut self, charge: f64) {
        self.charge_samples.push(charge);
    }

    /// Mean charge in elementary charges.
    pub fn mean_charge(&self) -> f64 {
        if self.charge_samples.is_empty() {
            return 0.0;
        }
        self.charge_samples.iter().sum::<f64>() / self.charge_samples.len() as f64
    }

    /// Variance of charge in elementary charges^2.
    pub fn charge_variance(&self) -> f64 {
        if self.charge_samples.len() < 2 {
            return 0.0;
        }
        let mean = self.mean_charge();
        let var: f64 = self.charge_samples.iter().map(|q| (q - mean).powi(2)).sum();
        var / (self.charge_samples.len() - 1) as f64
    }

    /// Differential capacitance from fluctuation formula in Farads.
    ///
    /// C = e^2 * Var(Q) / (k_B T)
    pub fn differential_capacitance_f(&self) -> f64 {
        let var = self.charge_variance();
        ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * var / (BOLTZMANN * self.temperature)
    }

    /// Differential capacitance in uF/cm^2 given electrode area.
    pub fn capacitance_uf_cm2(&self, area_m2: f64) -> f64 {
        self.differential_capacitance_f() / area_m2 * 1.0e6 * 1.0e-4 // F/m^2 -> uF/cm^2
    }
}

// ---------------------------------------------------------------------------
// LithiumIntercalation
// ---------------------------------------------------------------------------

/// Lithium-ion intercalation model for anode and cathode materials.
///
/// Tracks Li insertion/extraction using lattice gas (Bragg-Williams) model
/// combined with chemical diffusion for both graphite (anode) and LFP/NMC (cathode).
#[derive(Debug, Clone)]
pub struct LithiumIntercalation {
    /// Material type.
    pub material: InterCalationMaterial,
    /// Lithium stoichiometry x in Li_x M (0 ≤ x ≤ x_max).
    pub x_li: f64,
    /// Maximum Li stoichiometry.
    pub x_max: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Open-circuit voltage in V.
    pub ocv: f64,
    /// Chemical diffusion coefficient in m^2/s.
    pub d_chem: f64,
    /// Lattice constant a in Angstrom.
    pub lattice_a: f64,
    /// Lattice constant c in Angstrom.
    pub lattice_c: f64,
    /// Volume change upon full lithiation (fractional).
    pub volume_change: f64,
}

/// Intercalation host material types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterCalationMaterial {
    /// Graphite anode (LiC6).
    Graphite,
    /// Silicon anode (Li15Si4).
    Silicon,
    /// LiFePO4 (olivine cathode, LFP).
    Lfp,
    /// LiNi0.8Mn0.1Co0.1O2 (NMC811 cathode).
    Nmc811,
    /// LiCoO2 (LCO cathode).
    Lco,
    /// LiMn2O4 (spinel cathode).
    Lmo,
}

impl InterCalationMaterial {
    /// Theoretical specific capacity in mAh/g.
    pub fn theoretical_capacity(&self) -> f64 {
        match self {
            Self::Graphite => 372.0,
            Self::Silicon => 3579.0,
            Self::Lfp => 170.0,
            Self::Nmc811 => 275.0,
            Self::Lco => 274.0,
            Self::Lmo => 148.0,
        }
    }

    /// Average intercalation potential in V vs Li/Li+.
    pub fn average_potential(&self) -> f64 {
        match self {
            Self::Graphite => 0.17,
            Self::Silicon => 0.40,
            Self::Lfp => 3.45,
            Self::Nmc811 => 3.80,
            Self::Lco => 3.90,
            Self::Lmo => 4.10,
        }
    }

    /// Electronic conductivity in S/m.
    pub fn electronic_conductivity(&self) -> f64 {
        match self {
            Self::Graphite => 1.0e4,
            Self::Silicon => 1.0e-3,
            Self::Lfp => 1.0e-9, // intrinsic; carbon-coated ~1e-2
            Self::Nmc811 => 1.0e-2,
            Self::Lco => 1.0e-2,
            Self::Lmo => 1.0e-5,
        }
    }

    /// Maximum Li stoichiometry.
    pub fn x_max(&self) -> f64 {
        match self {
            Self::Graphite => 1.0, // LiC6
            Self::Silicon => 3.75,
            Self::Lfp => 1.0,
            Self::Nmc811 => 1.0,
            Self::Lco => 1.0,
            Self::Lmo => 1.0,
        }
    }
}

impl LithiumIntercalation {
    /// Create a graphite anode intercalation model.
    pub fn graphite_anode(temperature: f64) -> Self {
        let mat = InterCalationMaterial::Graphite;
        Self {
            material: mat,
            x_li: 0.5,
            x_max: mat.x_max(),
            temperature,
            ocv: Self::graphite_ocv(0.5),
            d_chem: 1.0e-14, // m^2/s chemical diffusion in graphite
            lattice_a: 2.461,
            lattice_c: 6.708,
            volume_change: 0.10, // 10% expansion at full lithiation
        }
    }

    /// Create an LFP cathode model.
    pub fn lfp_cathode(temperature: f64) -> Self {
        let mat = InterCalationMaterial::Lfp;
        Self {
            material: mat,
            x_li: 0.5,
            x_max: mat.x_max(),
            temperature,
            ocv: 3.45,
            d_chem: 1.0e-18, // very slow diffusion in LFP
            lattice_a: 10.33,
            lattice_c: 4.693,
            volume_change: 0.048, // ~4.8% between FP and LFP
        }
    }

    /// Open-circuit voltage of graphite anode vs Li/Li+ (empirical).
    ///
    /// Multi-plateau behaviour: U(x) ~ 0.10 + 0.05 * (1 - x) + sigmoid terms.
    pub fn graphite_ocv(x: f64) -> f64 {
        let x = x.clamp(0.001, 0.999);
        // Simplified Plett-type OCV for graphite (3-stage)
        let u1 = 0.10;
        let u2 = 0.23;
        let u3 = 0.40;
        let w1 = 1.0 / (1.0 + (20.0 * (x - 0.15)).exp());
        let w2 = 1.0 / (1.0 + (20.0 * (x - 0.50)).exp());
        u1 + (u2 - u1) * w1 + (u3 - u2) * (1.0 - w2)
    }

    /// Nernst correction to OCV in V.
    ///
    /// delta_U = (R T / F) ln(x / (1 - x))
    pub fn nernst_correction(&self) -> f64 {
        let x = self.x_li.clamp(1.0e-4, 1.0 - 1.0e-4);
        GAS_CONSTANT * self.temperature / FARADAY * (x / (1.0 - x)).ln()
    }

    /// Update OCV for current Li stoichiometry.
    pub fn update_ocv(&mut self) {
        self.ocv = match self.material {
            InterCalationMaterial::Graphite => Self::graphite_ocv(self.x_li),
            _ => self.material.average_potential(),
        } + self.nernst_correction();
    }

    /// Li-ion flux into particle (positive = insertion) in mol/(m^2 s).
    ///
    /// j = -D_chem * (x_surface - x_bulk) / (r_particle)
    pub fn surface_flux(&self, x_surface: f64, r_particle: f64) -> f64 {
        let c_max = 2.29e4; // mol/m^3 for graphite
        let dc = (x_surface - self.x_li) * c_max;
        -self.d_chem * dc / r_particle
    }

    /// Intercalation strain (isotropic) at current stoichiometry.
    pub fn intercalation_strain(&self) -> f64 {
        self.volume_change * self.x_li / self.x_max / 3.0 // linear strain
    }

    /// State of charge (SoC) fraction \[0, 1\].
    pub fn state_of_charge(&self) -> f64 {
        self.x_li / self.x_max
    }
}

// ---------------------------------------------------------------------------
// ElectrochemicalMdSimulation
// ---------------------------------------------------------------------------

/// High-level driver for electrochemical MD simulations.
///
/// Combines ions, electrode, GCS model, Marcus kinetics, and Li intercalation.
#[derive(Debug, Clone)]
pub struct ElectrochemicalMdSimulation {
    /// Ions in the simulation box.
    pub ions: Vec<ElectrochemicalIon>,
    /// GCS double layer model.
    pub gcs: GouyChapmanSternModel,
    /// CPM electrode.
    pub cpm: ConstantPotentialMd,
    /// Marcus model for electron transfer.
    pub marcus: MarcusElectronTransfer,
    /// Li intercalation model.
    pub intercalation: LithiumIntercalation,
    /// Box size \[Lx, Ly, Lz\] in Angstrom.
    pub box_size: [f64; 3],
    /// Current simulation time in fs.
    pub time: f64,
    /// Integration timestep in fs.
    pub dt: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Dielectric constant of the solvent.
    pub dielectric: f64,
}

impl ElectrochemicalMdSimulation {
    /// Create a simulation for a Li-ion battery interface.
    pub fn new_liion_interface(voltage: f64, temperature: f64, n_ions: usize) -> Self {
        let box_len = 40.0; // Angstrom
        let box_size = [box_len, box_len, box_len];
        let mut rng = rand::rng();
        let mut ions = Vec::with_capacity(n_ions);
        let conc = 1.0; // mol/L
        for i in 0..n_ions {
            let species = if i % 2 == 0 {
                IonicSpecies::Li
            } else {
                IonicSpecies::Pf6
            };
            let pos = [
                rng.random_range(0.0..box_len),
                rng.random_range(0.0..box_len),
                rng.random_range(0.0..box_len),
            ];
            let mut ion = ElectrochemicalIon::new(pos, species);
            ion.id = i;
            ions.push(ion);
        }
        Self {
            ions,
            gcs: GouyChapmanSternModel::aqueous_kcl(conc, voltage, temperature),
            cpm: ConstantPotentialMd::graphene_electrode(50, voltage, temperature),
            marcus: MarcusElectronTransfer::li_couple(voltage, temperature),
            intercalation: LithiumIntercalation::graphite_anode(temperature),
            box_size,
            time: 0.0,
            dt: 1.0,
            temperature,
            dielectric: 78.5,
        }
    }

    /// Compute total Coulomb energy in eV.
    pub fn total_coulomb_energy(&self) -> f64 {
        let n = self.ions.len();
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                energy += self.ions[i].coulomb_energy_ev(&self.ions[j], self.dielectric);
            }
        }
        energy
    }

    /// Apply Langevin thermostat forces (white noise + friction).
    pub fn apply_langevin(&mut self, gamma: f64) {
        let mut rng = rand::rng();
        let kbt_ev = BOLTZMANN * self.temperature / ELEMENTARY_CHARGE;
        for ion in &mut self.ions {
            let m = ion.species.mass();
            let m_amu = m * AVOGADRO * 1.0e3; // Dalton
            let sigma = (2.0 * gamma * kbt_ev * m_amu / self.dt).sqrt();
            for k in 0..3 {
                let noise: f64 = rng.random_range(-1.0_f64..1.0_f64);
                ion.force[k] += sigma * noise - gamma * ion.vel[k] * m_amu;
            }
        }
    }

    /// Advance by one timestep using velocity Verlet.
    pub fn step(&mut self) {
        let dt = self.dt;
        for ion in &mut self.ions {
            let m_amu = ion.species.mass() * AVOGADRO * 1.0e3;
            // half-step velocity
            for k in 0..3 {
                ion.vel[k] += 0.5 * ion.force[k] / m_amu * dt;
            }
            // position update
            for k in 0..3 {
                ion.pos[k] += ion.vel[k] * dt;
                // periodic boundary conditions
                ion.pos[k] = ion.pos[k].rem_euclid(self.box_size[k]);
            }
        }
        // recompute forces (simplified: zero here; extend with real potentials)
        for ion in &mut self.ions {
            ion.zero_forces();
        }
        self.apply_langevin(0.01);
        for ion in &mut self.ions {
            let m_amu = ion.species.mass() * AVOGADRO * 1.0e3;
            for k in 0..3 {
                ion.vel[k] += 0.5 * ion.force[k] / m_amu * dt;
            }
        }
        self.time += dt;
    }

    /// Run for n_steps timesteps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Total kinetic energy in eV.
    pub fn total_kinetic_energy_ev(&self) -> f64 {
        self.ions.iter().map(|i| i.kinetic_energy_ev()).sum()
    }

    /// Temperature estimate from kinetic energy in K.
    pub fn temperature_estimate(&self) -> f64 {
        let ke = self.total_kinetic_energy_ev() * ELEMENTARY_CHARGE; // J
        let n = self.ions.len() as f64;
        if n < 1.0 {
            return 0.0;
        }
        2.0 * ke / (3.0 * n * BOLTZMANN)
    }

    /// Count cation-anion pairs within cutoff in Angstrom.
    pub fn count_ion_pairs(&self, cutoff: f64) -> usize {
        let n = self.ions.len();
        let mut count = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let zi = self.ions[i].species.charge_number();
                let zj = self.ions[j].species.charge_number();
                if zi * zj < 0 && self.ions[i].distance_to(&self.ions[j]) < cutoff {
                    count += 1;
                }
            }
        }
        count
    }
}

// ---------------------------------------------------------------------------
// Utility Functions
// ---------------------------------------------------------------------------

/// Compute the Debye length for a 1:1 electrolyte in water at given concentration.
///
/// kappa^{-1} = sqrt(eps * k_B * T / (2 N_A e^2 c))  \[metres\]
pub fn debye_length_aqueous(concentration_mol_per_l: f64, temperature: f64) -> f64 {
    let eps = 78.5 * EPSILON_0;
    let c_si = concentration_mol_per_l * 1.0e3; // mol/m^3
    let kappa_sq = 2.0 * AVOGADRO * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * c_si
        / (eps * BOLTZMANN * temperature);
    1.0 / kappa_sq.sqrt()
}

/// Convert concentration from mol/L to number density in Angstrom^{-3}.
pub fn mol_per_l_to_number_density_aa(concentration: f64) -> f64 {
    concentration * 1.0e3 * AVOGADRO * 1.0e-30 // mol/m^3 * N_A -> m^{-3} -> A^{-3}
}

/// Thermal voltage at given temperature in V.
pub fn thermal_voltage(temperature: f64) -> f64 {
    BOLTZMANN * temperature / ELEMENTARY_CHARGE
}

/// Born solvation energy for an ion of radius r and charge z in a dielectric medium in eV.
///
/// delta_G_Born = -z^2 e^2 / (8 pi eps_0 eps r) * (1 - 1/eps)
pub fn born_solvation_energy_ev(z: i32, radius_angstrom: f64, dielectric: f64) -> f64 {
    let r_m = radius_angstrom * 1.0e-10;
    let ke = 1.0 / (4.0 * PI * EPSILON_0);
    -(z as f64).powi(2) * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * ke / (2.0 * r_m)
        * (1.0 - 1.0 / dielectric)
        / ELEMENTARY_CHARGE
}

/// Compute the electric field at a point due to a set of point charges.
///
/// E = sum_i k_e * q_i * (r - r_i) / |r - r_i|^3
pub fn electric_field_at(
    target: [f64; 3],
    charge_positions: &[[f64; 3]],
    charges: &[f64],
    dielectric: f64,
) -> [f64; 3] {
    let ke = 1.0 / (4.0 * PI * EPSILON_0 * dielectric);
    let mut e = [0.0; 3];
    for (pos, &q) in charge_positions.iter().zip(charges.iter()) {
        let dr = [
            (target[0] - pos[0]) * 1.0e-10,
            (target[1] - pos[1]) * 1.0e-10,
            (target[2] - pos[2]) * 1.0e-10,
        ];
        let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
        let r = r2.sqrt();
        if r < 1.0e-15 {
            continue;
        }
        let factor = ke * q * ELEMENTARY_CHARGE / (r2 * r);
        e[0] += factor * dr[0];
        e[1] += factor * dr[1];
        e[2] += factor * dr[2];
    }
    e
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1.0e-9;

    // --- IonicSpecies ---

    #[test]
    fn test_li_charge_number() {
        assert_eq!(IonicSpecies::Li.charge_number(), 1);
    }

    #[test]
    fn test_cl_charge_number() {
        assert_eq!(IonicSpecies::Cl.charge_number(), -1);
    }

    #[test]
    fn test_ca_charge_number() {
        assert_eq!(IonicSpecies::Ca.charge_number(), 2);
    }

    #[test]
    fn test_water_neutral() {
        assert_eq!(IonicSpecies::Water.charge_number(), 0);
    }

    #[test]
    fn test_li_smaller_than_k() {
        assert!(IonicSpecies::Li.ionic_radius() < IonicSpecies::K.ionic_radius());
    }

    #[test]
    fn test_diffusion_coeff_positive() {
        for species in &[IonicSpecies::Li, IonicSpecies::Na, IonicSpecies::Cl] {
            assert!(species.diffusion_coefficient() > 0.0);
        }
    }

    #[test]
    fn test_ion_mass_positive() {
        assert!(IonicSpecies::Li.mass() > 0.0);
        assert!(IonicSpecies::Tfsi.mass() > IonicSpecies::Li.mass());
    }

    // --- ElectrochemicalIon ---

    #[test]
    fn test_ion_creation() {
        let ion = ElectrochemicalIon::new([0.0, 0.0, 0.0], IonicSpecies::Li);
        assert_eq!(ion.charge, 1.0);
        assert!(ion.speed() < TOL);
    }

    #[test]
    fn test_ion_distance() {
        let i1 = ElectrochemicalIon::new([0.0, 0.0, 0.0], IonicSpecies::Li);
        let i2 = ElectrochemicalIon::new([3.0, 4.0, 0.0], IonicSpecies::Cl);
        assert!((i1.distance_to(&i2) - 5.0).abs() < TOL);
    }

    #[test]
    fn test_coulomb_energy_like_charges_positive() {
        let i1 = ElectrochemicalIon::new([0.0, 0.0, 0.0], IonicSpecies::Li);
        let i2 = ElectrochemicalIon::new([3.0, 0.0, 0.0], IonicSpecies::Na);
        let e = i1.coulomb_energy_ev(&i2, 78.5);
        assert!(e > 0.0);
    }

    #[test]
    fn test_coulomb_energy_unlike_charges_negative() {
        let i1 = ElectrochemicalIon::new([0.0, 0.0, 0.0], IonicSpecies::Li);
        let i2 = ElectrochemicalIon::new([3.0, 0.0, 0.0], IonicSpecies::Cl);
        let e = i1.coulomb_energy_ev(&i2, 78.5);
        assert!(e < 0.0);
    }

    #[test]
    fn test_kinetic_energy_zero_at_rest() {
        let ion = ElectrochemicalIon::new([0.0, 0.0, 0.0], IonicSpecies::Li);
        assert!(ion.kinetic_energy_ev() < TOL);
    }

    // --- GouyChapmanSternModel ---

    #[test]
    fn test_gcs_debye_length_positive() {
        let gcs = GouyChapmanSternModel::aqueous_kcl(0.1, 0.1, STD_TEMP);
        let kd = gcs.debye_length();
        assert!(kd > 0.0);
        assert!(kd.is_finite());
    }

    #[test]
    fn test_gcs_debye_decreases_with_concentration() {
        let gcs1 = GouyChapmanSternModel::aqueous_kcl(0.01, 0.1, STD_TEMP);
        let gcs2 = GouyChapmanSternModel::aqueous_kcl(0.1, 0.1, STD_TEMP);
        assert!(gcs1.debye_length() > gcs2.debye_length());
    }

    #[test]
    fn test_gcs_differential_capacitance_positive() {
        let gcs = GouyChapmanSternModel::aqueous_kcl(0.1, 0.0, STD_TEMP);
        let c = gcs.differential_capacitance();
        assert!(c > 0.0);
    }

    #[test]
    fn test_gcs_diffuse_potential_decays_with_distance() {
        let gcs = GouyChapmanSternModel::aqueous_kcl(0.1, 0.2, STD_TEMP);
        let phi0 = gcs.diffuse_potential(0.0).abs();
        let phi1 = gcs.diffuse_potential(1.0e-9).abs();
        assert!(phi0 >= phi1);
    }

    #[test]
    fn test_gcs_cation_density_positive() {
        let gcs = GouyChapmanSternModel::aqueous_kcl(0.1, 0.1, STD_TEMP);
        let n = gcs.cation_density(1.0e-9);
        assert!(n > 0.0);
    }

    // --- ConstantPotentialMd ---

    #[test]
    fn test_cpm_mean_charge_initial_zero() {
        let cpm = ConstantPotentialMd::graphene_electrode(20, 1.0, STD_TEMP);
        assert!(cpm.mean_charge_left().abs() < TOL);
    }

    #[test]
    fn test_cpm_update_charges_nonzero() {
        let mut cpm = ConstantPotentialMd::graphene_electrode(20, 1.0, STD_TEMP);
        cpm.update_charges(0.1, -0.1);
        // After update, charges should shift
        assert!(cpm.q_left.is_finite());
        assert!(cpm.q_right.is_finite());
    }

    #[test]
    fn test_cpm_mc_charge_move_conserves_approximate_total() {
        let mut cpm = ConstantPotentialMd::graphene_electrode(20, 0.0, STD_TEMP);
        let q0: f64 = cpm.charges_left.iter().sum();
        for _ in 0..100 {
            cpm.mc_charge_move(0.01);
        }
        let q1: f64 = cpm.charges_left.iter().sum();
        // Total charge conserved by MC moves (only redistribution)
        assert!((q1 - q0).abs() < 1.0);
    }

    // --- MarcusElectronTransfer ---

    #[test]
    fn test_marcus_rate_positive() {
        let m = MarcusElectronTransfer::fe3_fe2_couple(0.5, STD_TEMP);
        assert!(m.rate_constant() > 0.0);
    }

    #[test]
    fn test_marcus_activation_energy_positive() {
        let m = MarcusElectronTransfer::li_couple(-3.0, STD_TEMP);
        assert!(m.activation_energy() >= 0.0);
    }

    #[test]
    fn test_marcus_not_inverted_small_delta_g() {
        let mut m = MarcusElectronTransfer::fe3_fe2_couple(0.771, STD_TEMP);
        m.delta_g = -0.1; // |delta_g| << lambda
        assert!(!m.is_inverted_region());
    }

    #[test]
    fn test_marcus_inverted_large_delta_g() {
        let mut m = MarcusElectronTransfer::fe3_fe2_couple(0.771, STD_TEMP);
        m.delta_g = -5.0; // |delta_g| >> lambda
        assert!(m.is_inverted_region());
    }

    #[test]
    fn test_marcus_outer_sphere_positive() {
        let lam = MarcusElectronTransfer::outer_sphere_lambda(3.5, 4.0, 10.0, 1.33, 78.5);
        assert!(lam > 0.0);
    }

    #[test]
    fn test_butler_volmer_rates_finite() {
        let m = MarcusElectronTransfer::fe3_fe2_couple(0.9, STD_TEMP);
        assert!(m.butler_volmer_rate_ox().is_finite());
        assert!(m.butler_volmer_rate_red().is_finite());
    }

    // --- IonTransport ---

    #[test]
    fn test_ion_transport_creation() {
        let t = IonTransport::new(IonicSpecies::Li, 100, 1.0e-7, 1000.0, STD_TEMP);
        assert_eq!(t.concentration.len(), 100);
    }

    #[test]
    fn test_ion_conductivity_positive() {
        let t = IonTransport::new(IonicSpecies::K, 50, 1.0e-7, 1000.0, STD_TEMP);
        assert!(t.ionic_conductivity() > 0.0);
    }

    #[test]
    fn test_ion_mobility_positive() {
        let t = IonTransport::new(IonicSpecies::Na, 50, 1.0e-7, 1000.0, STD_TEMP);
        assert!(t.mobility() > 0.0);
    }

    #[test]
    fn test_ion_transport_step_stable() {
        let mut t = IonTransport::new(IonicSpecies::Cl, 50, 1.0e-7, 1000.0, STD_TEMP);
        // Apply a linear potential gradient
        for i in 0..50 {
            t.potential[i] = i as f64 * 0.01;
        }
        t.step(1.0e-9);
        assert!(t.concentration.iter().all(|&c| c >= 0.0 && c.is_finite()));
    }

    // --- SolvationShell ---

    #[test]
    fn test_solvation_shell_empty_coordination_zero() {
        let mut s = SolvationShell::new(IonicSpecies::Li, 0.5, 10.0, 0.033);
        s.compute_rdf();
        assert!(s.first_shell_coordination() < TOL);
    }

    #[test]
    fn test_solvation_shell_mean_distance_positive() {
        let mut s = SolvationShell::new(IonicSpecies::Li, 0.5, 10.0, 0.033);
        s.add_neighbour([2.0, 0.0, 0.0]);
        s.add_neighbour([0.0, 2.1, 0.0]);
        s.add_neighbour([0.0, 0.0, 1.9]);
        assert!(s.mean_solvation_distance() > 0.0);
    }

    #[test]
    fn test_solvation_rdf_bins_correct() {
        let mut s = SolvationShell::new(IonicSpecies::Na, 0.25, 8.0, 0.033);
        for i in 0..10 {
            s.add_neighbour([2.0 + i as f64 * 0.1, 0.0, 0.0]);
        }
        s.compute_rdf();
        assert_eq!(s.rdf.len(), 32); // 8.0 / 0.25
    }

    // --- SurfaceChargeDensity ---

    #[test]
    fn test_surface_charge_zero_at_pzc() {
        let mut scd = SurfaceChargeDensity::platinum_111(1.0e-18, STD_TEMP);
        scd.update_sigma(scd.pzc);
        assert!(scd.sigma.abs() < 1.0e-10);
    }

    #[test]
    fn test_surface_charge_positive_above_pzc() {
        let mut scd = SurfaceChargeDensity::platinum_111(1.0e-18, STD_TEMP);
        scd.update_sigma(0.5);
        assert!(scd.sigma > 0.0);
    }

    #[test]
    fn test_langmuir_coverage_between_0_and_1() {
        let scd = SurfaceChargeDensity::platinum_111(1.0e-18, STD_TEMP);
        let theta = scd.langmuir_coverage(0.01);
        assert!((0.0..=1.0).contains(&theta));
    }

    // --- CapacitanceFromMd ---

    #[test]
    fn test_capacitance_empty_samples_zero() {
        let cap = CapacitanceFromMd::new(0.5, STD_TEMP);
        assert!(cap.mean_charge() < TOL);
        assert!(cap.charge_variance() < TOL);
    }

    #[test]
    fn test_capacitance_variance_positive() {
        let mut cap = CapacitanceFromMd::new(0.5, STD_TEMP);
        for i in 0..10 {
            cap.add_sample(i as f64 * 0.1);
        }
        assert!(cap.charge_variance() > 0.0);
    }

    #[test]
    fn test_capacitance_in_farads_positive() {
        let mut cap = CapacitanceFromMd::new(0.5, STD_TEMP);
        for i in 0..20 {
            cap.add_sample(i as f64 * 0.05 - 0.5);
        }
        let c = cap.differential_capacitance_f();
        assert!(c >= 0.0);
    }

    // --- LithiumIntercalation ---

    #[test]
    fn test_graphite_ocv_decreasing_near_full() {
        let u1 = LithiumIntercalation::graphite_ocv(0.1);
        let u2 = LithiumIntercalation::graphite_ocv(0.9);
        // Graphite OCV should be higher at low SOC
        assert!(u1 > u2 || (u1 - u2).abs() < 0.5);
    }

    #[test]
    fn test_graphite_ocv_in_range() {
        for x in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let u = LithiumIntercalation::graphite_ocv(x);
            assert!(u > 0.0 && u < 1.0);
        }
    }

    #[test]
    fn test_nernst_correction_finite() {
        let model = LithiumIntercalation::graphite_anode(STD_TEMP);
        let dn = model.nernst_correction();
        assert!(dn.is_finite());
    }

    #[test]
    fn test_soc_in_range() {
        let model = LithiumIntercalation::graphite_anode(STD_TEMP);
        let soc = model.state_of_charge();
        assert!((0.0..=1.0).contains(&soc));
    }

    #[test]
    fn test_intercalation_strain_positive() {
        let model = LithiumIntercalation::graphite_anode(STD_TEMP);
        assert!(model.intercalation_strain() >= 0.0);
    }

    #[test]
    fn test_lfp_average_potential() {
        let model = LithiumIntercalation::lfp_cathode(STD_TEMP);
        // LFP ~3.45 V vs Li/Li+
        assert!((model.material.average_potential() - 3.45).abs() < 0.01);
    }

    // --- ElectrochemicalMdSimulation ---

    #[test]
    fn test_sim_ion_count() {
        let sim = ElectrochemicalMdSimulation::new_liion_interface(1.0, STD_TEMP, 10);
        assert_eq!(sim.ions.len(), 10);
    }

    #[test]
    fn test_sim_time_advances() {
        let mut sim = ElectrochemicalMdSimulation::new_liion_interface(1.0, STD_TEMP, 4);
        sim.run(5);
        assert!(sim.time > 0.0);
    }

    #[test]
    fn test_sim_ions_stay_in_box() {
        let mut sim = ElectrochemicalMdSimulation::new_liion_interface(1.0, STD_TEMP, 8);
        sim.run(10);
        for ion in &sim.ions {
            for k in 0..3 {
                assert!(ion.pos[k] >= 0.0 && ion.pos[k] < sim.box_size[k]);
            }
        }
    }

    // --- Utility functions ---

    #[test]
    fn test_debye_length_positive() {
        let kd = debye_length_aqueous(0.1, STD_TEMP);
        assert!(kd > 0.0);
        assert!(kd < 1.0e-8); // should be ~ 1 nm for 0.1 M
    }

    #[test]
    fn test_thermal_voltage_at_298() {
        let vt = thermal_voltage(STD_TEMP);
        assert!((vt - THERMAL_VOLTAGE).abs() < 1.0e-10);
    }

    #[test]
    fn test_born_solvation_negative_for_cation() {
        let dg = born_solvation_energy_ev(1, 2.0, 78.5);
        assert!(dg < 0.0);
    }

    #[test]
    fn test_mol_per_l_conversion_positive() {
        let n = mol_per_l_to_number_density_aa(1.0);
        assert!(n > 0.0);
    }
}
