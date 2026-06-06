// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Battery molecular dynamics module for Li-ion battery simulation.
//!
//! This module provides MD-level models for simulating:
//! - Lattice atoms with species identification (Li, C, O, P, Fe, Mn)
//! - Electrolyte molecules with solvent properties
//! - Intercalation sites in host lattices
//! - Solid Electrolyte Interface (SEI) layer growth
//! - Butler-Volmer charge-transfer kinetics
//! - Full battery cell state tracking
//! - Cycle degradation and capacity fade
//! - Battery analysis (Ragone plot, coulombic efficiency, MSD diffusivity)

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Faraday constant in C/mol.
const FARADAY: f64 = 96_485.332_9;

/// Universal gas constant in J/(mol K).
const GAS_CONSTANT: f64 = 8.314;

/// Boltzmann constant in J/K.
const BOLTZMANN: f64 = 1.380649e-23;

/// Elementary charge in C.
const ELEMENTARY_CHARGE: f64 = 1.602176634e-19;

/// Avogadro's number.
const AVOGADRO: f64 = 6.02214076e23;

/// Vacuum permittivity in F/m.
const EPSILON_0: f64 = 8.854187817e-12;

/// Reference temperature in K (25 C).
const REF_TEMPERATURE: f64 = 298.15;

/// Li atomic mass in kg.
const LI_MASS: f64 = 6.941e-3 / AVOGADRO;

/// Carbon atomic mass in kg.
const C_MASS: f64 = 12.011e-3 / AVOGADRO;

/// Oxygen atomic mass in kg.
const O_MASS: f64 = 15.999e-3 / AVOGADRO;

/// Phosphorus atomic mass in kg.
const P_MASS: f64 = 30.974e-3 / AVOGADRO;

/// Iron atomic mass in kg.
const FE_MASS: f64 = 55.845e-3 / AVOGADRO;

/// Manganese atomic mass in kg.
const MN_MASS: f64 = 54.938e-3 / AVOGADRO;

// ---------------------------------------------------------------------------
// AtomSpecies
// ---------------------------------------------------------------------------

/// Atomic species in the battery system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomSpecies {
    /// Lithium.
    Li,
    /// Carbon (graphite anode).
    C,
    /// Oxygen.
    O,
    /// Phosphorus.
    P,
    /// Iron.
    Fe,
    /// Manganese.
    Mn,
}

impl AtomSpecies {
    /// Atomic mass in kg.
    pub fn mass(&self) -> f64 {
        match self {
            Self::Li => LI_MASS,
            Self::C => C_MASS,
            Self::O => O_MASS,
            Self::P => P_MASS,
            Self::Fe => FE_MASS,
            Self::Mn => MN_MASS,
        }
    }

    /// Default partial charge in elementary charge units.
    pub fn default_charge(&self) -> f64 {
        match self {
            Self::Li => 1.0,
            Self::C => 0.0,
            Self::O => -2.0,
            Self::P => 5.0,
            Self::Fe => 2.0,
            Self::Mn => 2.0,
        }
    }

    /// Lennard-Jones sigma in Angstrom (approximate).
    pub fn lj_sigma(&self) -> f64 {
        match self {
            Self::Li => 2.18,
            Self::C => 3.40,
            Self::O => 3.12,
            Self::P => 3.56,
            Self::Fe => 2.91,
            Self::Mn => 2.96,
        }
    }

    /// Lennard-Jones epsilon in eV (approximate).
    pub fn lj_epsilon(&self) -> f64 {
        match self {
            Self::Li => 0.002,
            Self::C => 0.003,
            Self::O => 0.008,
            Self::P => 0.010,
            Self::Fe => 0.013,
            Self::Mn => 0.012,
        }
    }
}

// ---------------------------------------------------------------------------
// LatticeAtom
// ---------------------------------------------------------------------------

/// An atom in the battery lattice.
#[derive(Debug, Clone)]
pub struct LatticeAtom {
    /// Position \[x, y, z\] in Angstrom.
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\] in Angstrom/fs.
    pub vel: [f64; 3],
    /// Force \[fx, fy, fz\] in eV/Angstrom.
    pub force: [f64; 3],
    /// Atomic species.
    pub species: AtomSpecies,
    /// Partial charge in elementary charge units.
    pub charge: f64,
    /// Layer index (0 = anode, 1 = electrolyte, 2 = cathode).
    pub layer: u8,
    /// Whether the atom is frozen (rigid lattice site).
    pub frozen: bool,
    /// Atom ID.
    pub id: usize,
}

impl LatticeAtom {
    /// Create a new lattice atom.
    pub fn new(pos: [f64; 3], species: AtomSpecies, layer: u8) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            force: [0.0; 3],
            species,
            charge: species.default_charge(),
            layer,
            frozen: false,
            id: 0,
        }
    }

    /// Kinetic energy of this atom in eV.
    pub fn kinetic_energy(&self) -> f64 {
        let v2: f64 = self.vel.iter().map(|v| v * v).sum();
        // 0.5 * m * v^2, but we need consistent units
        // mass in kg, vel in A/fs = 1e-10 m / 1e-15 s = 1e5 m/s
        let v_si = v2 * 1.0e10; // (m/s)^2
        let ke_j = 0.5 * self.species.mass() * v_si;
        ke_j / ELEMENTARY_CHARGE // convert J to eV
    }

    /// Distance to another atom in Angstrom.
    pub fn distance_to(&self, other: &LatticeAtom) -> f64 {
        let dx = self.pos[0] - other.pos[0];
        let dy = self.pos[1] - other.pos[1];
        let dz = self.pos[2] - other.pos[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Distance to another atom with periodic boundary conditions.
    pub fn distance_pbc(&self, other: &LatticeAtom, box_size: [f64; 3]) -> f64 {
        let mut r2 = 0.0;
        for (d, &bs) in box_size.iter().enumerate() {
            let mut dr = self.pos[d] - other.pos[d];
            dr -= (dr / bs).round() * bs;
            r2 += dr * dr;
        }
        r2.sqrt()
    }

    /// Reset forces to zero.
    pub fn zero_forces(&mut self) {
        self.force = [0.0; 3];
    }
}

// ---------------------------------------------------------------------------
// SolventType
// ---------------------------------------------------------------------------

/// Electrolyte solvent types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolventType {
    /// Ethylene carbonate.
    EC,
    /// Dimethyl carbonate.
    DMC,
    /// Diethyl carbonate.
    DEC,
    /// Propylene carbonate.
    PC,
    /// Fluoroethylene carbonate.
    FEC,
}

impl SolventType {
    /// Dielectric constant at 25C.
    pub fn dielectric_constant(&self) -> f64 {
        match self {
            Self::EC => 89.78,
            Self::DMC => 3.11,
            Self::DEC => 2.82,
            Self::PC => 64.92,
            Self::FEC => 78.4,
        }
    }

    /// Viscosity in mPa.s at 25C.
    pub fn viscosity(&self) -> f64 {
        match self {
            Self::EC => 1.90,
            Self::DMC => 0.59,
            Self::DEC => 0.75,
            Self::PC => 2.53,
            Self::FEC => 4.10,
        }
    }

    /// Molecular mass in g/mol.
    pub fn molecular_mass(&self) -> f64 {
        match self {
            Self::EC => 88.06,
            Self::DMC => 90.08,
            Self::DEC => 118.13,
            Self::PC => 102.09,
            Self::FEC => 106.05,
        }
    }
}

// ---------------------------------------------------------------------------
// ElectrolyteMolecule
// ---------------------------------------------------------------------------

/// Represents an electrolyte molecule or ion in the simulation.
#[derive(Debug, Clone)]
pub struct ElectrolyteMolecule {
    /// Position \[x, y, z\] in Angstrom.
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\] in Angstrom/fs.
    pub vel: [f64; 3],
    /// Solvent type.
    pub solvent: SolventType,
    /// Ionic conductivity in S/m.
    pub conductivity: f64,
    /// Local dielectric constant (effective).
    pub dielectric: f64,
    /// Salt concentration in mol/L.
    pub salt_concentration: f64,
    /// Molecule ID.
    pub id: usize,
}

impl ElectrolyteMolecule {
    /// Create a new electrolyte molecule.
    pub fn new(pos: [f64; 3], solvent: SolventType) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            solvent,
            conductivity: 1.0,
            dielectric: solvent.dielectric_constant(),
            salt_concentration: 1.0,
            id: 0,
        }
    }

    /// Nernst-Einstein ionic conductivity estimate.
    ///
    /// sigma = n * q^2 * D / (kB * T)
    pub fn nernst_einstein_conductivity(
        concentration_mol_l: f64,
        diffusivity: f64,
        temperature: f64,
    ) -> f64 {
        let n = concentration_mol_l * 1000.0 * AVOGADRO; // ions/m^3
        n * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * diffusivity / (BOLTZMANN * temperature)
    }

    /// Debye screening length in Angstrom.
    pub fn debye_length(&self, temperature: f64) -> f64 {
        let c = self.salt_concentration * 1000.0 * AVOGADRO; // ions/m^3
        if c < 1.0e-10 {
            return f64::INFINITY;
        }
        let eps = self.dielectric * EPSILON_0;
        let lambda_m =
            (eps * BOLTZMANN * temperature / (2.0 * c * ELEMENTARY_CHARGE.powi(2))).sqrt();
        lambda_m * 1.0e10 // metres to Angstrom
    }

    /// Stokes-Einstein diffusion coefficient estimate in m^2/s.
    pub fn stokes_einstein_diffusivity(
        radius_angstrom: f64,
        viscosity_pa_s: f64,
        temperature: f64,
    ) -> f64 {
        let r = radius_angstrom * 1.0e-10;
        BOLTZMANN * temperature / (6.0 * PI * viscosity_pa_s * r)
    }
}

// ---------------------------------------------------------------------------
// InterCalationSite
// ---------------------------------------------------------------------------

/// A site in the host lattice that can accommodate a Li ion.
#[derive(Debug, Clone)]
pub struct IntercalationSite {
    /// Position in the host lattice \[x, y, z\] in Angstrom.
    pub pos: [f64; 3],
    /// Whether the site is occupied by a Li ion.
    pub occupied: bool,
    /// Activation energy for Li hopping in eV.
    pub activation_energy: f64,
    /// Diffusion barrier height in eV.
    pub diffusion_barrier: f64,
    /// Intercalation voltage vs Li/Li+ in V.
    pub voltage: f64,
    /// Site energy in eV (relative to reference).
    pub site_energy: f64,
    /// Host lattice type identifier.
    pub host: HostLattice,
}

/// Host lattice type for intercalation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HostLattice {
    /// Graphite (LiC6).
    Graphite,
    /// Lithium iron phosphate (LiFePO4).
    LFP,
    /// Lithium cobalt oxide (LiCoO2).
    LCO,
    /// Lithium manganese oxide (LiMn2O4).
    LMO,
    /// Lithium nickel manganese cobalt oxide.
    NMC,
    /// Silicon anode.
    Silicon,
}

impl IntercalationSite {
    /// Create a new intercalation site.
    pub fn new(pos: [f64; 3], host: HostLattice) -> Self {
        let (activation, barrier, voltage) = match host {
            HostLattice::Graphite => (0.30, 0.22, 0.1),
            HostLattice::LFP => (0.55, 0.30, 3.4),
            HostLattice::LCO => (0.40, 0.28, 3.9),
            HostLattice::LMO => (0.45, 0.33, 4.0),
            HostLattice::NMC => (0.42, 0.26, 3.7),
            HostLattice::Silicon => (0.50, 0.40, 0.3),
        };
        Self {
            pos,
            occupied: false,
            activation_energy: activation,
            diffusion_barrier: barrier,
            voltage,
            site_energy: -activation * 0.5,
            host,
        }
    }

    /// Hop rate from this site at a given temperature (Arrhenius).
    ///
    /// k = nu_0 * exp(-E_a / (kB * T))
    pub fn hop_rate(&self, temperature: f64) -> f64 {
        let nu_0 = 1.0e13; // attempt frequency in Hz
        let kbt = BOLTZMANN * temperature / ELEMENTARY_CHARGE; // in eV
        nu_0 * (-self.activation_energy / kbt).exp()
    }

    /// Diffusivity estimate from hop rate and lattice spacing.
    ///
    /// D = (1/6) * a^2 * k   (for 3D random walk)
    pub fn diffusivity(&self, temperature: f64, lattice_spacing_angstrom: f64) -> f64 {
        let k = self.hop_rate(temperature);
        let a = lattice_spacing_angstrom * 1.0e-10; // to m
        (1.0 / 6.0) * a * a * k
    }

    /// Occupation probability at equilibrium (Fermi-Dirac).
    pub fn occupation_probability(&self, chemical_potential: f64, temperature: f64) -> f64 {
        let kbt = BOLTZMANN * temperature / ELEMENTARY_CHARGE; // eV
        if kbt < 1.0e-15 {
            return if self.site_energy < chemical_potential {
                1.0
            } else {
                0.0
            };
        }
        1.0 / (1.0 + ((self.site_energy - chemical_potential) / kbt).exp())
    }
}

// ---------------------------------------------------------------------------
// SeiLayer
// ---------------------------------------------------------------------------

/// Solid Electrolyte Interface (SEI) layer model.
///
/// The SEI forms on the anode surface during the first charge cycles
/// and controls battery aging.
#[derive(Debug, Clone)]
pub struct SeiLayer {
    /// SEI thickness in nm.
    pub thickness: f64,
    /// SEI growth rate constant in nm/sqrt(s).
    pub growth_rate_constant: f64,
    /// Ionic conductivity of SEI in S/m.
    pub ionic_conductivity: f64,
    /// SEI density in kg/m^3.
    pub density: f64,
    /// Activation energy for SEI growth in eV.
    pub activation_energy: f64,
    /// Total charge consumed in SEI formation in C/m^2.
    pub charge_consumed: f64,
    /// Whether SEI has passivated (reached steady state).
    pub passivated: bool,
    /// Maximum thickness for passivation in nm.
    pub passivation_thickness: f64,
}

impl SeiLayer {
    /// Create a new SEI layer model.
    pub fn new() -> Self {
        Self {
            thickness: 0.0,
            growth_rate_constant: 0.1,
            ionic_conductivity: 1.0e-6,
            density: 2100.0,
            activation_energy: 0.35,
            charge_consumed: 0.0,
            passivated: false,
            passivation_thickness: 50.0,
        }
    }

    /// SEI growth rate at a given temperature (parabolic growth law).
    ///
    /// dL/dt = k / (2*L) * exp(-Ea / (kB*T))
    pub fn growth_rate(&self, temperature: f64) -> f64 {
        let kbt = BOLTZMANN * temperature / ELEMENTARY_CHARGE; // eV
        let arrhenius = (-self.activation_energy / kbt).exp();
        if self.thickness < 0.01 {
            // Initial rapid growth
            return self.growth_rate_constant * arrhenius;
        }
        self.growth_rate_constant * arrhenius / (2.0 * self.thickness)
    }

    /// Grow SEI for a time interval dt (in seconds).
    pub fn grow(&mut self, dt: f64, temperature: f64) {
        if self.passivated {
            return;
        }
        let rate = self.growth_rate(temperature);
        self.thickness += rate * dt;
        // Charge consumed: assume 1 electron per Li consumed
        let volume_change = rate * dt * 1.0e-9; // nm to m
        self.charge_consumed += volume_change * self.density * FARADAY / 6.941e-3;

        if self.thickness >= self.passivation_thickness {
            self.passivated = true;
        }
    }

    /// Resistance contribution of SEI layer in Ohm m^2.
    pub fn resistance(&self) -> f64 {
        if self.ionic_conductivity < 1.0e-20 {
            return f64::INFINITY;
        }
        let l = self.thickness * 1.0e-9; // nm to m
        l / self.ionic_conductivity
    }

    /// Overpotential due to SEI layer at a given current density (A/m^2).
    pub fn overpotential(&self, current_density: f64) -> f64 {
        self.resistance() * current_density
    }

    /// Parabolic growth: thickness = k * sqrt(t).
    pub fn parabolic_thickness(k: f64, time: f64) -> f64 {
        k * time.sqrt()
    }

    /// Capacity loss due to SEI in Ah/m^2.
    pub fn capacity_loss(&self) -> f64 {
        self.charge_consumed / 3600.0
    }
}

impl Default for SeiLayer {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ChargeTransfer
// ---------------------------------------------------------------------------

/// Butler-Volmer charge-transfer kinetics model.
#[derive(Debug, Clone)]
pub struct ChargeTransfer {
    /// Exchange current density in A/m^2.
    pub exchange_current_density: f64,
    /// Anodic transfer coefficient (alpha_a).
    pub alpha_a: f64,
    /// Cathodic transfer coefficient (alpha_c).
    pub alpha_c: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Number of electrons transferred.
    pub n_electrons: f64,
}

impl ChargeTransfer {
    /// Create a new Butler-Volmer charge-transfer model.
    pub fn new(exchange_current_density: f64, temperature: f64) -> Self {
        Self {
            exchange_current_density,
            alpha_a: 0.5,
            alpha_c: 0.5,
            temperature,
            n_electrons: 1.0,
        }
    }

    /// Butler-Volmer current density in A/m^2.
    ///
    /// j = j0 * \[exp(alpha_a * F * eta / (RT)) - exp(-alpha_c * F * eta / (RT))\]
    pub fn current_density(&self, overpotential: f64) -> f64 {
        let f_rt = FARADAY / (GAS_CONSTANT * self.temperature);
        let anodic = (self.alpha_a * self.n_electrons * f_rt * overpotential).exp();
        let cathodic = (-self.alpha_c * self.n_electrons * f_rt * overpotential).exp();
        self.exchange_current_density * (anodic - cathodic)
    }

    /// Linearized Butler-Volmer for small overpotentials.
    ///
    /// j ~ j0 * (alpha_a + alpha_c) * F * eta / (RT)
    pub fn current_density_linear(&self, overpotential: f64) -> f64 {
        let f_rt = FARADAY / (GAS_CONSTANT * self.temperature);
        self.exchange_current_density
            * (self.alpha_a + self.alpha_c)
            * self.n_electrons
            * f_rt
            * overpotential
    }

    /// Charge-transfer resistance in Ohm m^2.
    ///
    /// Rct = RT / (n * F * j0)
    pub fn charge_transfer_resistance(&self) -> f64 {
        GAS_CONSTANT * self.temperature
            / (self.n_electrons * FARADAY * self.exchange_current_density)
    }

    /// Tafel slope in V/decade.
    pub fn tafel_slope(&self) -> f64 {
        2.303 * GAS_CONSTANT * self.temperature / (self.alpha_a * self.n_electrons * FARADAY)
    }

    /// Overpotential from current density (inverse Butler-Volmer, Tafel approx.).
    pub fn overpotential_from_current(&self, current: f64) -> f64 {
        if self.exchange_current_density.abs() < 1.0e-30 {
            return 0.0;
        }
        let ratio = current / self.exchange_current_density;
        let f_rt = FARADAY / (GAS_CONSTANT * self.temperature);
        // Simplified: use sinh^{-1} for symmetric case
        if (self.alpha_a - self.alpha_c).abs() < 0.01 {
            let alpha = self.alpha_a;
            (ratio / 2.0).asinh() / (alpha * self.n_electrons * f_rt)
        } else {
            // Newton iteration for asymmetric case
            let mut eta = 0.0;
            for _ in 0..20 {
                let j = self.current_density(eta);
                let dj = self.exchange_current_density
                    * f_rt
                    * (self.alpha_a
                        * self.n_electrons
                        * (self.alpha_a * self.n_electrons * f_rt * eta).exp()
                        + self.alpha_c
                            * self.n_electrons
                            * (-self.alpha_c * self.n_electrons * f_rt * eta).exp());
                if dj.abs() < 1.0e-30 {
                    break;
                }
                eta -= (j - current) / dj;
            }
            eta
        }
    }
}

// ---------------------------------------------------------------------------
// BatteryCell
// ---------------------------------------------------------------------------

/// Electrode type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElectrodeType {
    /// Anode (e.g., graphite).
    Anode,
    /// Cathode (e.g., LFP, NMC).
    Cathode,
}

/// Full battery cell model.
#[derive(Debug, Clone)]
pub struct BatteryCell {
    /// Anode host lattice.
    pub anode: HostLattice,
    /// Cathode host lattice.
    pub cathode: HostLattice,
    /// Separator thickness in um.
    pub separator_thickness: f64,
    /// Electrolyte solvent.
    pub electrolyte: SolventType,
    /// State of charge \[0, 1\].
    pub soc: f64,
    /// Cell voltage in V.
    pub voltage: f64,
    /// Nominal capacity in Ah.
    pub nominal_capacity: f64,
    /// Internal resistance in Ohm.
    pub internal_resistance: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Current in A (positive = discharge).
    pub current: f64,
    /// SEI layer on anode.
    pub sei: SeiLayer,
    /// Charge-transfer kinetics.
    pub charge_transfer: ChargeTransfer,
    /// Cycle count.
    pub cycle_count: u32,
}

impl BatteryCell {
    /// Create a new battery cell.
    pub fn new(anode: HostLattice, cathode: HostLattice, electrolyte: SolventType) -> Self {
        Self {
            anode,
            cathode,
            separator_thickness: 25.0,
            electrolyte,
            soc: 0.5,
            voltage: 3.7,
            nominal_capacity: 3.0,
            internal_resistance: 0.05,
            temperature: REF_TEMPERATURE,
            current: 0.0,
            sei: SeiLayer::new(),
            charge_transfer: ChargeTransfer::new(10.0, REF_TEMPERATURE),
            cycle_count: 0,
        }
    }

    /// Open-circuit voltage as function of SOC (simplified polynomial).
    pub fn open_circuit_voltage(&self) -> f64 {
        let x = self.soc;
        // Polynomial fit for typical NMC/graphite cell
        let v_cathode = 4.2 - 0.8 * (1.0 - x).powi(2) + 0.1 * x;
        let v_anode = 0.1 + 0.05 * x.powi(2);
        v_cathode - v_anode
    }

    /// Terminal voltage under load.
    pub fn terminal_voltage(&self) -> f64 {
        let ocv = self.open_circuit_voltage();
        let eta = self.current * self.internal_resistance;
        let sei_eta = self.sei.overpotential(self.current / 0.01); // rough area
        ocv - eta - sei_eta
    }

    /// Update state for a time step dt (seconds).
    pub fn step(&mut self, dt: f64) {
        // Update SOC
        let capacity_c = self.nominal_capacity * 3600.0; // Ah to C
        if capacity_c > 0.0 {
            self.soc -= self.current * dt / capacity_c;
            self.soc = self.soc.clamp(0.0, 1.0);
        }
        // Update voltage
        self.voltage = self.terminal_voltage();
        // Grow SEI
        self.sei.grow(dt, self.temperature);
    }

    /// Set discharge current (positive = discharge).
    pub fn set_current(&mut self, current: f64) {
        self.current = current;
    }

    /// C-rate for a given current.
    pub fn c_rate(&self) -> f64 {
        if self.nominal_capacity < 1.0e-10 {
            return 0.0;
        }
        self.current.abs() / self.nominal_capacity
    }

    /// Energy stored in Wh.
    pub fn stored_energy(&self) -> f64 {
        self.nominal_capacity * self.voltage * self.soc
    }

    /// Power output in W.
    pub fn power(&self) -> f64 {
        self.voltage * self.current
    }

    /// Discharge the cell at constant current until cutoff voltage.
    ///
    /// Returns (time, voltage) history.
    pub fn constant_current_discharge(
        &mut self,
        current: f64,
        cutoff_voltage: f64,
        dt: f64,
        max_steps: usize,
    ) -> Vec<[f64; 2]> {
        self.set_current(current);
        let mut history = Vec::new();
        let mut time = 0.0;
        for _ in 0..max_steps {
            self.step(dt);
            time += dt;
            history.push([time, self.voltage]);
            if self.voltage < cutoff_voltage || self.soc <= 0.0 {
                break;
            }
        }
        history
    }
}

// ---------------------------------------------------------------------------
// CycleDegradation
// ---------------------------------------------------------------------------

/// Cycle degradation model for capacity fade and impedance rise.
#[derive(Debug, Clone)]
pub struct CycleDegradation {
    /// Initial capacity in Ah.
    pub initial_capacity: f64,
    /// Current capacity in Ah.
    pub current_capacity: f64,
    /// Capacity fade rate per cycle (fraction).
    pub fade_rate: f64,
    /// Calendar aging rate per sqrt(day).
    pub calendar_rate: f64,
    /// Internal impedance in Ohm.
    pub impedance: f64,
    /// Impedance growth rate per cycle in Ohm.
    pub impedance_growth_rate: f64,
    /// Lithium plating rate in mol/m^2/cycle.
    pub lithium_plating_rate: f64,
    /// Cumulative lithium plated in mol/m^2.
    pub lithium_plated: f64,
    /// Number of cycles completed.
    pub cycles: u32,
    /// Calendar days elapsed.
    pub calendar_days: f64,
}

impl CycleDegradation {
    /// Create a new cycle degradation model.
    pub fn new(initial_capacity: f64) -> Self {
        Self {
            initial_capacity,
            current_capacity: initial_capacity,
            fade_rate: 2.0e-4,
            calendar_rate: 5.0e-4,
            impedance: 0.05,
            impedance_growth_rate: 1.0e-4,
            lithium_plating_rate: 0.0,
            lithium_plated: 0.0,
            cycles: 0,
            calendar_days: 0.0,
        }
    }

    /// Apply one charge/discharge cycle.
    pub fn apply_cycle(&mut self) {
        self.cycles += 1;
        // Capacity fade: linear + SEI-related sqrt(N)
        let cycle_fade = self.initial_capacity * self.fade_rate;
        let sqrt_fade =
            self.initial_capacity * self.calendar_rate * (1.0 / (self.cycles as f64).sqrt());
        self.current_capacity -= cycle_fade + sqrt_fade;
        self.current_capacity = self.current_capacity.max(0.0);

        // Impedance rise
        self.impedance += self.impedance_growth_rate;

        // Lithium plating
        self.lithium_plated += self.lithium_plating_rate;
    }

    /// Apply calendar aging for a given number of days.
    pub fn apply_calendar_aging(&mut self, days: f64) {
        self.calendar_days += days;
        let loss = self.initial_capacity * self.calendar_rate * self.calendar_days.sqrt();
        self.current_capacity = (self.initial_capacity - loss).max(0.0);
    }

    /// State of health (SOH) as fraction of initial capacity.
    pub fn state_of_health(&self) -> f64 {
        if self.initial_capacity < 1.0e-15 {
            return 0.0;
        }
        self.current_capacity / self.initial_capacity
    }

    /// Remaining useful life estimate (cycles until 80% SOH).
    pub fn remaining_useful_life(&self) -> Option<u32> {
        if self.fade_rate < 1.0e-15 {
            return None;
        }
        let soh = self.state_of_health();
        if soh <= 0.8 {
            return Some(0);
        }
        let capacity_to_lose = self.current_capacity - 0.8 * self.initial_capacity;
        let per_cycle = self.initial_capacity * self.fade_rate;
        if per_cycle < 1.0e-15 {
            return None;
        }
        Some((capacity_to_lose / per_cycle) as u32)
    }

    /// Capacity retention as a fraction after N cycles (empirical model).
    pub fn capacity_retention_model(n_cycles: u32, alpha: f64, beta: f64) -> f64 {
        let n = n_cycles as f64;
        1.0 - alpha * n - beta * n.sqrt()
    }
}

// ---------------------------------------------------------------------------
// BatteryAnalysis
// ---------------------------------------------------------------------------

/// Analysis tools for battery simulations.
#[derive(Debug, Clone)]
pub struct BatteryAnalysis {
    /// Discharge energy data: (C-rate, specific_energy_Wh_kg, specific_power_W_kg).
    pub ragone_data: Vec<[f64; 3]>,
    /// Coulombic efficiency per cycle.
    pub coulombic_efficiency: Vec<f64>,
    /// Mean squared displacement data: (time, msd).
    pub msd_data: Vec<[f64; 2]>,
    /// Voltage profiles: Vec of (time, voltage) arrays.
    pub voltage_profiles: Vec<Vec<[f64; 2]>>,
}

impl BatteryAnalysis {
    /// Create a new analysis container.
    pub fn new() -> Self {
        Self {
            ragone_data: Vec::new(),
            coulombic_efficiency: Vec::new(),
            msd_data: Vec::new(),
            voltage_profiles: Vec::new(),
        }
    }

    /// Record a Ragone point.
    pub fn add_ragone_point(&mut self, c_rate: f64, specific_energy: f64, specific_power: f64) {
        self.ragone_data
            .push([c_rate, specific_energy, specific_power]);
    }

    /// Record coulombic efficiency for a cycle.
    pub fn add_coulombic_efficiency(&mut self, efficiency: f64) {
        self.coulombic_efficiency.push(efficiency);
    }

    /// Record MSD data point.
    pub fn add_msd_point(&mut self, time: f64, msd: f64) {
        self.msd_data.push([time, msd]);
    }

    /// Compute diffusivity from MSD using Einstein relation: D = MSD / (6t).
    pub fn diffusivity_from_msd(&self) -> Option<f64> {
        if self.msd_data.len() < 2 {
            return None;
        }
        // Linear fit of MSD vs time
        let n = self.msd_data.len() as f64;
        let sum_t: f64 = self.msd_data.iter().map(|d| d[0]).sum();
        let sum_m: f64 = self.msd_data.iter().map(|d| d[1]).sum();
        let sum_tt: f64 = self.msd_data.iter().map(|d| d[0] * d[0]).sum();
        let sum_tm: f64 = self.msd_data.iter().map(|d| d[0] * d[1]).sum();
        let denom = n * sum_tt - sum_t * sum_t;
        if denom.abs() < 1.0e-30 {
            return None;
        }
        let slope = (n * sum_tm - sum_t * sum_m) / denom;
        // D = slope / 6  (3D diffusion)
        Some(slope / 6.0)
    }

    /// Mean coulombic efficiency.
    pub fn mean_coulombic_efficiency(&self) -> f64 {
        if self.coulombic_efficiency.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.coulombic_efficiency.iter().sum();
        sum / self.coulombic_efficiency.len() as f64
    }

    /// Compute MSD for a set of atoms between two snapshots.
    pub fn compute_msd(initial_positions: &[[f64; 3]], current_positions: &[[f64; 3]]) -> f64 {
        if initial_positions.len() != current_positions.len() || initial_positions.is_empty() {
            return 0.0;
        }
        let n = initial_positions.len() as f64;
        let sum: f64 = initial_positions
            .iter()
            .zip(current_positions.iter())
            .map(|(p0, p1)| {
                let dx = p1[0] - p0[0];
                let dy = p1[1] - p0[1];
                let dz = p1[2] - p0[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        sum / n
    }

    /// Specific energy from discharge data (Wh/kg).
    pub fn specific_energy(voltage_time: &[[f64; 2]], current: f64, cell_mass_kg: f64) -> f64 {
        if voltage_time.len() < 2 || cell_mass_kg < 1.0e-10 {
            return 0.0;
        }
        let mut energy_j = 0.0;
        for i in 1..voltage_time.len() {
            let dt = voltage_time[i][0] - voltage_time[i - 1][0];
            let v_avg = 0.5 * (voltage_time[i][1] + voltage_time[i - 1][1]);
            energy_j += v_avg * current * dt;
        }
        energy_j / (3600.0 * cell_mass_kg) // J -> Wh, per kg
    }

    /// Specific power from energy and discharge time.
    pub fn specific_power(specific_energy_wh_kg: f64, discharge_time_h: f64) -> f64 {
        if discharge_time_h < 1.0e-15 {
            return 0.0;
        }
        specific_energy_wh_kg / discharge_time_h
    }

    /// Arrhenius rate for diffusion: D = D0 * exp(-Ea / (kB * T)).
    pub fn arrhenius_diffusivity(d0: f64, activation_energy_ev: f64, temperature: f64) -> f64 {
        let kbt = BOLTZMANN * temperature / ELEMENTARY_CHARGE;
        d0 * (-activation_energy_ev / kbt).exp()
    }

    /// Generate Ragone plot data for a battery cell at multiple C-rates.
    pub fn generate_ragone_data(
        cell: &mut BatteryCell,
        c_rates: &[f64],
        cutoff_voltage: f64,
        cell_mass_kg: f64,
        dt: f64,
    ) -> Vec<[f64; 3]> {
        let mut data = Vec::new();
        for &c_rate in c_rates {
            let mut cell_copy = cell.clone();
            cell_copy.soc = 1.0;
            let current = c_rate * cell_copy.nominal_capacity;
            let max_steps = (3600.0 * 5.0 / dt) as usize;
            let history =
                cell_copy.constant_current_discharge(current, cutoff_voltage, dt, max_steps);
            if history.len() < 2 {
                continue;
            }
            let se = Self::specific_energy(&history, current, cell_mass_kg);
            let discharge_time_h = history.last().map_or(0.0, |h| h[0]) / 3600.0;
            let sp = Self::specific_power(se, discharge_time_h);
            data.push([c_rate, se, sp]);
        }
        data
    }
}

impl Default for BatteryAnalysis {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Lennard-Jones force between battery atoms
// ---------------------------------------------------------------------------

/// Compute LJ force magnitude between two atoms.
///
/// Uses Lorentz-Berthelot combining rules.
pub fn lj_force_magnitude(r: f64, sigma1: f64, epsilon1: f64, sigma2: f64, epsilon2: f64) -> f64 {
    let sigma = 0.5 * (sigma1 + sigma2);
    let epsilon = (epsilon1 * epsilon2).sqrt();
    if r < 0.1 * sigma || r > 3.0 * sigma {
        return 0.0;
    }
    let sr6 = (sigma / r).powi(6);
    let sr12 = sr6 * sr6;
    24.0 * epsilon * (2.0 * sr12 - sr6) / r
}

/// Compute Coulomb force magnitude between two charged atoms.
///
/// F = k_e * q1 * q2 / r^2   (in eV/Angstrom, with charges in e).
pub fn coulomb_force(r: f64, q1: f64, q2: f64, dielectric: f64) -> f64 {
    if r < 0.1 {
        return 0.0;
    }
    // Coulomb constant in eV*Angstrom / e^2
    let k_e = ELEMENTARY_CHARGE / (4.0 * PI * EPSILON_0 * 1.0e-10) / ELEMENTARY_CHARGE;
    k_e * q1 * q2 / (dielectric * r * r)
}

/// Compute pairwise forces on a list of battery atoms.
///
/// Updates force arrays in-place.
pub fn compute_battery_forces(atoms: &mut [LatticeAtom], cutoff: f64, dielectric: f64) {
    let n = atoms.len();
    // Zero forces
    for atom in atoms.iter_mut() {
        atom.zero_forces();
    }
    // O(N^2) pairwise loop (for small systems)
    for i in 0..n {
        for j in (i + 1)..n {
            let mut dr = [0.0f64; 3];
            let mut r2 = 0.0;
            for (d, dr_val) in dr.iter_mut().enumerate() {
                *dr_val = atoms[j].pos[d] - atoms[i].pos[d];
                r2 += *dr_val * *dr_val;
            }
            let r = r2.sqrt();
            if r > cutoff || r < 0.01 {
                continue;
            }
            // LJ force
            let flj = lj_force_magnitude(
                r,
                atoms[i].species.lj_sigma(),
                atoms[i].species.lj_epsilon(),
                atoms[j].species.lj_sigma(),
                atoms[j].species.lj_epsilon(),
            );
            // Coulomb force
            let fc = coulomb_force(r, atoms[i].charge, atoms[j].charge, dielectric);

            let f_total = flj + fc;
            for (d, &dr_val) in dr.iter().enumerate() {
                let f_d = f_total * dr_val / r;
                atoms[i].force[d] += f_d;
                atoms[j].force[d] -= f_d;
            }
        }
    }
}

/// Velocity Verlet integration step for battery atoms.
pub fn velocity_verlet_step(atoms: &mut [LatticeAtom], dt: f64) {
    // Convert eV/Angstrom to acceleration:
    // a = F / m, with F in eV/A, m in kg
    // 1 eV/A = 1.602e-19 / 1e-10 = 1.602e-9 N
    // a in A/fs^2: multiply by 1e-10 / (1e-15)^2 ... use direct unit factor
    let force_to_acc = ELEMENTARY_CHARGE / 1.0e-10; // N per eV/A

    for atom in atoms.iter_mut() {
        if atom.frozen {
            continue;
        }
        let m = atom.species.mass();
        for d in 0..3 {
            let a = atom.force[d] * force_to_acc / m; // m/s^2
            // Convert to A/fs^2: a * 1e-10 / (1e-15)^2 = a * 1e20
            // Actually we work in A/fs consistently
            let a_aft2 = a * 1.0e-10 * 1.0e-30; // A/fs^2
            atom.vel[d] += 0.5 * a_aft2 * dt;
            atom.pos[d] += atom.vel[d] * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1.0e-6;

    #[test]
    fn test_atom_species_mass_positive() {
        assert!(AtomSpecies::Li.mass() > 0.0);
        assert!(AtomSpecies::C.mass() > AtomSpecies::Li.mass());
    }

    #[test]
    fn test_lattice_atom_creation() {
        let a = LatticeAtom::new([0.0, 0.0, 0.0], AtomSpecies::Li, 0);
        assert_eq!(a.species, AtomSpecies::Li);
        assert!((a.charge - 1.0).abs() < TOL);
        assert_eq!(a.layer, 0);
    }

    #[test]
    fn test_lattice_atom_distance() {
        let a = LatticeAtom::new([0.0, 0.0, 0.0], AtomSpecies::Li, 0);
        let b = LatticeAtom::new([3.0, 4.0, 0.0], AtomSpecies::C, 0);
        assert!((a.distance_to(&b) - 5.0).abs() < TOL);
    }

    #[test]
    fn test_lattice_atom_distance_pbc() {
        let a = LatticeAtom::new([0.5, 0.0, 0.0], AtomSpecies::Li, 0);
        let b = LatticeAtom::new([9.5, 0.0, 0.0], AtomSpecies::Li, 0);
        let d = a.distance_pbc(&b, [10.0, 10.0, 10.0]);
        assert!((d - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_solvent_dielectric() {
        assert!(SolventType::EC.dielectric_constant() > SolventType::DMC.dielectric_constant());
    }

    #[test]
    fn test_electrolyte_debye_length() {
        let e = ElectrolyteMolecule::new([0.0; 3], SolventType::EC);
        let dl = e.debye_length(REF_TEMPERATURE);
        // Debye length should be a few Angstrom for 1M solution
        assert!(dl > 1.0);
        assert!(dl < 100.0);
    }

    #[test]
    fn test_nernst_einstein() {
        let sigma =
            ElectrolyteMolecule::nernst_einstein_conductivity(1.0, 1.0e-10, REF_TEMPERATURE);
        assert!(sigma > 0.0);
        assert!(sigma.is_finite());
    }

    #[test]
    fn test_stokes_einstein() {
        let d = ElectrolyteMolecule::stokes_einstein_diffusivity(2.0, 0.001, REF_TEMPERATURE);
        assert!(d > 0.0);
        assert!(d < 1.0e-8);
    }

    #[test]
    fn test_intercalation_site_graphite() {
        let s = IntercalationSite::new([0.0; 3], HostLattice::Graphite);
        assert!((s.activation_energy - 0.30).abs() < TOL);
        assert!(!s.occupied);
    }

    #[test]
    fn test_hop_rate_increases_with_temperature() {
        let s = IntercalationSite::new([0.0; 3], HostLattice::Graphite);
        let k1 = s.hop_rate(300.0);
        let k2 = s.hop_rate(400.0);
        assert!(k2 > k1);
    }

    #[test]
    fn test_intercalation_diffusivity_positive() {
        let s = IntercalationSite::new([0.0; 3], HostLattice::Graphite);
        let d = s.diffusivity(300.0, 3.35);
        assert!(d > 0.0);
    }

    #[test]
    fn test_occupation_probability_limits() {
        let s = IntercalationSite::new([0.0; 3], HostLattice::Graphite);
        // Very high chemical potential => occupied
        let p_high = s.occupation_probability(10.0, 300.0);
        assert!(p_high > 0.99);
        // Very low => empty
        let p_low = s.occupation_probability(-10.0, 300.0);
        assert!(p_low < 0.01);
    }

    #[test]
    fn test_sei_growth() {
        let mut sei = SeiLayer::new();
        assert!((sei.thickness - 0.0).abs() < TOL);
        sei.grow(3600.0, REF_TEMPERATURE);
        assert!(sei.thickness > 0.0);
    }

    #[test]
    fn test_sei_parabolic() {
        let t = SeiLayer::parabolic_thickness(0.1, 100.0);
        assert!((t - 1.0).abs() < TOL);
    }

    #[test]
    fn test_sei_resistance_increases_with_thickness() {
        let mut sei = SeiLayer::new();
        sei.thickness = 10.0;
        let r1 = sei.resistance();
        sei.thickness = 20.0;
        let r2 = sei.resistance();
        assert!(r2 > r1);
    }

    #[test]
    fn test_sei_passivation() {
        let mut sei = SeiLayer::new();
        sei.passivation_thickness = 1.0;
        sei.growth_rate_constant = 10.0; // speed up for test
        for _ in 0..100_000 {
            sei.grow(1000.0, REF_TEMPERATURE);
            if sei.passivated {
                break;
            }
        }
        assert!(sei.passivated);
    }

    #[test]
    fn test_butler_volmer_zero_overpotential() {
        let ct = ChargeTransfer::new(10.0, REF_TEMPERATURE);
        let j = ct.current_density(0.0);
        assert!(j.abs() < 1.0e-10);
    }

    #[test]
    fn test_butler_volmer_antisymmetric() {
        let ct = ChargeTransfer::new(10.0, REF_TEMPERATURE);
        let j_pos = ct.current_density(0.1);
        let j_neg = ct.current_density(-0.1);
        // Antisymmetric for symmetric alpha
        assert!((j_pos + j_neg).abs() < 1.0e-6);
    }

    #[test]
    fn test_butler_volmer_positive_overpotential() {
        let ct = ChargeTransfer::new(10.0, REF_TEMPERATURE);
        let j = ct.current_density(0.1);
        assert!(j > 0.0);
    }

    #[test]
    fn test_charge_transfer_resistance() {
        let ct = ChargeTransfer::new(10.0, REF_TEMPERATURE);
        let r = ct.charge_transfer_resistance();
        assert!(r > 0.0);
        // R = RT/(nF*j0)
        let expected = GAS_CONSTANT * REF_TEMPERATURE / (1.0 * FARADAY * 10.0);
        assert!((r - expected).abs() / expected < 1.0e-10);
    }

    #[test]
    fn test_tafel_slope() {
        let ct = ChargeTransfer::new(10.0, REF_TEMPERATURE);
        let ts = ct.tafel_slope();
        // ~0.118 V/decade at 25C for alpha=0.5, n=1
        assert!(ts > 0.1);
        assert!(ts < 0.15);
    }

    #[test]
    fn test_battery_cell_creation() {
        let cell = BatteryCell::new(HostLattice::Graphite, HostLattice::NMC, SolventType::EC);
        assert!((cell.soc - 0.5).abs() < TOL);
        assert!(cell.voltage > 3.0);
    }

    #[test]
    fn test_open_circuit_voltage_range() {
        let mut cell = BatteryCell::new(HostLattice::Graphite, HostLattice::NMC, SolventType::EC);
        cell.soc = 0.0;
        let v_low = cell.open_circuit_voltage();
        cell.soc = 1.0;
        let v_high = cell.open_circuit_voltage();
        assert!(v_high > v_low);
    }

    #[test]
    fn test_battery_cell_discharge_step() {
        let mut cell = BatteryCell::new(HostLattice::Graphite, HostLattice::NMC, SolventType::EC);
        cell.soc = 1.0;
        cell.set_current(1.0);
        cell.step(1.0);
        assert!(cell.soc < 1.0);
    }

    #[test]
    fn test_cycle_degradation() {
        let mut deg = CycleDegradation::new(3.0);
        for _ in 0..100 {
            deg.apply_cycle();
        }
        assert!(deg.state_of_health() < 1.0);
        assert!(deg.current_capacity < deg.initial_capacity);
    }

    #[test]
    fn test_capacity_retention_model() {
        let ret = CycleDegradation::capacity_retention_model(0, 1.0e-4, 1.0e-3);
        assert!((ret - 1.0).abs() < TOL);
        let ret100 = CycleDegradation::capacity_retention_model(100, 1.0e-4, 1.0e-3);
        assert!(ret100 < 1.0);
    }

    #[test]
    fn test_remaining_useful_life() {
        let mut deg = CycleDegradation::new(3.0);
        let rul = deg.remaining_useful_life();
        assert!(rul.is_some());
        assert!(rul.unwrap() > 0);
        // After many cycles
        for _ in 0..5000 {
            deg.apply_cycle();
        }
        // May be 0
        let _rul2 = deg.remaining_useful_life();
    }

    #[test]
    fn test_msd_diffusivity() {
        let mut analysis = BatteryAnalysis::new();
        // Linear MSD: MSD = 6 * D * t, with D = 1.0e-10
        for i in 0..100 {
            let t = i as f64 * 0.1;
            let msd = 6.0 * 1.0e-10 * t;
            analysis.add_msd_point(t, msd);
        }
        let d = analysis.diffusivity_from_msd().unwrap();
        assert!((d - 1.0e-10).abs() / 1.0e-10 < 0.01);
    }

    #[test]
    fn test_compute_msd() {
        let p0 = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let p1 = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let msd = BatteryAnalysis::compute_msd(&p0, &p1);
        assert!((msd - 1.0).abs() < TOL);
    }

    #[test]
    fn test_coulombic_efficiency_mean() {
        let mut a = BatteryAnalysis::new();
        a.add_coulombic_efficiency(0.99);
        a.add_coulombic_efficiency(0.995);
        let mean = a.mean_coulombic_efficiency();
        assert!((mean - 0.9925).abs() < 1.0e-10);
    }

    #[test]
    fn test_specific_energy_calculation() {
        let data = vec![[0.0, 3.7], [3600.0, 3.5]];
        let se = BatteryAnalysis::specific_energy(&data, 1.0, 0.1);
        assert!(se > 0.0);
    }

    #[test]
    fn test_arrhenius_diffusivity() {
        let d1 = BatteryAnalysis::arrhenius_diffusivity(1.0e-4, 0.3, 300.0);
        let d2 = BatteryAnalysis::arrhenius_diffusivity(1.0e-4, 0.3, 400.0);
        assert!(d2 > d1);
    }

    #[test]
    fn test_lj_force_repulsive_at_close_range() {
        let f = lj_force_magnitude(2.0, 3.4, 0.003, 3.4, 0.003);
        // At r < sigma, should be repulsive (positive force pushing apart)
        assert!(f > 0.0);
    }

    #[test]
    fn test_coulomb_force_like_charges() {
        let f = coulomb_force(5.0, 1.0, 1.0, 1.0);
        // Like charges repel (positive force)
        assert!(f > 0.0);
    }
}
