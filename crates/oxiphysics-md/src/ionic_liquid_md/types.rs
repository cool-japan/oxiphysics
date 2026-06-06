//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::R_GAS;
use super::functions::*;

/// Parameters for the Born-Mayer-Huggins (BMH) potential.
///
/// V(r) = A · exp(B · (σ − r)) − C / r⁶ − D / r⁸
#[derive(Debug, Clone)]
pub struct BornMayerHugginsParams {
    /// Repulsive pre-exponential A (kJ/mol).
    pub a: f64,
    /// Repulsive steepness B (Å⁻¹).
    pub b: f64,
    /// Effective ion size σ (Å).
    pub sigma: f64,
    /// Dispersion coefficient C (kJ/mol·Å⁶).
    pub c: f64,
    /// Quadrupole dispersion D (kJ/mol·Å⁸).
    pub d: f64,
}
impl BornMayerHugginsParams {
    /// Create a new BMH parameter set.
    pub fn new(a: f64, b: f64, sigma: f64, c: f64, d: f64) -> Self {
        Self { a, b, sigma, c, d }
    }
    /// Example parameters for NaCl-like ionic liquid.
    pub fn nacl_like() -> Self {
        Self::new(4117.9, 3.1545, 2.34, 101.2, 48.2)
    }
    /// Evaluate potential energy at distance r.
    pub fn energy(&self, r: f64) -> f64 {
        if r < 0.1 {
            return 1e10;
        }
        let r6 = r.powi(6);
        let r8 = r.powi(8);
        self.a * (self.b * (self.sigma - r)).exp() - self.c / r6 - self.d / r8
    }
    /// Evaluate force magnitude at distance r (positive = repulsive).
    pub fn force(&self, r: f64) -> f64 {
        if r < 0.1 {
            return 1e10;
        }
        let r7 = r.powi(7);
        let r9 = r.powi(9);
        self.a * self.b * (self.b * (self.sigma - r)).exp() - 6.0 * self.c / r7 - 8.0 * self.d / r9
    }
}
/// Parameters for the Buckingham potential.
///
/// V(r) = A · exp(−r/ρ) − C / r⁶
#[derive(Debug, Clone)]
pub struct BuckinghamParams {
    /// Repulsive amplitude A (kJ/mol).
    pub a: f64,
    /// Repulsive length scale ρ (Å).
    pub rho: f64,
    /// Dispersion coefficient C (kJ/mol·Å⁶).
    pub c: f64,
}
impl BuckinghamParams {
    /// Create a new Buckingham parameter set.
    pub fn new(a: f64, rho: f64, c: f64) -> Self {
        Self { a, rho, c }
    }
    /// Example parameters for a generic ionic pair.
    pub fn generic_ionic() -> Self {
        Self::new(18000.0, 0.30, 120.0)
    }
    /// Evaluate potential energy at distance r.
    pub fn energy(&self, r: f64) -> f64 {
        if r < 0.1 {
            return 1e10;
        }
        self.a * (-r / self.rho).exp() - self.c / r.powi(6)
    }
    /// Evaluate force magnitude at distance r (positive = repulsive).
    pub fn force(&self, r: f64) -> f64 {
        if r < 0.1 {
            return 1e10;
        }
        (self.a / self.rho) * (-r / self.rho).exp() - 6.0 * self.c / r.powi(7)
    }
}
/// State of all ions in the simulation.
#[derive(Debug, Clone)]
pub struct IonState {
    /// Number of ions.
    pub n: usize,
    /// Positions (Å).
    pub positions: Vec<[f64; 3]>,
    /// Velocities (Å/ps).
    pub velocities: Vec<[f64; 3]>,
    /// Forces (kJ/(mol·Å)).
    pub forces: Vec<[f64; 3]>,
    /// Ion type for each particle.
    pub ion_types: Vec<IonType>,
    /// Charge of each ion (elementary charge units).
    pub charges: Vec<f64>,
    /// Mass of each ion (u).
    pub masses: Vec<f64>,
    /// LJ sigma for each ion (Å).
    pub lj_sigmas: Vec<f64>,
    /// LJ epsilon for each ion (kJ/mol).
    pub lj_epsilons: Vec<f64>,
}
impl IonState {
    /// Create a new empty ion state.
    pub fn new() -> Self {
        Self {
            n: 0,
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            ion_types: Vec::new(),
            charges: Vec::new(),
            masses: Vec::new(),
            lj_sigmas: Vec::new(),
            lj_epsilons: Vec::new(),
        }
    }
    /// Add an ion to the state.
    pub fn add_ion(&mut self, model: &IonModel, pos: [f64; 3], vel: [f64; 3]) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.forces.push([0.0; 3]);
        self.ion_types.push(model.ion_type);
        self.charges.push(model.charge);
        self.masses.push(model.mass);
        self.lj_sigmas.push(model.lj_sigma);
        self.lj_epsilons.push(model.lj_epsilon);
        self.n += 1;
    }
    /// Total number of cations.
    pub fn n_cations(&self) -> usize {
        self.ion_types
            .iter()
            .filter(|&&t| t == IonType::Cation)
            .count()
    }
    /// Total number of anions.
    pub fn n_anions(&self) -> usize {
        self.ion_types
            .iter()
            .filter(|&&t| t == IonType::Anion)
            .count()
    }
    /// Total charge of the system.
    pub fn total_charge(&self) -> f64 {
        self.charges.iter().sum()
    }
    /// Kinetic energy (kJ/mol).
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.n {
            let v2 = dot3(self.velocities[i], self.velocities[i]);
            ke += 0.5 * self.masses[i] * v2;
        }
        ke
    }
    /// Temperature from kinetic energy (K).
    pub fn temperature(&self) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        let ke = self.kinetic_energy();
        2.0 * ke * 1000.0 / (3.0 * self.n as f64 * R_GAS)
    }
    /// Center of mass position.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let mut com = [0.0; 3];
        let mut total_mass = 0.0;
        for i in 0..self.n {
            com = add3(com, scale3(self.positions[i], self.masses[i]));
            total_mass += self.masses[i];
        }
        if total_mass > 0.0 {
            scale3(com, 1.0 / total_mass)
        } else {
            com
        }
    }
    /// Zero all forces.
    pub fn zero_forces(&mut self) {
        for f in &mut self.forces {
            *f = [0.0; 3];
        }
    }
}
/// Result of a radial distribution function calculation.
#[derive(Debug, Clone)]
pub struct RdfResult {
    /// Bin centers (Å).
    pub r: Vec<f64>,
    /// g(r) values.
    pub g_r: Vec<f64>,
    /// Number of bins.
    pub n_bins: usize,
    /// Bin width (Å).
    pub dr: f64,
}
/// Result of ion cage analysis.
#[derive(Debug, Clone)]
pub struct CageDynamics {
    /// Average cage radius (Å) — mean distance to nearest neighbors.
    pub avg_cage_radius: f64,
    /// Cage rattling amplitude (Å) — RMS displacement within cage.
    pub rattling_amplitude: f64,
    /// Cage persistence time (ps) — how long ion stays in same cage.
    pub persistence_time: f64,
}
/// High-level ionic liquid simulation.
#[derive(Debug, Clone)]
pub struct IonicLiquidSimulation {
    /// Current ion state.
    pub state: IonState,
    /// Simulation parameters.
    pub params: IlSimParams,
    /// Current simulation time (ps).
    pub time: f64,
    /// Step counter.
    pub step: usize,
    /// Stored trajectory for analysis.
    pub trajectory: Vec<Vec<[f64; 3]>>,
    /// Stored current vectors for Green-Kubo.
    pub current_history: Vec<[f64; 3]>,
    /// Stored stress for viscosity.
    pub stress_history: Vec<f64>,
}
impl IonicLiquidSimulation {
    /// Create a new ionic liquid simulation.
    pub fn new(state: IonState, params: IlSimParams) -> Self {
        Self {
            state,
            params,
            time: 0.0,
            step: 0,
            trajectory: Vec::new(),
            current_history: Vec::new(),
            stress_history: Vec::new(),
        }
    }
    /// Run the simulation for a given number of steps.
    pub fn run(&mut self, n_steps: usize, save_interval: usize) {
        compute_forces(&mut self.state, &self.params);
        for _step in 0..n_steps {
            velocity_verlet_step(&mut self.state, &self.params);
            berendsen_thermostat(&mut self.state, &self.params);
            self.time += self.params.dt;
            self.step += 1;
            if save_interval > 0 && self.step.is_multiple_of(save_interval) {
                self.trajectory.push(self.state.positions.clone());
                let j = electric_current(&self.state.charges, &self.state.velocities);
                self.current_history.push(j);
                let vol = self.params.box_len.powi(3);
                let sxy = stress_tensor_element(&self.state, vol, 0, 1);
                self.stress_history.push(sxy);
            }
        }
    }
    /// Get the current potential energy.
    pub fn potential_energy(&self) -> f64 {
        compute_potential_energy(&self.state, &self.params)
    }
    /// Get the current kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.state.kinetic_energy()
    }
    /// Get the current temperature.
    pub fn temperature(&self) -> f64 {
        self.state.temperature()
    }
    /// Compute cation-anion RDF.
    pub fn compute_cation_anion_rdf(&self, n_bins: usize, r_max: f64) -> RdfResult {
        compute_rdf(
            &self.state,
            &self.params,
            PairType::CationAnion,
            n_bins,
            r_max,
        )
    }
    /// Compute conductivity via Green-Kubo.
    pub fn compute_conductivity(&self) -> f64 {
        if self.current_history.len() < 2 {
            return 0.0;
        }
        let max_lag = self.current_history.len() / 2;
        let cac = current_autocorrelation(&self.current_history, max_lag);
        let vol = self.params.box_len.powi(3);
        green_kubo_conductivity(&cac, self.params.dt, vol, self.params.temperature)
    }
    /// Compute viscosity via Green-Kubo.
    pub fn compute_viscosity(&self) -> f64 {
        if self.stress_history.len() < 2 {
            return 0.0;
        }
        let max_lag = self.stress_history.len() / 2;
        let sac = stress_autocorrelation(&self.stress_history, max_lag);
        let vol = self.params.box_len.powi(3);
        green_kubo_viscosity(&sac, self.params.dt, vol, self.params.temperature)
    }
}
/// Coarse-grained ion model.
#[derive(Debug, Clone)]
pub struct IonModel {
    /// Name or label for the ion.
    pub name: String,
    /// Ion type (cation or anion).
    pub ion_type: IonType,
    /// Formal charge in units of elementary charge.
    pub charge: f64,
    /// Mass in atomic mass units (u).
    pub mass: f64,
    /// Effective ionic radius (Å).
    pub radius: f64,
    /// Lennard-Jones sigma parameter (Å).
    pub lj_sigma: f64,
    /// Lennard-Jones epsilon parameter (kJ/mol).
    pub lj_epsilon: f64,
}
impl IonModel {
    /// Create a new ion model.
    pub fn new(
        name: &str,
        ion_type: IonType,
        charge: f64,
        mass: f64,
        radius: f64,
        lj_sigma: f64,
        lj_epsilon: f64,
    ) -> Self {
        Self {
            name: name.to_string(),
            ion_type,
            charge,
            mass,
            radius,
            lj_sigma,
            lj_epsilon,
        }
    }
    /// Create a generic imidazolium cation model (e.g. EMIM+).
    pub fn emim_cation() -> Self {
        Self::new("EMIM+", IonType::Cation, 1.0, 111.17, 3.0, 4.55, 1.05)
    }
    /// Create a generic BF4- anion model.
    pub fn bf4_anion() -> Self {
        Self::new("BF4-", IonType::Anion, -1.0, 86.81, 2.5, 3.95, 0.85)
    }
    /// Create a generic PF6- anion model.
    pub fn pf6_anion() -> Self {
        Self::new("PF6-", IonType::Anion, -1.0, 144.96, 2.7, 4.20, 0.95)
    }
    /// Create a lithium cation model.
    pub fn lithium_cation() -> Self {
        Self::new("Li+", IonType::Cation, 1.0, 6.941, 0.76, 2.13, 0.30)
    }
    /// Create a chloride anion model.
    pub fn chloride_anion() -> Self {
        Self::new("Cl-", IonType::Anion, -1.0, 35.453, 1.81, 4.42, 0.42)
    }
}
/// Parameters for Wolf summation of Coulomb interactions.
///
/// The Wolf method is a damped, shifted real-space Coulomb sum that converges
/// without Ewald reciprocal space.
#[derive(Debug, Clone)]
pub struct WolfParams {
    /// Damping parameter α (Å⁻¹).
    pub alpha: f64,
    /// Cutoff radius (Å).
    pub r_cut: f64,
    /// Conversion factor for charge-charge to kJ/mol·Å.
    pub coulomb_conv: f64,
}
impl WolfParams {
    /// Create Wolf summation parameters.
    pub fn new(alpha: f64, r_cut: f64) -> Self {
        Self {
            alpha,
            r_cut,
            coulomb_conv: 138.935_485,
        }
    }
    /// Default parameters for ionic liquids.
    pub fn default_il() -> Self {
        Self::new(0.2, 12.0)
    }
    /// Wolf potential between charges q_i and q_j at distance r.
    pub fn energy(&self, q_i: f64, q_j: f64, r: f64) -> f64 {
        if r >= self.r_cut || r < 0.01 {
            return 0.0;
        }
        let erfc_r = erfc_approx(self.alpha * r);
        let erfc_rc = erfc_approx(self.alpha * self.r_cut);
        self.coulomb_conv * q_i * q_j * (erfc_r / r - erfc_rc / self.r_cut)
    }
    /// Wolf force magnitude between charges q_i and q_j at distance r.
    pub fn force(&self, q_i: f64, q_j: f64, r: f64) -> f64 {
        if r >= self.r_cut || r < 0.01 {
            return 0.0;
        }
        let erfc_r = erfc_approx(self.alpha * r);
        let exp_r = (-self.alpha * self.alpha * r * r).exp();
        let erfc_rc = erfc_approx(self.alpha * self.r_cut);
        let force_r = erfc_r / (r * r) + 2.0 * self.alpha * exp_r / (PI.sqrt() * r);
        let shift = erfc_rc / (self.r_cut * self.r_cut);
        self.coulomb_conv * q_i * q_j * (force_r - shift)
    }
}
/// Pair type for RDF calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairType {
    /// Cation-anion pairs.
    CationAnion,
    /// Cation-cation pairs.
    CationCation,
    /// Anion-anion pairs.
    AnionAnion,
    /// All pairs.
    All,
}
/// Single data point for a Walden plot.
#[derive(Debug, Clone)]
pub struct WaldenPoint {
    /// Log₁₀ of molar conductivity (S·cm²/mol).
    pub log_lambda: f64,
    /// Log₁₀ of inverse viscosity (poise⁻¹).
    pub log_inv_eta: f64,
    /// Temperature (K).
    pub temperature: f64,
}
/// Simulation parameters for ionic liquid MD.
#[derive(Debug, Clone)]
pub struct IlSimParams {
    /// Box length for cubic periodic cell (Å).
    pub box_len: f64,
    /// LJ cutoff distance (Å).
    pub lj_cutoff: f64,
    /// Wolf summation parameters.
    pub wolf: WolfParams,
    /// Time step (ps).
    pub dt: f64,
    /// Temperature for thermostat (K).
    pub temperature: f64,
    /// Berendsen coupling constant (ps).
    pub tau_t: f64,
}
impl IlSimParams {
    /// Create simulation parameters with sensible defaults.
    pub fn default_params() -> Self {
        Self {
            box_len: 30.0,
            lj_cutoff: 10.0,
            wolf: WolfParams::default_il(),
            dt: 0.002,
            temperature: 400.0,
            tau_t: 0.5,
        }
    }
}
/// Type of ion in the ionic liquid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IonType {
    /// Positively charged ion.
    Cation,
    /// Negatively charged ion.
    Anion,
}
