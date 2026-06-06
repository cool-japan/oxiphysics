//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use rand::RngExt;
use std::f64::consts::PI;

/// High-level driver for cryogenic SPH simulations.
///
/// Combines particle dynamics with phase transition, two-fluid model,
/// Kapitza resistance, and vortex filament tracking.
#[derive(Debug, Clone)]
pub struct CryogenicSphSimulation {
    /// SPH particles.
    pub particles: Vec<CryogenicParticle>,
    /// Equation of state.
    pub eos: TaitEquationCryogenic,
    /// Phase transition model.
    pub phase_model: CryogenicPhaseTransition,
    /// Kapitza resistance model.
    pub kapitza: KapitzaResistance,
    /// Active vortex filaments (for He-II).
    pub vortex_filaments: Vec<VortexFilament>,
    /// Current simulation time in seconds.
    pub time: f64,
    /// Time step in seconds.
    pub dt: f64,
    /// Fluid type.
    pub fluid: CryogenicFluid,
}
impl CryogenicSphSimulation {
    /// Create a simulation filled with liquid nitrogen particles on a lattice.
    pub fn new_ln2_box(n_particles: usize, box_size: f64) -> Self {
        let n_side = (n_particles as f64).cbrt().ceil() as usize;
        let spacing = box_size / n_side as f64;
        let mass = LN2_DENSITY * spacing.powi(3);
        let mut particles = Vec::with_capacity(n_particles);
        let mut id = 0;
        'outer: for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    if id >= n_particles {
                        break 'outer;
                    }
                    let mut p = CryogenicParticle::new(
                        [
                            (ix as f64 + 0.5) * spacing,
                            (iy as f64 + 0.5) * spacing,
                            (iz as f64 + 0.5) * spacing,
                        ],
                        mass,
                        LN2_NBP,
                        CryogenicFluid::LiquidNitrogen,
                    );
                    p.h = spacing * 1.2;
                    p.id = id;
                    particles.push(p);
                    id += 1;
                }
            }
        }
        Self {
            particles,
            eos: TaitEquationCryogenic::liquid_nitrogen(),
            phase_model: CryogenicPhaseTransition::liquid_nitrogen(),
            kapitza: KapitzaResistance::new(),
            vortex_filaments: Vec::new(),
            time: 0.0,
            dt: 1.0e-5,
            fluid: CryogenicFluid::LiquidNitrogen,
        }
    }
    /// Create a simulation for superfluid helium-4 with quantised vortices.
    pub fn new_he2_box(n_particles: usize, box_size: f64, temperature: f64) -> Self {
        let n_side = (n_particles as f64).cbrt().ceil() as usize;
        let spacing = box_size / n_side as f64;
        let mass = LHE_DENSITY * spacing.powi(3);
        let mut particles = Vec::with_capacity(n_particles);
        let mut id = 0;
        'outer: for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    if id >= n_particles {
                        break 'outer;
                    }
                    let fluid = if temperature < HE4_LAMBDA {
                        CryogenicFluid::SuperfluidHelium4
                    } else {
                        CryogenicFluid::LiquidHelium4
                    };
                    let mut p = CryogenicParticle::new(
                        [
                            (ix as f64 + 0.5) * spacing,
                            (iy as f64 + 0.5) * spacing,
                            (iz as f64 + 0.5) * spacing,
                        ],
                        mass,
                        temperature,
                        fluid,
                    );
                    p.h = spacing * 1.2;
                    p.id = id;
                    particles.push(p);
                    id += 1;
                }
            }
        }
        let vortex = VortexFilament::straight_line(box_size, 10);
        Self {
            particles,
            eos: TaitEquationCryogenic::liquid_helium4(),
            phase_model: CryogenicPhaseTransition::liquid_helium4(),
            kapitza: KapitzaResistance::new(),
            vortex_filaments: vec![vortex],
            time: 0.0,
            dt: 1.0e-6,
            fluid: CryogenicFluid::SuperfluidHelium4,
        }
    }
    /// Add thermal velocity perturbation to particles using Maxwell-Boltzmann.
    pub fn add_thermal_velocities(&mut self, temperature: f64) {
        let mut rng = rand::rng();
        let kb_per_m = BOLTZMANN / self.particles.first().map(|p| p.mass).unwrap_or(N2_MASS);
        let sigma = (kb_per_m * temperature).sqrt();
        for p in &mut self.particles {
            p.vel[0] = rng.random_range(-sigma..sigma);
            p.vel[1] = rng.random_range(-sigma..sigma);
            p.vel[2] = rng.random_range(-sigma..sigma);
        }
    }
    /// Update densities via kernel summation (simplified SPH).
    ///
    /// Uses the cubic spline kernel W(r, h).
    pub fn update_densities(&mut self) {
        let n = self.particles.len();
        let mut densities = vec![0.0_f64; n];
        for (i, d) in densities.iter_mut().enumerate() {
            let mut rho = 0.0;
            let hi = self.particles[i].h;
            for j in 0..n {
                let r = self.particles[i].distance_to(&self.particles[j]);
                let w = cubic_spline_kernel(r, hi);
                rho += self.particles[j].mass * w;
            }
            *d = rho.max(1.0);
        }
        for (i, p) in self.particles.iter_mut().enumerate() {
            p.density = densities[i];
            p.pressure = self.eos.pressure(densities[i]).max(0.0);
        }
    }
    /// Apply gravity and pressure gradient forces.
    pub fn apply_forces(&mut self) {
        let n = self.particles.len();
        let mut accs = vec![[0.0_f64; 3]; n];
        for (i, acc) in accs.iter_mut().enumerate() {
            acc[2] -= GRAVITY;
            let rho_i = self.particles[i].density;
            let p_i = self.particles[i].pressure;
            let hi = self.particles[i].h;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dx = self.particles[i].pos[0] - self.particles[j].pos[0];
                let dy = self.particles[i].pos[1] - self.particles[j].pos[1];
                let dz = self.particles[i].pos[2] - self.particles[j].pos[2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-10);
                let dw = cubic_spline_gradient(r, hi);
                let rho_j = self.particles[j].density;
                let p_j = self.particles[j].pressure;
                let pressure_term = p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j);
                let factor = -self.particles[j].mass * pressure_term * dw / r;
                acc[0] += factor * dx;
                acc[1] += factor * dy;
                acc[2] += factor * dz;
            }
        }
        for (i, p) in self.particles.iter_mut().enumerate() {
            p.acc = accs[i];
        }
    }
    /// Integrate positions and velocities using leapfrog.
    pub fn integrate(&mut self) {
        let dt = self.dt;
        for p in &mut self.particles {
            p.vel[0] += p.acc[0] * dt;
            p.vel[1] += p.acc[1] * dt;
            p.vel[2] += p.acc[2] * dt;
            p.pos[0] += p.vel[0] * dt;
            p.pos[1] += p.vel[1] * dt;
            p.pos[2] += p.vel[2] * dt;
        }
        self.time += dt;
    }
    /// Run N time steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.update_densities();
            self.apply_forces();
            self.integrate();
        }
    }
    /// Total kinetic energy of all particles in J.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }
    /// Mean temperature estimate from kinetic energy.
    ///
    /// T = 2 KE / (3 N k_B)
    pub fn mean_temperature(&self) -> f64 {
        let ke = self.total_kinetic_energy();
        let n = self.particles.len() as f64;
        if n < 1.0 {
            return 0.0;
        }
        2.0 * ke / (3.0 * n * BOLTZMANN)
    }
    /// Number of particles in superfluid state.
    pub fn superfluid_count(&self) -> usize {
        self.particles.iter().filter(|p| p.is_superfluid()).count()
    }
}
/// Nucleate and film boiling heat transfer model for cryogenic fluids.
///
/// Implements the Rohsenow correlation for nucleate boiling and the
/// Leidenfrost criterion for film boiling.
#[derive(Debug, Clone)]
pub struct BoilingHeatTransfer {
    /// Fluid type.
    pub fluid: CryogenicFluid,
    /// Saturation temperature in K.
    pub t_sat: f64,
    /// Saturation pressure in Pa.
    pub p_sat: f64,
    /// Rohsenow constant (depends on surface-fluid combination).
    pub rohsenow_c: f64,
    /// Rohsenow exponent (typically 1/3).
    pub rohsenow_n: f64,
    /// Wall superheat in K (T_wall - T_sat).
    pub wall_superheat: f64,
}
impl BoilingHeatTransfer {
    /// Create a boiling model for liquid nitrogen on a smooth copper surface.
    pub fn ln2_copper_surface() -> Self {
        Self {
            fluid: CryogenicFluid::LiquidNitrogen,
            t_sat: LN2_NBP,
            p_sat: 101_325.0,
            rohsenow_c: 0.013,
            rohsenow_n: 1.0 / 3.0,
            wall_superheat: 5.0,
        }
    }
    /// Rohsenow nucleate boiling correlation.
    ///
    /// q" = mu_l * h_fg * \[g(rho_l - rho_v) / sigma\]^0.5 * \[Cp * DeltaT / (C_s * h_fg * Pr^n)\]^3
    /// Simplified form used here: q" = h_nucleate * delta_T
    pub fn nucleate_boiling_flux(&self) -> f64 {
        let delta_t = self.wall_superheat;
        let cp = self.fluid.specific_heat();
        let h_fg = self.fluid.latent_heat();
        let mu_l = self.fluid.dynamic_viscosity().max(1.0e-12);
        let rho_l = self.fluid.reference_density();
        let rho_v = rho_l * 0.01;
        let sigma = LN2_SURFACE_TENSION;
        let bracket = (GRAVITY * (rho_l - rho_v) / sigma).sqrt();
        let pr = mu_l * cp / self.fluid.thermal_conductivity().max(1.0e-20);
        let factor = cp * delta_t / (self.rohsenow_c * h_fg * pr.powf(self.rohsenow_n));
        mu_l * h_fg * bracket * factor * factor * factor
    }
    /// Critical heat flux (Zuber correlation) in W/m^2.
    ///
    /// q"_max = 0.131 * h_fg * rho_v * \[sigma * g * (rho_l - rho_v) / rho_v^2\]^0.25
    pub fn critical_heat_flux(&self) -> f64 {
        let h_fg = self.fluid.latent_heat();
        let rho_l = self.fluid.reference_density();
        let rho_v = rho_l * 0.01;
        let sigma = LN2_SURFACE_TENSION;
        0.131 * h_fg * rho_v * (sigma * GRAVITY * (rho_l - rho_v) / (rho_v * rho_v)).powf(0.25)
    }
    /// Film boiling Leidenfrost point.
    ///
    /// Minimum temperature for stable film boiling (Berenson correlation):
    /// delta_T_min ~ 0.127 * (rho_l - rho_v) / rho_v * h_fg * \[sigma / g(rho_l - rho_v)\]^0.5
    ///   / (rho_l^0.5 * cp_v^{-1} * k_v + ...)
    /// Simplified empirical minimum superheat:
    pub fn leidenfrost_superheat(&self) -> f64 {
        let rho_l = self.fluid.reference_density();
        let rho_v = rho_l * 0.01;
        let sigma = LN2_SURFACE_TENSION;
        let h_fg = self.fluid.latent_heat();
        0.127 * (rho_l - rho_v) / rho_v * h_fg * (sigma / (GRAVITY * (rho_l - rho_v))).sqrt()
            / (rho_l * self.fluid.specific_heat() * 1.0e-3)
    }
}
/// Tait equation of state for cryogenic fluids.
///
/// Uses the modified Tait EOS: P = P0 + B * \[(rho/rho0)^gamma - 1\]
#[derive(Debug, Clone)]
pub struct TaitEquationCryogenic {
    /// Reference pressure P0 in Pa.
    pub p0: f64,
    /// Stiffness parameter B in Pa.
    pub bulk_modulus: f64,
    /// Polytropic exponent gamma.
    pub gamma: f64,
    /// Reference density rho0 in kg/m^3.
    pub rho0: f64,
    /// Speed of sound at reference state in m/s.
    pub c0: f64,
}
impl TaitEquationCryogenic {
    /// Create a Tait EOS tuned for liquid nitrogen.
    pub fn liquid_nitrogen() -> Self {
        let rho0 = LN2_DENSITY;
        let c0 = 855.0;
        let gamma = 7.0;
        let bulk_modulus = rho0 * c0 * c0 / gamma;
        Self {
            p0: 101_325.0,
            bulk_modulus,
            gamma,
            rho0,
            c0,
        }
    }
    /// Create a Tait EOS tuned for liquid helium-4.
    pub fn liquid_helium4() -> Self {
        let rho0 = LHE_DENSITY;
        let c0 = 238.0;
        let gamma = 7.0;
        let bulk_modulus = rho0 * c0 * c0 / gamma;
        Self {
            p0: 101_325.0,
            bulk_modulus,
            gamma,
            rho0,
            c0,
        }
    }
    /// Pressure from density via Tait EOS.
    pub fn pressure(&self, density: f64) -> f64 {
        self.p0 + self.bulk_modulus * ((density / self.rho0).powf(self.gamma) - 1.0)
    }
    /// Speed of sound at given density.
    pub fn sound_speed(&self, density: f64) -> f64 {
        let dp_drho = self.bulk_modulus * self.gamma / self.rho0
            * (density / self.rho0).powf(self.gamma - 1.0);
        dp_drho.sqrt()
    }
}
/// Kapitza (acoustic mismatch) thermal boundary resistance.
///
/// The Kapitza resistance R_K describes the thermal resistance at an
/// interface between He-II and a solid material.
///
/// R_K ≈ R0 * T^{-3} (Khalatnikov acoustic mismatch theory)
#[derive(Debug, Clone)]
pub struct KapitzaResistance {
    /// Prefactor R0 in m^2 K^4 / W.
    pub r0: f64,
    /// Temperature exponent (typically -3).
    pub exponent: f64,
}
impl KapitzaResistance {
    /// Create a Kapitza resistance model with default parameters.
    pub fn new() -> Self {
        Self::default()
    }
    /// Create with custom prefactor r0.
    pub fn with_prefactor(r0: f64) -> Self {
        Self { r0, exponent: -3.0 }
    }
    /// Thermal boundary resistance in m^2 K / W at temperature T.
    pub fn resistance(&self, temperature: f64) -> f64 {
        self.r0 * temperature.powf(self.exponent)
    }
    /// Kapitza conductance (inverse resistance) in W / (m^2 K).
    pub fn conductance(&self, temperature: f64) -> f64 {
        1.0 / self.resistance(temperature)
    }
    /// Heat flux across the interface given a temperature jump delta_T.
    ///
    /// Q = G_K * delta_T  where G_K = 1/R_K.
    pub fn heat_flux(&self, temperature: f64, delta_t: f64) -> f64 {
        self.conductance(temperature) * delta_t
    }
}
/// Cryogenic liquid-vapour phase transition model.
///
/// Uses a generalised van der Waals equation of state to track
/// liquid-vapour coexistence at cryogenic temperatures.
#[derive(Debug, Clone)]
pub struct CryogenicPhaseTransition {
    /// Critical temperature in K.
    pub t_critical: f64,
    /// Critical pressure in Pa.
    pub p_critical: f64,
    /// Critical density in kg/m^3.
    pub rho_critical: f64,
    /// van der Waals attractive constant a in J m^3 / kg^2.
    pub vdw_a: f64,
    /// van der Waals volume constant b in m^3 / kg.
    pub vdw_b: f64,
    /// Molecular mass in kg/mol.
    pub molar_mass: f64,
}
impl CryogenicPhaseTransition {
    /// Create a phase transition model for liquid nitrogen.
    pub fn liquid_nitrogen() -> Self {
        let t_c = 126.19;
        let p_c = 3.396e6;
        let rho_c = 313.3;
        let m = N2_MASS * AVOGADRO;
        let a_molar = 27.0 * GAS_CONSTANT * GAS_CONSTANT * t_c * t_c / (64.0 * p_c);
        let b_molar = GAS_CONSTANT * t_c / (8.0 * p_c);
        Self {
            t_critical: t_c,
            p_critical: p_c,
            rho_critical: rho_c,
            vdw_a: a_molar / (m * m),
            vdw_b: b_molar / m,
            molar_mass: m,
        }
    }
    /// Create a phase transition model for liquid helium-4.
    pub fn liquid_helium4() -> Self {
        let t_c = 5.20;
        let p_c = 2.274e5;
        let rho_c = 69.3;
        let m = HE4_MASS * AVOGADRO;
        let a_molar = 27.0 * GAS_CONSTANT * GAS_CONSTANT * t_c * t_c / (64.0 * p_c);
        let b_molar = GAS_CONSTANT * t_c / (8.0 * p_c);
        Self {
            t_critical: t_c,
            p_critical: p_c,
            rho_critical: rho_c,
            vdw_a: a_molar / (m * m),
            vdw_b: b_molar / m,
            molar_mass: m,
        }
    }
    /// Pressure from van der Waals EOS at given specific volume and temperature.
    ///
    /// P = RT/(v-b) - a/v^2  where v = 1/rho is specific volume.
    pub fn pressure_vdw(&self, density: f64, temperature: f64) -> f64 {
        let v = 1.0 / density;
        let r_specific = GAS_CONSTANT / self.molar_mass;

        r_specific * temperature / (v - self.vdw_b) - self.vdw_a / (v * v)
    }
    /// Saturation pressure by Clausius-Clapeyron approximation.
    ///
    /// P_sat = P_c * exp(-L * M / R * (1/T - 1/T_c))
    pub fn saturation_pressure(&self, temperature: f64, latent_heat: f64) -> f64 {
        let l_molar = latent_heat * self.molar_mass;
        let exponent = -l_molar / GAS_CONSTANT * (1.0 / temperature - 1.0 / self.t_critical);
        self.p_critical * exponent.exp()
    }
    /// Reduced temperature T* = T / T_c.
    pub fn reduced_temperature(&self, temperature: f64) -> f64 {
        temperature / self.t_critical
    }
    /// Reduced density rho* = rho / rho_c.
    pub fn reduced_density(&self, density: f64) -> f64 {
        density / self.rho_critical
    }
    /// Phase indicator: returns 0 (liquid), 1 (vapour), 0.5 (two-phase).
    pub fn phase_indicator(&self, density: f64, temperature: f64) -> f64 {
        let t_star = self.reduced_temperature(temperature);
        if t_star > 1.0 {
            return 0.5;
        }
        let rho_star = self.reduced_density(density);
        if rho_star > 1.2 {
            0.0
        } else if rho_star < 0.3 {
            1.0
        } else {
            0.5
        }
    }
}
/// Zero-boil-off (ZBO) insulation model for cryogenic storage tanks.
///
/// Models the multi-layer insulation (MLI) system combined with an
/// active cryocooler that removes the residual heat load.
#[derive(Debug, Clone)]
pub struct ZeroBoilOffInsulation {
    /// MLI total thickness in metres.
    pub mli_thickness: f64,
    /// Number of MLI shield layers.
    pub n_layers: u32,
    /// Tank outer diameter in metres.
    pub outer_diameter: f64,
    /// Tank inner diameter in metres.
    pub inner_diameter: f64,
    /// Warm boundary temperature in K.
    pub warm_temp: f64,
    /// Cold boundary (fluid) temperature in K.
    pub cold_temp: f64,
    /// Effective MLI thermal conductivity in W/(m K).
    pub mli_k_eff: f64,
    /// Cryocooler coefficient of performance (COP).
    pub cryocooler_cop: f64,
    /// Residual support structure heat leak in W.
    pub support_heat_leak: f64,
}
impl ZeroBoilOffInsulation {
    /// Create a ZBO system for a liquid hydrogen tank.
    pub fn liquid_hydrogen_tank(outer_diam: f64, inner_diam: f64) -> Self {
        Self {
            mli_thickness: 0.02,
            n_layers: 30,
            outer_diameter: outer_diam,
            inner_diameter: inner_diam,
            warm_temp: 295.0,
            cold_temp: 20.28,
            mli_k_eff: 5.0e-5,
            cryocooler_cop: 0.05,
            support_heat_leak: 1.0,
        }
    }
    /// Total MLI heat flux in W/m^2 (Fourier law).
    pub fn mli_heat_flux(&self) -> f64 {
        self.mli_k_eff * (self.warm_temp - self.cold_temp) / self.mli_thickness
    }
    /// Tank outer surface area in m^2.
    pub fn outer_area(&self) -> f64 {
        PI * self.outer_diameter * self.outer_diameter
    }
    /// Total heat load on the cryogenic fluid in W.
    pub fn total_heat_load(&self) -> f64 {
        self.mli_heat_flux() * self.outer_area() + self.support_heat_leak
    }
    /// Cryocooler electrical power required to maintain ZBO in W.
    pub fn cryocooler_power(&self) -> f64 {
        self.total_heat_load() / self.cryocooler_cop
    }
    /// Effective thermal resistance of the MLI in K/W.
    pub fn thermal_resistance(&self) -> f64 {
        let area = self.outer_area();
        self.mli_thickness / (self.mli_k_eff * area)
    }
    /// Radiation heat leak through MLI using Stefan-Boltzmann law in W/m^2.
    pub fn radiation_heat_flux(&self) -> f64 {
        let sigma = 5.670374419e-8;
        let emissivity = 0.03 / self.n_layers as f64;
        emissivity * sigma * (self.warm_temp.powi(4) - self.cold_temp.powi(4))
    }
}
/// Enumeration of supported cryogenic fluids.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CryogenicFluid {
    /// Liquid helium-4.
    LiquidHelium4,
    /// Superfluid helium-4 (He-II, below lambda point).
    SuperfluidHelium4,
    /// Liquid nitrogen.
    LiquidNitrogen,
    /// Liquid hydrogen.
    LiquidHydrogen,
    /// Liquid neon.
    LiquidNeon,
}
impl CryogenicFluid {
    /// Reference density in kg/m^3 at normal boiling point.
    pub fn reference_density(&self) -> f64 {
        match self {
            Self::LiquidHelium4 | Self::SuperfluidHelium4 => LHE_DENSITY,
            Self::LiquidNitrogen => LN2_DENSITY,
            Self::LiquidHydrogen => 70.9,
            Self::LiquidNeon => 1207.0,
        }
    }
    /// Normal boiling point in K.
    pub fn boiling_point(&self) -> f64 {
        match self {
            Self::LiquidHelium4 | Self::SuperfluidHelium4 => LHE_NBP,
            Self::LiquidNitrogen => LN2_NBP,
            Self::LiquidHydrogen => 20.28,
            Self::LiquidNeon => 27.07,
        }
    }
    /// Specific heat capacity at constant pressure in J/(kg K).
    pub fn specific_heat(&self) -> f64 {
        match self {
            Self::LiquidHelium4 | Self::SuperfluidHelium4 => LHE_CP,
            Self::LiquidNitrogen => LN2_CP,
            Self::LiquidHydrogen => 9_690.0,
            Self::LiquidNeon => 1_950.0,
        }
    }
    /// Latent heat of vaporisation in J/kg.
    pub fn latent_heat(&self) -> f64 {
        match self {
            Self::LiquidHelium4 | Self::SuperfluidHelium4 => LHE_LHV,
            Self::LiquidNitrogen => LN2_LHV,
            Self::LiquidHydrogen => 445_000.0,
            Self::LiquidNeon => 86_000.0,
        }
    }
    /// Thermal conductivity in W/(m K).
    pub fn thermal_conductivity(&self) -> f64 {
        match self {
            Self::LiquidHelium4 => 0.0186,
            Self::SuperfluidHelium4 => 1_000.0,
            Self::LiquidNitrogen => LN2_THERMAL_COND,
            Self::LiquidHydrogen => 0.0989,
            Self::LiquidNeon => 0.113,
        }
    }
    /// Dynamic viscosity in Pa s.
    pub fn dynamic_viscosity(&self) -> f64 {
        match self {
            Self::LiquidHelium4 => 3.5e-6,
            Self::SuperfluidHelium4 => 0.0,
            Self::LiquidNitrogen => 1.58e-4,
            Self::LiquidHydrogen => 1.32e-5,
            Self::LiquidNeon => 1.14e-4,
        }
    }
    /// Returns true if this fluid can exhibit superfluidity.
    pub fn is_superfluid(&self) -> bool {
        matches!(self, Self::SuperfluidHelium4)
    }
}
/// A single SPH particle representing a parcel of cryogenic fluid.
#[derive(Debug, Clone)]
pub struct CryogenicParticle {
    /// Position \[x, y, z\] in metres.
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\] in m/s.
    pub vel: [f64; 3],
    /// Acceleration \[ax, ay, az\] in m/s^2.
    pub acc: [f64; 3],
    /// Temperature in Kelvin.
    pub temperature: f64,
    /// Pressure in Pa.
    pub pressure: f64,
    /// Density in kg/m^3.
    pub density: f64,
    /// Mass in kg.
    pub mass: f64,
    /// Smoothing length in metres.
    pub h: f64,
    /// Thermal energy per unit mass in J/kg.
    pub internal_energy: f64,
    /// Fluid type.
    pub fluid: CryogenicFluid,
    /// Phase indicator: 0 = liquid, 1 = vapour, intermediate = two-phase.
    pub phase: f64,
    /// Particle ID.
    pub id: usize,
}
impl CryogenicParticle {
    /// Create a new cryogenic particle.
    pub fn new(pos: [f64; 3], mass: f64, temperature: f64, fluid: CryogenicFluid) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            acc: [0.0; 3],
            temperature,
            pressure: 101_325.0,
            density: fluid.reference_density(),
            mass,
            h: 0.01,
            internal_energy: fluid.specific_heat() * temperature,
            fluid,
            phase: 0.0,
            id: 0,
        }
    }
    /// Kinetic energy of this particle in J.
    pub fn kinetic_energy(&self) -> f64 {
        let v2: f64 = self.vel.iter().map(|v| v * v).sum();
        0.5 * self.mass * v2
    }
    /// Speed in m/s.
    pub fn speed(&self) -> f64 {
        let v2: f64 = self.vel.iter().map(|v| v * v).sum();
        v2.sqrt()
    }
    /// Distance to another particle in metres.
    pub fn distance_to(&self, other: &CryogenicParticle) -> f64 {
        let dx = self.pos[0] - other.pos[0];
        let dy = self.pos[1] - other.pos[1];
        let dz = self.pos[2] - other.pos[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Check if this particle is in a superfluid state (He-4 below lambda).
    pub fn is_superfluid(&self) -> bool {
        self.fluid.is_superfluid()
            || (matches!(self.fluid, CryogenicFluid::LiquidHelium4)
                && self.temperature < HE4_LAMBDA)
    }
    /// Normal fluid fraction (Landau two-fluid model).
    ///
    /// rho_n / rho = (T / T_lambda)^(5.6) for T < T_lambda.
    pub fn normal_fraction(&self) -> f64 {
        if self.temperature >= HE4_LAMBDA {
            return 1.0;
        }
        (self.temperature / HE4_LAMBDA).powf(5.6).min(1.0)
    }
    /// Superfluid fraction (rho_s / rho).
    pub fn superfluid_fraction(&self) -> f64 {
        1.0 - self.normal_fraction()
    }
}
/// Gross-Pitaevskii equation (GPE) solver for Bose-Einstein condensate dynamics.
///
/// The dimensionless GPE in 3D is:
///   i hbar d/dt psi = \[-hbar^2 / (2m) nabla^2 + V_ext + g |psi|^2\] psi
///
/// Here we solve on a uniform 1D spatial grid using split-step Fourier method
/// (time-split operator). The 1D version captures the essential physics.
#[derive(Debug, Clone)]
pub struct GrossPitaevskiiSolver {
    /// Number of grid points.
    pub n_grid: usize,
    /// Grid spacing in metres.
    pub dx: f64,
    /// Particle mass in kg.
    pub mass: f64,
    /// s-wave scattering length in metres.
    pub scattering_length: f64,
    /// Total number of particles.
    pub n_particles: f64,
    /// Real part of wave function psi on grid.
    pub psi_re: Vec<f64>,
    /// Imaginary part of wave function psi on grid.
    pub psi_im: Vec<f64>,
    /// External potential V_ext on grid in J.
    pub v_ext: Vec<f64>,
    /// Chemical potential in J.
    pub chemical_potential: f64,
    /// Coupling constant g = 4 pi hbar^2 a_s / m in J m^3.
    pub g_coupling: f64,
}
impl GrossPitaevskiiSolver {
    /// Create a new GPE solver on a uniform grid.
    ///
    /// Initialises the condensate in the Thomas-Fermi ground state.
    pub fn new(
        n_grid: usize,
        length: f64,
        mass: f64,
        scattering_length: f64,
        n_particles: f64,
    ) -> Self {
        let dx = length / n_grid as f64;
        let g = 4.0 * PI * HBAR * HBAR * scattering_length / mass;
        let mu =
            0.5 * (15.0 * n_particles * g * mass * 1.0 / (8.0 * PI)).powf(0.4) / (length * length);
        let mut psi_re = vec![0.0; n_grid];
        let mut psi_im = vec![0.0; n_grid];
        let mut v_ext = vec![0.0_f64; n_grid];
        for i in 0..n_grid {
            let x = (i as f64 - n_grid as f64 / 2.0) * dx;
            let v = 0.5 * mass * (2.0 * PI * 1000.0).powi(2) * x * x;
            v_ext[i] = v;
            let n_tf = ((mu - v) / g).max(0.0);
            psi_re[i] = n_tf.sqrt();
            psi_im[i] = 0.0;
        }
        let norm2: f64 = psi_re.iter().map(|r| r * r).sum::<f64>() * dx;
        let norm = (norm2 / n_particles).sqrt().max(1.0e-30);
        for r in &mut psi_re {
            *r /= norm;
        }
        Self {
            n_grid,
            dx,
            mass,
            scattering_length,
            n_particles,
            psi_re,
            psi_im,
            v_ext,
            chemical_potential: mu,
            g_coupling: g,
        }
    }
    /// Density n(x) = |psi(x)|^2 on the grid.
    pub fn density(&self) -> Vec<f64> {
        self.psi_re
            .iter()
            .zip(self.psi_im.iter())
            .map(|(r, i)| r * r + i * i)
            .collect()
    }
    /// Total norm N = integral |psi|^2 dx (should equal n_particles after normalisation).
    pub fn norm(&self) -> f64 {
        self.density().iter().sum::<f64>() * self.dx
    }
    /// Total energy of the condensate in J.
    ///
    /// E = integral \[hbar^2/(2m) |grad psi|^2 + V|psi|^2 + g/2 |psi|^4\] dx
    pub fn total_energy(&self) -> f64 {
        let n = self.n_grid;
        let mut energy = 0.0;
        for i in 0..n {
            let rho = self.psi_re[i] * self.psi_re[i] + self.psi_im[i] * self.psi_im[i];
            let im1 = if i > 0 { i - 1 } else { 0 };
            let ip1 = if i < n - 1 { i + 1 } else { n - 1 };
            let dpsi_re = (self.psi_re[ip1] - self.psi_re[im1]) / (2.0 * self.dx);
            let dpsi_im = (self.psi_im[ip1] - self.psi_im[im1]) / (2.0 * self.dx);
            let kin = HBAR * HBAR / (2.0 * self.mass) * (dpsi_re * dpsi_re + dpsi_im * dpsi_im);
            let pot = self.v_ext[i] * rho;
            let interact = 0.5 * self.g_coupling * rho * rho;
            energy += (kin + pot + interact) * self.dx;
        }
        energy
    }
    /// Advance the GPE by one imaginary time step (for ground-state relaxation).
    ///
    /// Uses a split-step method: potential step then kinetic step.
    pub fn imaginary_time_step(&mut self, dtau: f64) {
        let n = self.n_grid;
        let g = self.g_coupling;
        let v_effs: Vec<f64> = (0..n)
            .map(|i| {
                let rho = self.psi_re[i] * self.psi_re[i] + self.psi_im[i] * self.psi_im[i];
                self.v_ext[i] + g * rho
            })
            .collect();
        let v_mean = v_effs.iter().copied().sum::<f64>() / n as f64;
        for (i, &v_eff) in v_effs.iter().enumerate() {
            let exp_arg = -(v_eff - v_mean) * dtau / HBAR;
            let decay = exp_arg.clamp(-500.0, 500.0).exp();
            self.psi_re[i] *= decay;
            self.psi_im[i] *= decay;
        }
        let coeff = HBAR * HBAR / (2.0 * self.mass) * dtau / (self.dx * self.dx);
        let psi_re_old = self.psi_re.clone();
        let psi_im_old = self.psi_im.clone();
        for i in 1..(n - 1) {
            let lap_re = psi_re_old[i + 1] - 2.0 * psi_re_old[i] + psi_re_old[i - 1];
            let lap_im = psi_im_old[i + 1] - 2.0 * psi_im_old[i] + psi_im_old[i - 1];
            self.psi_re[i] += coeff * lap_re;
            self.psi_im[i] += coeff * lap_im;
        }
        let norm2: f64 = (self.psi_re.iter().map(|r| r * r).sum::<f64>()
            + self.psi_im.iter().map(|i| i * i).sum::<f64>())
            * self.dx;
        let norm = (norm2 / self.n_particles).sqrt().max(1.0e-30);
        for r in &mut self.psi_re {
            *r /= norm;
        }
        for im in &mut self.psi_im {
            *im /= norm;
        }
    }
    /// Healing length xi = hbar / sqrt(2 m g n0) in metres.
    ///
    /// Characterises the length scale over which the order parameter varies.
    pub fn healing_length(&self) -> f64 {
        let n0 = self.density().iter().cloned().fold(0.0_f64, f64::max);
        if n0 < 1.0e-30 {
            return f64::INFINITY;
        }
        HBAR / (2.0 * self.mass * self.g_coupling * n0).sqrt()
    }
    /// Speed of sound in the BEC: c = sqrt(g n0 / m) in m/s.
    pub fn speed_of_sound(&self) -> f64 {
        let n0 = self.density().iter().cloned().fold(0.0_f64, f64::max);
        (self.g_coupling * n0 / self.mass).sqrt()
    }
}
/// Landau two-fluid model for superfluid helium-4 (He-II).
///
/// He-II is modelled as two interpenetrating components:
/// - Normal fluid (viscous, carries entropy)
/// - Superfluid (inviscid, zero entropy)
///
/// The superfluid velocity is governed by the Euler equation with a
/// chemical potential gradient, while the normal fluid obeys the
/// Navier-Stokes equations.
#[derive(Debug, Clone)]
pub struct SuperfluidTwoFluidModel {
    /// Total density in kg/m^3.
    pub total_density: f64,
    /// Normal fluid fraction rho_n / rho.
    pub normal_fraction: f64,
    /// Normal fluid velocity \[vx, vy, vz\] in m/s.
    pub v_normal: [f64; 3],
    /// Superfluid velocity \[vx, vy, vz\] in m/s.
    pub v_super: [f64; 3],
    /// Temperature in K.
    pub temperature: f64,
    /// Specific entropy in J/(kg K).
    pub specific_entropy: f64,
    /// Second sound speed in m/s.
    pub second_sound_speed: f64,
}
impl SuperfluidTwoFluidModel {
    /// Create a new two-fluid model at a given temperature.
    ///
    /// Temperature must be below the lambda point (2.172 K).
    pub fn new(temperature: f64, total_density: f64) -> Self {
        let normal_fraction = if temperature >= HE4_LAMBDA {
            1.0
        } else {
            (temperature / HE4_LAMBDA).powf(5.6).min(1.0)
        };
        let second_sound = Self::compute_second_sound(temperature, normal_fraction, total_density);
        Self {
            total_density,
            normal_fraction,
            v_normal: [0.0; 3],
            v_super: [0.0; 3],
            temperature,
            specific_entropy: 756.0 * (temperature / HE4_LAMBDA).powf(3.0),
            second_sound_speed: second_sound,
        }
    }
    /// Second sound speed in He-II (Tisza formula).
    ///
    /// c2 = sqrt(rho_s * s^2 * T / (rho_n * cv))
    pub fn compute_second_sound(
        temperature: f64,
        normal_fraction: f64,
        _total_density: f64,
    ) -> f64 {
        let superfluid_fraction = 1.0 - normal_fraction;
        if normal_fraction < 1.0e-10 || superfluid_fraction < 1.0e-10 {
            return 0.0;
        }
        let entropy = 756.0 * (temperature / HE4_LAMBDA).powf(3.0);
        let cv = 3.0 * entropy;
        let numerator = superfluid_fraction * entropy * entropy * temperature;
        let denominator = normal_fraction * cv;
        if denominator < 1.0e-30 {
            return 0.0;
        }
        (numerator / denominator).sqrt()
    }
    /// Relative velocity (counterflow) in m/s.
    pub fn relative_velocity(&self) -> [f64; 3] {
        [
            self.v_normal[0] - self.v_super[0],
            self.v_normal[1] - self.v_super[1],
            self.v_normal[2] - self.v_super[2],
        ]
    }
    /// Mass flux of normal fluid in kg/(m^2 s).
    pub fn normal_flux(&self) -> [f64; 3] {
        let rho_n = self.total_density * self.normal_fraction;
        [
            rho_n * self.v_normal[0],
            rho_n * self.v_normal[1],
            rho_n * self.v_normal[2],
        ]
    }
    /// Entropy flux (heat flux divided by temperature) in W/(m^2 K).
    pub fn entropy_flux(&self) -> [f64; 3] {
        let rho_n = self.total_density * self.normal_fraction;
        let s_rho_n = self.specific_entropy * rho_n;
        [
            s_rho_n * self.v_normal[0],
            s_rho_n * self.v_normal[1],
            s_rho_n * self.v_normal[2],
        ]
    }
    /// Normal fluid density in kg/m^3.
    pub fn normal_density(&self) -> f64 {
        self.total_density * self.normal_fraction
    }
    /// Superfluid density in kg/m^3.
    pub fn superfluid_density(&self) -> f64 {
        self.total_density * (1.0 - self.normal_fraction)
    }
}
/// Available magnetocaloric materials for ADR cooling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MagnetocaloricMaterial {
    /// Ferric ammonium alum – commonly used for milli-Kelvin cooling.
    FerricAmmoniumAlum,
    /// Chrome potassium alum.
    ChromePotassiumAlum,
    /// Gadolinium – near room temperature MCE.
    Gadolinium,
    /// Gadolinium gallium garnet (GGG) – sub-Kelvin ADR.
    GadoliniumGalliumGarnet,
}
impl MagnetocaloricMaterial {
    /// Debye temperature in K.
    pub fn debye_temperature(&self) -> f64 {
        match self {
            Self::FerricAmmoniumAlum => 0.027,
            Self::ChromePotassiumAlum => 0.010,
            Self::Gadolinium => 169.0,
            Self::GadoliniumGalliumGarnet => 5.0,
        }
    }
    /// Molar spin quantum number J.
    pub fn spin_j(&self) -> f64 {
        match self {
            Self::FerricAmmoniumAlum => 2.5,
            Self::ChromePotassiumAlum => 1.5,
            Self::Gadolinium => 3.5,
            Self::GadoliniumGalliumGarnet => 3.5,
        }
    }
    /// Molar mass in kg/mol.
    pub fn molar_mass(&self) -> f64 {
        match self {
            Self::FerricAmmoniumAlum => 0.482,
            Self::ChromePotassiumAlum => 0.499,
            Self::Gadolinium => 0.15725,
            Self::GadoliniumGalliumGarnet => 0.9876,
        }
    }
}
/// Cryogenic tank slosh dynamics model.
///
/// Models the liquid motion in a partially filled cryogenic tank
/// under external accelerations (launch, manoeuvring). Uses a
/// simplified pendulum analogy for the dominant slosh mode.
#[derive(Debug, Clone)]
pub struct TankSloshDynamics {
    /// Tank inner radius in metres.
    pub tank_radius: f64,
    /// Fill level fraction \[0, 1\].
    pub fill_level: f64,
    /// Fluid type.
    pub fluid: CryogenicFluid,
    /// External acceleration vector \[ax, ay, az\] in m/s^2.
    pub external_acc: [f64; 3],
    /// Pendulum mass fraction.
    pub slosh_mass_fraction: f64,
    /// Pendulum length (effective) in metres.
    pub pendulum_length: f64,
    /// Slosh displacement in metres.
    pub slosh_displacement: f64,
    /// Slosh velocity in m/s.
    pub slosh_velocity: f64,
    /// Damping coefficient in s^{-1}.
    pub damping: f64,
}
impl TankSloshDynamics {
    /// Create a new tank slosh model for a cylindrical tank.
    pub fn cylindrical(radius: f64, fill_level: f64, fluid: CryogenicFluid) -> Self {
        let h = 2.0 * radius * fill_level;
        let pendulum_len = radius / (1.841_f64.tanh() * 1.841) * (1.841_f64 * h / radius).tanh();
        let xi = 1.841_f64 * h / radius;
        let slosh_fraction = (xi.tanh() / xi).min(1.0);
        Self {
            tank_radius: radius,
            fill_level,
            fluid,
            external_acc: [0.0, 0.0, -GRAVITY],
            slosh_mass_fraction: slosh_fraction,
            pendulum_length: pendulum_len.max(0.01),
            slosh_displacement: 0.0,
            slosh_velocity: 0.0,
            damping: 0.01,
        }
    }
    /// Natural slosh frequency in rad/s.
    pub fn natural_frequency(&self) -> f64 {
        let g_eff = self.external_acc.iter().map(|a| a * a).sum::<f64>().sqrt();
        (g_eff / self.pendulum_length).sqrt()
    }
    /// Slosh force on the tank wall in N, given total fluid mass.
    pub fn slosh_force(&self, total_mass: f64) -> f64 {
        let m_slosh = total_mass * self.slosh_mass_fraction;
        let omega0 = self.natural_frequency();
        m_slosh * omega0 * omega0 * self.slosh_displacement
    }
    /// Advance the slosh model by one time step.
    ///
    /// Uses a damped harmonic oscillator driven by lateral acceleration.
    pub fn step(&mut self, dt: f64, lateral_acc: f64) {
        let omega0 = self.natural_frequency();
        let omega0_sq = omega0 * omega0;
        let acc = -omega0_sq * self.slosh_displacement
            - 2.0 * self.damping * omega0 * self.slosh_velocity
            + lateral_acc;
        self.slosh_velocity += acc * dt;
        self.slosh_displacement += self.slosh_velocity * dt;
    }
    /// Surface wave height (first mode) at the tank wall in metres.
    pub fn wave_height(&self) -> f64 {
        self.slosh_displacement.abs()
    }
}
/// A quantised vortex filament in superfluid helium-4.
///
/// Each filament carries a single quantum of circulation
/// kappa = h/m = 9.97e-8 m^2/s for He-4.
///
/// The Biot-Savart law is used to compute the velocity field induced
/// by a set of filament segments.
#[derive(Debug, Clone)]
pub struct VortexFilament {
    /// Control points \[x, y, z\] defining the filament geometry in metres.
    pub points: Vec<[f64; 3]>,
    /// Circulation quantum in m^2/s (kappa = h/m for He-4).
    pub circulation: f64,
    /// Vortex core radius in metres (regularisation).
    pub core_radius: f64,
    /// Whether the filament is closed (ring) or open.
    pub is_closed: bool,
    /// Tension coefficient (energy per unit length) in J/m.
    pub tension: f64,
}
impl VortexFilament {
    /// Create a straight vortex filament along the z-axis from z=0 to z=L.
    pub fn straight_line(length: f64, n_points: usize) -> Self {
        let mut points = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let z = length * i as f64 / (n_points - 1) as f64;
            points.push([0.0, 0.0, z]);
        }
        let kappa = HBAR / HE4_MASS * 2.0 * PI;
        let a = 1.0e-10;
        let rho_s = LHE_DENSITY;
        let tension = 0.5 * rho_s * kappa * kappa / (2.0 * PI) * (1.0_f64 / a).ln();
        Self {
            points,
            circulation: kappa,
            core_radius: a,
            is_closed: false,
            tension,
        }
    }
    /// Create a circular vortex ring of radius R centred at origin in the xy-plane.
    pub fn vortex_ring(radius: f64, n_points: usize) -> Self {
        let mut points = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let theta = 2.0 * PI * i as f64 / n_points as f64;
            points.push([radius * theta.cos(), radius * theta.sin(), 0.0]);
        }
        let kappa = HBAR / HE4_MASS * 2.0 * PI;
        let a = 1.0e-10;
        let rho_s = LHE_DENSITY;
        let tension = 0.5 * rho_s * kappa * kappa / (2.0 * PI) * (1.0_f64 / a).ln();
        Self {
            points,
            circulation: kappa,
            core_radius: a,
            is_closed: true,
            tension,
        }
    }
    /// Compute the Biot-Savart induced velocity at a given point.
    ///
    /// v(r) = (kappa / 4pi) * integral ds x (r - s) / |r - s|^3
    /// regularised with core radius `a0`.
    pub fn induced_velocity(&self, r: [f64; 3]) -> [f64; 3] {
        let n = self.points.len();
        if n < 2 {
            return [0.0; 3];
        }
        let mut vel = [0.0; 3];
        let n_segs = if self.is_closed { n } else { n - 1 };
        for seg in 0..n_segs {
            let p1 = self.points[seg];
            let p2 = self.points[(seg + 1) % n];
            let dl = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
            let r_mid = [
                r[0] - 0.5 * (p1[0] + p2[0]),
                r[1] - 0.5 * (p1[1] + p2[1]),
                r[2] - 0.5 * (p1[2] + p2[2]),
            ];
            let rmag = (r_mid[0] * r_mid[0] + r_mid[1] * r_mid[1] + r_mid[2] * r_mid[2]).sqrt();
            let rmag_reg = (rmag * rmag + self.core_radius * self.core_radius).sqrt();
            let cross = [
                dl[1] * r_mid[2] - dl[2] * r_mid[1],
                dl[2] * r_mid[0] - dl[0] * r_mid[2],
                dl[0] * r_mid[1] - dl[1] * r_mid[0],
            ];
            let factor = self.circulation / (4.0 * PI * rmag_reg * rmag_reg * rmag_reg);
            vel[0] += factor * cross[0];
            vel[1] += factor * cross[1];
            vel[2] += factor * cross[2];
        }
        vel
    }
    /// Total filament length in metres.
    pub fn total_length(&self) -> f64 {
        let n = self.points.len();
        let n_segs = if self.is_closed { n } else { n - 1 };
        let mut len = 0.0;
        for seg in 0..n_segs {
            let p1 = self.points[seg];
            let p2 = self.points[(seg + 1) % n];
            let dl = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
            len += (dl[0] * dl[0] + dl[1] * dl[1] + dl[2] * dl[2]).sqrt();
        }
        len
    }
    /// Total elastic energy of the filament in J.
    pub fn elastic_energy(&self) -> f64 {
        self.tension * self.total_length()
    }
    /// Self-induced velocity of a vortex ring (translational).
    ///
    /// v_ring = (kappa / 4 pi R) * (ln(8 R / a) - 0.5) (Biot-Savart)
    pub fn ring_self_velocity(radius: f64) -> f64 {
        let kappa = QUANTUM_CIRCULATION;
        let a = 1.0e-10_f64;
        kappa / (4.0 * PI * radius) * ((8.0 * radius / a).ln() - 0.5)
    }
}
/// Magnetocaloric cooling model for magnetic refrigeration.
///
/// The magnetocaloric effect (MCE) is an adiabatic temperature change
/// upon application/removal of a magnetic field. Used for cooling
/// below 1 K via adiabatic demagnetisation refrigeration (ADR).
#[derive(Debug, Clone)]
pub struct MagnetocaloricCooling {
    /// Magnetic material name.
    pub material: MagnetocaloricMaterial,
    /// Current magnetic field in Tesla.
    pub b_field: f64,
    /// Current temperature in K.
    pub temperature: f64,
    /// Specific heat at zero field in J/(kg K).
    pub cv_zero: f64,
    /// Magnetic Gruneisen parameter.
    pub gamma_mag: f64,
}
impl MagnetocaloricCooling {
    /// Create a new magnetocaloric cooling stage with specified material and field.
    pub fn new(material: MagnetocaloricMaterial, b_field: f64, temperature: f64) -> Self {
        let cv = 10.0 * BOLTZMANN / material.molar_mass();
        Self {
            material,
            b_field,
            temperature,
            cv_zero: cv,
            gamma_mag: 0.5,
        }
    }
    /// Magnetic entropy of a spin-J system in a field B at temperature T.
    ///
    /// S_mag = R * ln(2J+1) - ... (Brillouin function based).
    /// Simplified using high-T approximation.
    pub fn magnetic_entropy(&self) -> f64 {
        let j = self.material.spin_j();
        let g = 2.0;
        let mu_b = 9.2740100783e-24;
        let x = g * j * mu_b * self.b_field / (BOLTZMANN * self.temperature);
        let two_j = 2.0 * j;
        let s_max = GAS_CONSTANT * (two_j + 1.0).ln() / self.material.molar_mass();
        if x < 0.01 {
            s_max
        } else {
            let bj = brillouin_function(j, x);
            s_max * (1.0 - bj)
        }
    }
    /// Adiabatic temperature change upon field removal (B -> 0).
    ///
    /// Delta_T = -T / Cv * (dS_mag/dB) * DeltaB
    pub fn adiabatic_temperature_change(&self, delta_b: f64) -> f64 {
        let j = self.material.spin_j();
        let g = 2.0;
        let mu_b = 9.2740100783e-24;
        let x = g * j * mu_b * self.b_field.max(1.0e-6) / (BOLTZMANN * self.temperature);
        let dbj_dx = brillouin_derivative(j, x);
        let ds_db = -GAS_CONSTANT / self.material.molar_mass() * g * j * mu_b
            / (BOLTZMANN * self.temperature)
            * dbj_dx;
        -self.temperature / self.cv_zero * ds_db * delta_b
    }
    /// Minimum temperature achievable (limited by nuclear hyperfine interactions).
    pub fn minimum_temperature(&self) -> f64 {
        self.material.debye_temperature()
    }
}
