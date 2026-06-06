// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geothermal fluid flow simulation using SPH.
//!
//! Covers:
//! - Hot spring and hydrothermal vent dynamics
//! - Supercritical fluid properties near the H₂O critical point
//! - Buoyancy-driven flow in porous rock (Darcy-Brinkman coupling)
//! - Hydrothermal alteration and mineral-rock interaction
//! - Steam-water phase separation (boiling two-phase flow)
//! - Heat pipe effect in the vadose zone
//! - Darcy flow coupled with SPH momentum
//! - Mineral precipitation and dissolution kinetics
//! - Enhanced Geothermal System (EGS) energy extraction
//! - Crustal thermal gradient and conductive heat flow

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Gravitational acceleration (m s⁻²).
const G: f64 = 9.81;

/// Universal gas constant (J mol⁻¹ K⁻¹).
const R_GAS: f64 = 8.314_462_618;

/// Critical temperature of water (K).
const T_CRIT_WATER: f64 = 647.096;

/// Critical pressure of water (Pa).
const P_CRIT_WATER: f64 = 22.064e6;

/// Critical density of water (kg m⁻³).
const RHO_CRIT_WATER: f64 = 322.0;

/// Molar mass of water (kg mol⁻¹).
const M_WATER: f64 = 0.018_015;

/// Latent heat of vaporization of water at 100 °C (J kg⁻¹).
const L_VAP: f64 = 2.257e6;

/// Reference atmospheric pressure (Pa).
const P_ATM: f64 = 101_325.0;

/// Typical rock thermal conductivity (W m⁻¹ K⁻¹).
const K_ROCK: f64 = 2.5;

/// Typical crustal geothermal gradient (K m⁻¹).
const GRAD_T_CRUST: f64 = 0.03;

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

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Normalize a 3-vector (returns zero vector if near-zero input).
#[inline]
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1.0e-30 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}

// ---------------------------------------------------------------------------
// SPH kernel
// ---------------------------------------------------------------------------

/// Cubic spline SPH kernel value W(r, h).
///
/// Returns the kernel weight for distance `r` and smoothing length `h`.
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Gradient of the cubic spline kernel: dW/dr · (r_vec / r).
///
/// Returns the kernel gradient vector for displacement `r_vec` and smoothing length `h`.
pub fn cubic_kernel_grad(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    let r = norm3(r_vec);
    if r < 1.0e-12 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dq: f64 = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * (-0.75 * t * t)
    } else {
        0.0
    };
    let dw_dr = dw_dq / h;
    scale3(normalize3(r_vec), dw_dr)
}

// ---------------------------------------------------------------------------
// Phase classification
// ---------------------------------------------------------------------------

/// Phase state of a geothermal fluid particle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluidPhase {
    /// Liquid water.
    LiquidWater,
    /// Steam (vapor).
    Steam,
    /// Two-phase mixture (boiling).
    TwoPhase,
    /// Supercritical fluid (T > T_crit, P > P_crit).
    Supercritical,
}

impl FluidPhase {
    /// Determine the phase from temperature (K) and pressure (Pa).
    pub fn from_tp(temperature: f64, pressure: f64) -> Self {
        if temperature > T_CRIT_WATER && pressure > P_CRIT_WATER {
            FluidPhase::Supercritical
        } else if temperature > boiling_point(pressure) {
            FluidPhase::Steam
        } else {
            FluidPhase::LiquidWater
        }
    }
}

/// Approximate boiling point of water at pressure `p` (Pa) using Antoine equation.
///
/// Returns temperature in K.
pub fn boiling_point(p: f64) -> f64 {
    // Antoine constants for water (T in °C, P in mmHg): A=8.07131, B=1730.63, C=233.426
    // Valid range 1–100 °C.  Convert Pa → mmHg (1 Pa = 0.007500617 mmHg).
    let p_mmhg = p * 0.007_500_617;
    let a = 8.07131;
    let b = 1730.63;
    let c = 233.426;
    let log_p = p_mmhg.log10();
    let t_c = b / (a - log_p) - c;
    t_c + 273.15
}

// ---------------------------------------------------------------------------
// Water equation of state
// ---------------------------------------------------------------------------

/// Simplified equation of state for water/steam.
///
/// Computes density (kg m⁻³) given temperature `t` (K) and pressure `p` (Pa).
pub fn water_density(t: f64, p: f64) -> f64 {
    let phase = FluidPhase::from_tp(t, p);
    match phase {
        FluidPhase::LiquidWater => {
            // Modified Tait equation
            let rho_ref = 999.84;
            let beta = 4.5e-10; // isothermal compressibility Pa⁻¹
            let alpha_t = 2.5e-4; // thermal expansion K⁻¹
            let t_ref = 277.15;
            rho_ref * (1.0 - alpha_t * (t - t_ref)) * (1.0 + beta * (p - P_ATM))
        }
        FluidPhase::Steam => {
            // Ideal gas approximation
            p * M_WATER / (R_GAS * t)
        }
        FluidPhase::TwoPhase => {
            let x = steam_quality(t, p);
            let rho_l = water_density(t - 1.0, p);
            let rho_g = p * M_WATER / (R_GAS * t);
            1.0 / (x / rho_g + (1.0 - x) / rho_l)
        }
        FluidPhase::Supercritical => {
            // Near critical point: use reduced properties
            let tau = T_CRIT_WATER / t;
            let pi = p / P_CRIT_WATER;
            RHO_CRIT_WATER * tau * pi * (1.0 + 0.1 * (tau - 1.0).powi(2))
        }
    }
}

/// Steam quality (vapor mass fraction) in the two-phase region.
///
/// `t` in K, `p` in Pa.
pub fn steam_quality(t: f64, p: f64) -> f64 {
    let t_sat = boiling_point(p);
    // Approximate: quality grows with superheat
    let dt = t - t_sat;
    (0.5 + dt * 0.05).clamp(0.0, 1.0)
}

/// Dynamic viscosity of water (Pa·s) as a function of temperature (K).
pub fn water_viscosity(t: f64) -> f64 {
    // Simplified Vogel equation
    let a = 1.856e-14;
    let b = 4209.0;
    let c = 0.04527;
    let d = -3.376e-5;
    a * (b / t + c + d * t).exp()
}

/// Specific heat capacity of water (J kg⁻¹ K⁻¹) as function of T (K).
pub fn water_cp(t: f64) -> f64 {
    // Polynomial fit
    let tc = t - 273.15;
    4217.6 - 3.82 * tc + 0.0169 * tc * tc
}

/// Thermal conductivity of water (W m⁻¹ K⁻¹) as function of T (K).
pub fn water_thermal_conductivity(t: f64) -> f64 {
    let tc = t - 273.15;
    0.5560 + 1.878e-3 * tc - 7.11e-6 * tc * tc
}

// ---------------------------------------------------------------------------
// Supercritical fluid
// ---------------------------------------------------------------------------

/// Supercritical fluid properties near the H₂O critical point.
///
/// Stores reduced temperature, pressure, and resulting thermophysical properties.
#[derive(Debug, Clone)]
pub struct SupercriticalFluid {
    /// Temperature (K).
    pub temperature: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Density (kg m⁻³).
    pub density: f64,
    /// Dynamic viscosity (Pa s).
    pub viscosity: f64,
    /// Specific heat capacity (J kg⁻¹ K⁻¹).
    pub cp: f64,
    /// Thermal conductivity (W m⁻¹ K⁻¹).
    pub thermal_conductivity: f64,
    /// Isothermal compressibility (Pa⁻¹).
    pub compressibility: f64,
}

impl SupercriticalFluid {
    /// Create a new supercritical fluid state at given temperature and pressure.
    pub fn new(temperature: f64, pressure: f64) -> Self {
        let tau = T_CRIT_WATER / temperature;
        let pi = pressure / P_CRIT_WATER;
        let density = RHO_CRIT_WATER * tau * pi * (1.0 + 0.1 * (tau - 1.0).powi(2));
        // Enhanced transport near critical point
        let viscosity = 1.0e-4 * (1.0 + 0.3 * ((tau - 1.0).powi(2) + (pi - 1.0).powi(2)).sqrt());
        let cp = 4500.0 * (1.0 + 5.0 * (-(tau - 1.0).powi(2) / 0.01).exp());
        let thermal_conductivity = 0.6 + 1.5 * (-(tau - 1.0).powi(2) / 0.05).exp();
        let compressibility = 1.0 / (P_CRIT_WATER * pi * (1.0 + (tau - 1.0).abs()));
        SupercriticalFluid {
            temperature,
            pressure,
            density,
            viscosity,
            cp,
            thermal_conductivity,
            compressibility,
        }
    }

    /// Check whether this fluid is in the supercritical state.
    pub fn is_supercritical(&self) -> bool {
        self.temperature > T_CRIT_WATER && self.pressure > P_CRIT_WATER
    }

    /// Reduced temperature τ = T_c / T.
    pub fn reduced_temperature(&self) -> f64 {
        T_CRIT_WATER / self.temperature
    }

    /// Reduced pressure π = P / P_c.
    pub fn reduced_pressure(&self) -> f64 {
        self.pressure / P_CRIT_WATER
    }
}

// ---------------------------------------------------------------------------
// Porous rock model
// ---------------------------------------------------------------------------

/// Rock type for porous media.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RockType {
    /// Granite (low permeability).
    Granite,
    /// Sandstone (high permeability).
    Sandstone,
    /// Basalt (intermediate).
    Basalt,
    /// Fractured rock (enhanced permeability).
    FracturedRock,
    /// Limestone (carbonate, reactive).
    Limestone,
}

impl RockType {
    /// Intrinsic permeability (m²) of the rock type.
    pub fn permeability(self) -> f64 {
        match self {
            RockType::Granite => 1.0e-18,
            RockType::Sandstone => 1.0e-13,
            RockType::Basalt => 1.0e-15,
            RockType::FracturedRock => 1.0e-12,
            RockType::Limestone => 5.0e-14,
        }
    }

    /// Porosity of the rock type (dimensionless).
    pub fn porosity(self) -> f64 {
        match self {
            RockType::Granite => 0.01,
            RockType::Sandstone => 0.25,
            RockType::Basalt => 0.10,
            RockType::FracturedRock => 0.05,
            RockType::Limestone => 0.15,
        }
    }

    /// Thermal conductivity of rock matrix (W m⁻¹ K⁻¹).
    pub fn thermal_conductivity(self) -> f64 {
        match self {
            RockType::Granite => 3.2,
            RockType::Sandstone => 2.0,
            RockType::Basalt => 1.8,
            RockType::FracturedRock => 2.5,
            RockType::Limestone => 2.8,
        }
    }

    /// Rock density (kg m⁻³).
    pub fn density(self) -> f64 {
        match self {
            RockType::Granite => 2700.0,
            RockType::Sandstone => 2200.0,
            RockType::Basalt => 2900.0,
            RockType::FracturedRock => 2600.0,
            RockType::Limestone => 2500.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Darcy flow
// ---------------------------------------------------------------------------

/// Darcy flow velocity in porous media.
///
/// Returns the Darcy velocity (m s⁻¹) given permeability `k` (m²),
/// fluid viscosity `mu` (Pa·s), and pressure gradient `grad_p` (Pa m⁻¹).
pub fn darcy_velocity(k: f64, mu: f64, grad_p: [f64; 3]) -> [f64; 3] {
    scale3(grad_p, -k / mu)
}

/// Darcy-Brinkman extended Darcy velocity including viscous term.
///
/// Adds Brinkman correction for the effective viscosity `mu_eff` in porous media.
pub fn darcy_brinkman_velocity(
    k: f64,
    mu: f64,
    mu_eff: f64,
    grad_p: [f64; 3],
    lap_u: [f64; 3],
    porosity: f64,
) -> [f64; 3] {
    let darcy = darcy_velocity(k, mu, grad_p);
    let brinkman = scale3(lap_u, mu_eff * k / (mu * porosity));
    add3(darcy, brinkman)
}

/// Compute the buoyancy force per unit volume (N m⁻³) for hydrothermal flow.
///
/// Uses Boussinesq approximation with reference density `rho_ref` (kg m⁻³),
/// thermal expansion coefficient `beta_t` (K⁻¹), and temperature excess `dt` (K).
pub fn buoyancy_force(rho_ref: f64, beta_t: f64, dt: f64) -> [f64; 3] {
    // Boussinesq approximation: hot fluid (dt > 0) is lighter and rises (+y direction).
    // f = rho_ref * beta_t * dt * g  (upward)
    [0.0, rho_ref * G * beta_t * dt, 0.0]
}

// ---------------------------------------------------------------------------
// Geothermal particle
// ---------------------------------------------------------------------------

/// A single SPH particle in the geothermal simulation.
#[derive(Debug, Clone)]
pub struct GeoParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m s⁻¹).
    pub vel: [f64; 3],
    /// SPH force accumulator (N kg⁻¹ = m s⁻²).
    pub acc: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg m⁻³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Enthalpy (J kg⁻¹).
    pub enthalpy: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Phase state.
    pub phase: FluidPhase,
    /// Mineral concentration (kg kg⁻¹ fluid).
    pub mineral_conc: f64,
    /// Dissolved silica concentration.
    pub silica_conc: f64,
    /// Porosity of host rock at this particle location.
    pub porosity: f64,
    /// Steam quality (vapor fraction) if two-phase.
    pub steam_quality: f64,
}

impl GeoParticle {
    /// Create a new geothermal particle at rest.
    pub fn new(pos: [f64; 3], mass: f64, temperature: f64, pressure: f64) -> Self {
        let phase = FluidPhase::from_tp(temperature, pressure);
        let density = water_density(temperature, pressure);
        let enthalpy = water_cp(temperature) * temperature;
        GeoParticle {
            pos,
            vel: [0.0; 3],
            acc: [0.0; 3],
            mass,
            density,
            pressure,
            temperature,
            enthalpy,
            h: 0.05,
            phase,
            mineral_conc: 0.0,
            silica_conc: 0.0,
            porosity: 0.1,
            steam_quality: 0.0,
        }
    }

    /// Update the particle's phase and thermodynamic state.
    pub fn update_phase(&mut self) {
        self.phase = FluidPhase::from_tp(self.temperature, self.pressure);
        self.density = water_density(self.temperature, self.pressure);
        self.enthalpy = water_cp(self.temperature) * self.temperature;
        if self.phase == FluidPhase::TwoPhase {
            self.steam_quality = steam_quality(self.temperature, self.pressure);
        }
    }

    /// Kinetic energy (J) of this particle.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }

    /// Thermal energy (J) of this particle.
    pub fn thermal_energy(&self) -> f64 {
        water_cp(self.temperature) * self.mass * self.temperature
    }
}

// ---------------------------------------------------------------------------
// Thermal gradient / crustal heat
// ---------------------------------------------------------------------------

/// Crustal thermal gradient model.
///
/// Computes temperature at depth `z` (m below surface) given surface temperature
/// `t_surface` (K) and gradient `grad` (K m⁻¹).
#[derive(Debug, Clone)]
pub struct CrustalThermalGradient {
    /// Surface temperature (K).
    pub t_surface: f64,
    /// Geothermal gradient (K m⁻¹).
    pub gradient: f64,
    /// Crustal radiogenic heat production (W m⁻³).
    pub heat_production: f64,
}

impl CrustalThermalGradient {
    /// Create a default continental crustal geothermal model.
    pub fn continental() -> Self {
        CrustalThermalGradient {
            t_surface: 288.15,
            gradient: GRAD_T_CRUST,
            heat_production: 1.5e-6,
        }
    }

    /// Create a volcanic/rift zone gradient model.
    pub fn volcanic() -> Self {
        CrustalThermalGradient {
            t_surface: 288.15,
            gradient: 0.08,
            heat_production: 3.0e-6,
        }
    }

    /// Temperature at depth `z` (m, positive downward).
    pub fn temperature_at_depth(&self, z: f64) -> f64 {
        self.t_surface + self.gradient * z
    }

    /// Conductive heat flux (W m⁻²) at the surface.
    pub fn surface_heat_flux(&self) -> f64 {
        K_ROCK * self.gradient
    }

    /// Lithostatic pressure (Pa) at depth `z` (m).
    pub fn lithostatic_pressure(&self, z: f64, rho_rock: f64) -> f64 {
        P_ATM + rho_rock * G * z
    }
}

// ---------------------------------------------------------------------------
// Hydrothermal vent
// ---------------------------------------------------------------------------

/// Hydrothermal vent or hot spring source term.
#[derive(Debug, Clone)]
pub struct HydrothermalVent {
    /// Position of the vent orifice (m).
    pub position: [f64; 3],
    /// Exit fluid temperature (K).
    pub temperature: f64,
    /// Exit fluid velocity (m s⁻¹).
    pub exit_velocity: f64,
    /// Vent radius (m).
    pub radius: f64,
    /// Salinity of vent fluid (kg NaCl per kg fluid).
    pub salinity: f64,
    /// Mineral loading (kg mineral per kg fluid).
    pub mineral_loading: f64,
    /// Whether this is a black smoker (T > 350°C) or white smoker.
    pub is_black_smoker: bool,
}

impl HydrothermalVent {
    /// Create a typical mid-ocean ridge black smoker vent.
    pub fn black_smoker(position: [f64; 3]) -> Self {
        HydrothermalVent {
            position,
            temperature: 623.15, // ~350°C
            exit_velocity: 1.0,
            radius: 0.05,
            salinity: 0.034,
            mineral_loading: 0.001,
            is_black_smoker: true,
        }
    }

    /// Create a hot spring vent (lower temperature).
    pub fn hot_spring(position: [f64; 3], temperature: f64) -> Self {
        HydrothermalVent {
            position,
            temperature,
            exit_velocity: 0.05,
            radius: 0.5,
            salinity: 0.005,
            mineral_loading: 0.0001,
            is_black_smoker: false,
        }
    }

    /// Mass flow rate (kg s⁻¹) from this vent.
    pub fn mass_flow_rate(&self) -> f64 {
        let area = PI * self.radius * self.radius;
        let rho = water_density(self.temperature, P_ATM + 200.0e5);
        rho * self.exit_velocity * area
    }

    /// Heat flux (W) from this vent.
    pub fn heat_flux(&self) -> f64 {
        let mdot = self.mass_flow_rate();
        let cp = water_cp(self.temperature);
        let dt = self.temperature - 277.15; // relative to 4°C ocean bottom
        mdot * cp * dt
    }
}

// ---------------------------------------------------------------------------
// Mineral precipitation/dissolution
// ---------------------------------------------------------------------------

/// Mineral type in hydrothermal system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mineral {
    /// Quartz / amorphous silica.
    Silica,
    /// Calcite (CaCO₃).
    Calcite,
    /// Pyrite (FeS₂).
    Pyrite,
    /// Anhydrite (CaSO₄).
    Anhydrite,
    /// Halite (NaCl).
    Halite,
}

impl Mineral {
    /// Molar mass (kg mol⁻¹) of the mineral.
    pub fn molar_mass(self) -> f64 {
        match self {
            Mineral::Silica => 0.06008,
            Mineral::Calcite => 0.10009,
            Mineral::Pyrite => 0.11987,
            Mineral::Anhydrite => 0.13614,
            Mineral::Halite => 0.05844,
        }
    }

    /// Mineral density (kg m⁻³).
    pub fn density(self) -> f64 {
        match self {
            Mineral::Silica => 2200.0,
            Mineral::Calcite => 2710.0,
            Mineral::Pyrite => 5010.0,
            Mineral::Anhydrite => 2960.0,
            Mineral::Halite => 2165.0,
        }
    }
}

/// Solubility of silica (kg m⁻³) in water at temperature `t` (K).
///
/// Uses empirical Van Lier fit.
pub fn silica_solubility(t: f64) -> f64 {
    let tc = t - 273.15;
    // mg/L to kg/m³ (1 mg/L = 1e-3 kg/m³)
    let s_mg_l = 7.695 * (tc * 0.053).exp();
    s_mg_l * 1.0e-3
}

/// Calcite solubility (mol m⁻³) as function of T (K) and P (Pa).
///
/// Simplified pressure and temperature dependence.
pub fn calcite_solubility(t: f64, p: f64) -> f64 {
    let t_ref = 298.15;
    let p_ref = P_ATM;
    let s_ref = 0.013; // mol/L ~ 13 mol/m³
    let dt = t - t_ref;
    let dp = p - p_ref;
    s_ref * (1.0 - 0.002 * dt + 1.0e-8 * dp) * 1000.0
}

/// Precipitation/dissolution rate (kg m⁻³ s⁻¹) using first-order kinetics.
///
/// Positive = precipitation, negative = dissolution.
pub fn mineral_reaction_rate(
    concentration: f64,
    solubility: f64,
    rate_constant: f64,
    surface_area: f64,
) -> f64 {
    rate_constant * surface_area * (concentration - solubility)
}

// ---------------------------------------------------------------------------
// Heat pipe effect
// ---------------------------------------------------------------------------

/// Heat pipe effect in geothermal vadose zone.
///
/// Models the two-phase counterflow of steam upward and liquid water downward
/// that greatly enhances effective thermal conductivity.
#[derive(Debug, Clone)]
pub struct HeatPipeModel {
    /// Capillary pressure (Pa) driving reflux.
    pub capillary_pressure: f64,
    /// Permeability of medium (m²).
    pub permeability: f64,
    /// Porosity.
    pub porosity: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl HeatPipeModel {
    /// Create a new heat pipe model.
    pub fn new(
        capillary_pressure: f64,
        permeability: f64,
        porosity: f64,
        temperature: f64,
    ) -> Self {
        HeatPipeModel {
            capillary_pressure,
            permeability,
            porosity,
            temperature,
        }
    }

    /// Steam flux (kg m⁻² s⁻¹) upward.
    pub fn steam_flux(&self) -> f64 {
        let mu_steam = 1.5e-5; // Pa·s
        let rho_steam = water_density(self.temperature, P_ATM);
        rho_steam * self.permeability * self.capillary_pressure / (mu_steam * 0.1)
    }

    /// Effective thermal conductivity (W m⁻¹ K⁻¹) of heat pipe.
    pub fn effective_conductivity(&self) -> f64 {
        let flux = self.steam_flux();
        let k_background = K_ROCK * self.porosity;
        k_background + L_VAP * flux / 100.0
    }
}

// ---------------------------------------------------------------------------
// EGS (Enhanced Geothermal System)
// ---------------------------------------------------------------------------

/// Enhanced Geothermal System extraction model.
#[derive(Debug, Clone)]
pub struct EgsSystem {
    /// Injection well position (m).
    pub injection_pos: [f64; 3],
    /// Production well position (m).
    pub production_pos: [f64; 3],
    /// Injection temperature (K).
    pub injection_temp: f64,
    /// Injection pressure (Pa).
    pub injection_pressure: f64,
    /// Production pressure (Pa).
    pub production_pressure: f64,
    /// Injection flow rate (m³ s⁻¹).
    pub flow_rate: f64,
    /// Reservoir rock type.
    pub rock_type: RockType,
    /// Reservoir depth (m).
    pub depth: f64,
    /// Reservoir temperature (K).
    pub reservoir_temp: f64,
}

impl EgsSystem {
    /// Create a typical EGS configuration.
    pub fn new(depth: f64, rock: RockType) -> Self {
        let thermal_grad = CrustalThermalGradient::continental();
        let t_reservoir = thermal_grad.temperature_at_depth(depth);
        let p_inj = P_ATM + 2700.0 * G * depth + 5.0e6;
        EgsSystem {
            injection_pos: [0.0, -depth, 0.0],
            production_pos: [300.0, -depth, 0.0],
            injection_temp: 303.15,
            injection_pressure: p_inj,
            production_pressure: p_inj - 2.0e6,
            flow_rate: 0.05,
            rock_type: rock,
            depth,
            reservoir_temp: t_reservoir,
        }
    }

    /// Thermal power extraction (W).
    pub fn thermal_power(&self) -> f64 {
        let rho = water_density(self.reservoir_temp, self.injection_pressure);
        let mdot = rho * self.flow_rate;
        let cp = water_cp(self.reservoir_temp);
        let dt = self.reservoir_temp - self.injection_temp;
        mdot * cp * dt
    }

    /// Effective thermal efficiency (dimensionless).
    pub fn thermal_efficiency(&self) -> f64 {
        let t_high = self.reservoir_temp;
        let t_low = self.injection_temp;
        1.0 - t_low / t_high // Carnot limit
    }

    /// Estimated reservoir impedance (Pa s m⁻³).
    pub fn reservoir_impedance(&self) -> f64 {
        let k = self.rock_type.permeability();
        let mu = water_viscosity(self.reservoir_temp);
        let l = norm3(sub3(self.production_pos, self.injection_pos));
        mu * l / (k * 1.0) // unit cross-section approximation
    }
}

// ---------------------------------------------------------------------------
// Hydrothermal alteration
// ---------------------------------------------------------------------------

/// Hydrothermal alteration reaction tracking.
///
/// Tracks the conversion of primary rock minerals to secondary alteration
/// products (e.g., feldspar to clay minerals).
#[derive(Debug, Clone)]
pub struct AlterationTracker {
    /// Volume fraction of unaltered primary mineral.
    pub primary_fraction: f64,
    /// Volume fraction of secondary alteration product.
    pub secondary_fraction: f64,
    /// Reaction rate constant (s⁻¹).
    pub rate_constant: f64,
    /// Activation energy (J mol⁻¹).
    pub activation_energy: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl AlterationTracker {
    /// Create a feldspar-to-kaolinite alteration tracker.
    pub fn feldspar_kaolinite(temperature: f64) -> Self {
        AlterationTracker {
            primary_fraction: 1.0,
            secondary_fraction: 0.0,
            rate_constant: 1.0e-11,
            activation_energy: 65_000.0,
            temperature,
        }
    }

    /// Arrhenius reaction rate (s⁻¹) at current temperature.
    pub fn arrhenius_rate(&self) -> f64 {
        self.rate_constant * ((-self.activation_energy) / (R_GAS * self.temperature)).exp()
    }

    /// Advance the alteration reaction by time step `dt` (s).
    pub fn step(&mut self, dt: f64) {
        let rate = self.arrhenius_rate();
        let dprimary = rate * self.primary_fraction * dt;
        self.primary_fraction = (self.primary_fraction - dprimary).max(0.0);
        self.secondary_fraction = 1.0 - self.primary_fraction;
    }
}

// ---------------------------------------------------------------------------
// Phase separation: steam-water separator
// ---------------------------------------------------------------------------

/// Steam-water separator / flash vessel model.
#[derive(Debug, Clone)]
pub struct FlashSeparator {
    /// Feed temperature (K).
    pub feed_temp: f64,
    /// Feed pressure (Pa).
    pub feed_pressure: f64,
    /// Flash pressure (Pa) — lower pressure into which fluid is flashed.
    pub flash_pressure: f64,
    /// Feed mass flow rate (kg s⁻¹).
    pub feed_flow: f64,
}

impl FlashSeparator {
    /// Create a flash separator for given feed conditions.
    pub fn new(feed_temp: f64, feed_pressure: f64, flash_pressure: f64, feed_flow: f64) -> Self {
        FlashSeparator {
            feed_temp,
            feed_pressure,
            flash_pressure,
            feed_flow,
        }
    }

    /// Steam mass flow rate produced (kg s⁻¹).
    pub fn steam_flow(&self) -> f64 {
        let x = self.flash_quality();
        x * self.feed_flow
    }

    /// Brine (liquid) mass flow rate (kg s⁻¹).
    pub fn brine_flow(&self) -> f64 {
        (1.0 - self.flash_quality()) * self.feed_flow
    }

    /// Flash steam quality.
    pub fn flash_quality(&self) -> f64 {
        let t_sat = boiling_point(self.flash_pressure);
        let dt = (self.feed_temp - t_sat) / 50.0;
        dt.clamp(0.0, 1.0)
    }

    /// Power output if steam drives a turbine with isentropic efficiency `eta`.
    pub fn turbine_power(&self, eta: f64) -> f64 {
        let mdot_steam = self.steam_flow();
        let h_in = water_cp(self.feed_temp) * self.feed_temp;
        let t_out = boiling_point(P_ATM);
        let h_out = water_cp(t_out) * t_out;
        eta * mdot_steam * (h_in - h_out).max(0.0)
    }
}

// ---------------------------------------------------------------------------
// SPH pressure gradient and viscosity terms
// ---------------------------------------------------------------------------

/// Symmetric SPH pressure gradient contribution.
///
/// Returns acceleration contribution to particle `i` from neighbor `j`.
pub fn sph_pressure_accel(
    pos_i: [f64; 3],
    pos_j: [f64; 3],
    mass_j: f64,
    rho_i: f64,
    rho_j: f64,
    p_i: f64,
    p_j: f64,
    h: f64,
) -> [f64; 3] {
    let r_vec = sub3(pos_i, pos_j);
    let grad_w = cubic_kernel_grad(r_vec, h);
    let coeff = mass_j * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j));
    scale3(grad_w, -coeff)
}

/// Fluid-state parameters for the artificial viscosity pair computation.
#[derive(Debug, Clone, Copy)]
pub struct ArtViscPairParams {
    /// Position of particle i \[m\]
    pub pos_i: [f64; 3],
    /// Position of particle j \[m\]
    pub pos_j: [f64; 3],
    /// Velocity of particle i \[m/s\]
    pub vel_i: [f64; 3],
    /// Velocity of particle j \[m/s\]
    pub vel_j: [f64; 3],
    /// Density of particle i \[kg/m³\]
    pub rho_i: f64,
    /// Density of particle j \[kg/m³\]
    pub rho_j: f64,
    /// Sound speed of particle i \[m/s\]
    pub c_i: f64,
    /// Sound speed of particle j \[m/s\]
    pub c_j: f64,
}

/// SPH artificial viscosity acceleration (Monaghan 1992).
pub fn sph_artificial_viscosity(
    pair: ArtViscPairParams,
    mass_j: f64,
    h: f64,
    alpha: f64,
    beta: f64,
) -> [f64; 3] {
    let ArtViscPairParams {
        pos_i,
        pos_j,
        vel_i,
        vel_j,
        rho_i,
        rho_j,
        c_i,
        c_j,
    } = pair;
    let r_vec = sub3(pos_i, pos_j);
    let v_vec = sub3(vel_i, vel_j);
    let rv = dot3(r_vec, v_vec);
    if rv >= 0.0 {
        return [0.0; 3];
    }
    let r2 = dot3(r_vec, r_vec);
    let eta2 = 0.01 * h * h;
    let mu = h * rv / (r2 + eta2);
    let rho_ij = 0.5 * (rho_i + rho_j);
    let c_ij = 0.5 * (c_i + c_j);
    let pi_ij = (-alpha * c_ij * mu + beta * mu * mu) / rho_ij;
    let grad_w = cubic_kernel_grad(r_vec, h);
    scale3(grad_w, -mass_j * pi_ij)
}

/// SPH thermal diffusion (Brookshaw 1985 approximation).
pub fn sph_thermal_diffusion(
    pos_i: [f64; 3],
    pos_j: [f64; 3],
    mass_j: f64,
    rho_j: f64,
    t_i: f64,
    t_j: f64,
    kappa_i: f64,
    kappa_j: f64,
    cp_i: f64,
    h: f64,
) -> f64 {
    let r_vec = sub3(pos_i, pos_j);
    let r = norm3(r_vec);
    if r < 1.0e-12 {
        return 0.0;
    }
    let kappa_ij = 2.0 * kappa_i * kappa_j / (kappa_i + kappa_j);
    let dw_dr = {
        let q = r / h;
        let alpha = 1.0 / (PI * h * h * h);
        if q < 1.0 {
            alpha * (-3.0 * q + 2.25 * q * q) / h
        } else if q < 2.0 {
            let t = 2.0 - q;
            alpha * (-0.75 * t * t) / h
        } else {
            0.0
        }
    };
    let e_ij = dot3(normalize3(r_vec), r_vec) / (r + 1.0e-12);
    (mass_j / rho_j) * kappa_ij * (t_i - t_j) * dw_dr * e_ij / (cp_i * r)
}

// ---------------------------------------------------------------------------
// SPH density summation
// ---------------------------------------------------------------------------

/// SPH density estimate from neighbor particles.
///
/// Sums kernel contributions from all neighbors.
pub fn sph_density_sum(
    pos_i: [f64; 3],
    neighbors: &[(f64, [f64; 3])], // (mass_j, pos_j)
    h: f64,
) -> f64 {
    neighbors
        .iter()
        .map(|(m_j, pos_j)| {
            let r = norm3(sub3(pos_i, *pos_j));
            m_j * cubic_kernel(r, h)
        })
        .sum()
}

// ---------------------------------------------------------------------------
// Geothermal simulation
// ---------------------------------------------------------------------------

/// Configuration for a geothermal SPH simulation.
#[derive(Debug, Clone)]
pub struct GeothermalConfig {
    /// Smoothing length (m).
    pub smoothing_length: f64,
    /// Time step (s).
    pub dt: f64,
    /// Speed of sound reference (m s⁻¹).
    pub c_ref: f64,
    /// Artificial viscosity alpha parameter.
    pub alpha_visc: f64,
    /// Thermal diffusivity (m² s⁻¹).
    pub thermal_diffusivity: f64,
    /// Rock type for porous medium.
    pub rock_type: RockType,
    /// Geothermal gradient (K m⁻¹).
    pub geothermal_gradient: f64,
    /// Enable mineral precipitation.
    pub enable_chemistry: bool,
    /// Number of SPH steps between chemistry updates.
    pub chemistry_interval: usize,
}

impl Default for GeothermalConfig {
    fn default() -> Self {
        GeothermalConfig {
            smoothing_length: 0.05,
            dt: 1.0e-3,
            c_ref: 1500.0,
            alpha_visc: 0.1,
            thermal_diffusivity: 1.4e-7,
            rock_type: RockType::Granite,
            geothermal_gradient: GRAD_T_CRUST,
            enable_chemistry: true,
            chemistry_interval: 10,
        }
    }
}

/// Geothermal SPH simulation state.
#[derive(Debug, Clone)]
pub struct GeothermalSimulation {
    /// Configuration parameters.
    pub config: GeothermalConfig,
    /// All fluid particles.
    pub particles: Vec<GeoParticle>,
    /// Hydrothermal vents acting as sources.
    pub vents: Vec<HydrothermalVent>,
    /// EGS system (optional).
    pub egs: Option<EgsSystem>,
    /// Current simulation time (s).
    pub time: f64,
    /// Current step counter.
    pub step: usize,
    /// Alteration trackers for each particle.
    pub alteration: Vec<AlterationTracker>,
}

impl GeothermalSimulation {
    /// Create a new geothermal simulation.
    pub fn new(config: GeothermalConfig) -> Self {
        GeothermalSimulation {
            config,
            particles: Vec::new(),
            vents: Vec::new(),
            egs: None,
            time: 0.0,
            step: 0,
            alteration: Vec::new(),
        }
    }

    /// Add a fluid particle to the simulation.
    pub fn add_particle(&mut self, p: GeoParticle) {
        self.alteration
            .push(AlterationTracker::feldspar_kaolinite(p.temperature));
        self.particles.push(p);
    }

    /// Add a hydrothermal vent.
    pub fn add_vent(&mut self, v: HydrothermalVent) {
        self.vents.push(v);
    }

    /// Set an EGS system.
    pub fn set_egs(&mut self, egs: EgsSystem) {
        self.egs = Some(egs);
    }

    /// Total fluid kinetic energy (J).
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Total thermal energy (J).
    pub fn total_thermal_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.thermal_energy()).sum()
    }

    /// Mean temperature of all particles (K).
    pub fn mean_temperature(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.particles.iter().map(|p| p.temperature).sum();
        sum / self.particles.len() as f64
    }

    /// Advance the simulation by one time step.
    pub fn step_forward(&mut self) {
        let dt = self.config.dt;
        let h = self.config.smoothing_length;
        let n = self.particles.len();

        // 1. Density summation
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.pos).collect();
        let masses: Vec<f64> = self.particles.iter().map(|p| p.mass).collect();
        for i in 0..n {
            let neighbors: Vec<(f64, [f64; 3])> = (0..n)
                .filter(|&j| norm3(sub3(positions[i], positions[j])) < 2.0 * h)
                .map(|j| (masses[j], positions[j]))
                .collect();
            self.particles[i].density = sph_density_sum(positions[i], &neighbors, h);
        }

        // 2. Pressure from Tait EOS
        let rho0 = 1000.0;
        let c = self.config.c_ref;
        let gamma = 7.0;
        for p in &mut self.particles {
            p.pressure = rho0 * c * c / gamma * ((p.density / rho0).powf(gamma) - 1.0);
        }

        // 3. Accumulate accelerations
        let pressures: Vec<f64> = self.particles.iter().map(|p| p.pressure).collect();
        let densities: Vec<f64> = self.particles.iter().map(|p| p.density).collect();
        let velocities: Vec<[f64; 3]> = self.particles.iter().map(|p| p.vel).collect();

        for i in 0..n {
            let mut acc = [0.0f64; 3];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = norm3(sub3(positions[i], positions[j]));
                if r >= 2.0 * h {
                    continue;
                }
                let pa = sph_pressure_accel(
                    positions[i],
                    positions[j],
                    masses[j],
                    densities[i],
                    densities[j],
                    pressures[i],
                    pressures[j],
                    h,
                );
                let visc = sph_artificial_viscosity(
                    ArtViscPairParams {
                        pos_i: positions[i],
                        pos_j: positions[j],
                        vel_i: velocities[i],
                        vel_j: velocities[j],
                        rho_i: densities[i],
                        rho_j: densities[j],
                        c_i: c,
                        c_j: c,
                    },
                    masses[j],
                    h,
                    self.config.alpha_visc,
                    0.0,
                );
                acc = add3(acc, add3(pa, visc));
            }
            // Gravity
            acc[1] -= G;
            // Buoyancy correction
            let beta_t = 2.5e-4;
            let dt_temp = self.particles[i].temperature - 277.15;
            let buoy = buoyancy_force(densities[i], beta_t, dt_temp);
            acc = add3(acc, scale3(buoy, 1.0 / densities[i].max(1.0)));
            self.particles[i].acc = acc;
        }

        // 4. Integrate velocity and position (Symplectic Euler)
        for p in &mut self.particles {
            p.vel = add3(p.vel, scale3(p.acc, dt));
            p.pos = add3(p.pos, scale3(p.vel, dt));
        }

        // 5. Thermal diffusion
        let temps: Vec<f64> = self.particles.iter().map(|p| p.temperature).collect();
        for i in 0..n {
            let mut d_temp = 0.0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let r = norm3(sub3(positions[i], positions[j]));
                if r >= 2.0 * h {
                    continue;
                }
                let kappa_i = water_thermal_conductivity(temps[i]);
                let kappa_j = water_thermal_conductivity(temps[j]);
                let cp_i = water_cp(temps[i]);
                d_temp += sph_thermal_diffusion(
                    positions[i],
                    positions[j],
                    masses[j],
                    densities[j],
                    temps[i],
                    temps[j],
                    kappa_i,
                    kappa_j,
                    cp_i,
                    h,
                );
            }
            self.particles[i].temperature += d_temp * dt;
        }

        // 6. Update phase and alteration chemistry
        for p in &mut self.particles {
            p.update_phase();
        }

        if self.config.enable_chemistry && self.step.is_multiple_of(self.config.chemistry_interval)
        {
            let chem_dt = dt * self.config.chemistry_interval as f64;
            for alt in &mut self.alteration {
                alt.step(chem_dt);
            }
        }

        self.time += dt;
        self.step += 1;
    }
}

// ---------------------------------------------------------------------------
// Nusselt number correlations
// ---------------------------------------------------------------------------

/// Nusselt number for natural convection in porous media (Horton-Rogers-Lapwood).
///
/// `ra` is the Rayleigh-Darcy number.
pub fn nusselt_porous(ra: f64) -> f64 {
    // Combarnous & Bories (1975) correlation for porous-medium natural convection.
    // Nu = 1 for Ra < Ra_c ≈ 40; Nu = 0.069 * Ra^0.72 for Ra ≥ 40.
    if ra < 40.0 {
        1.0
    } else {
        0.069 * ra.powf(0.72)
    }
}

/// Rayleigh-Darcy number for hydrothermal convection.
///
/// `k` permeability (m²), `h_depth` layer thickness (m), `beta` thermal expansion (K⁻¹),
/// `delta_t` temperature difference (K), `kappa` thermal diffusivity (m² s⁻¹), `mu` viscosity (Pa·s).
pub fn rayleigh_darcy(k: f64, h_depth: f64, beta: f64, delta_t: f64, kappa: f64, mu: f64) -> f64 {
    let rho = 1000.0;
    rho * G * beta * delta_t * k * h_depth / (mu * kappa)
}

// ---------------------------------------------------------------------------
// Integral energy balance
// ---------------------------------------------------------------------------

/// Steady-state geothermal upflow temperature profile.
///
/// Solves 1-D advection-diffusion for fluid temperature rising from depth.
pub fn geothermal_temperature_profile(
    depth_max: f64,
    n_points: usize,
    flow_velocity: f64,
    kappa: f64,
    grad_t: f64,
) -> Vec<(f64, f64)> {
    let dz = depth_max / (n_points as f64 - 1.0);
    let t_surface = 288.15;
    let mut profile = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let z = i as f64 * dz;
        // Exponential boundary layer correction
        let t_cond = t_surface + grad_t * z;
        let pe = flow_velocity * z / kappa;
        let t_adv = t_surface + grad_t * z * (1.0 - (-pe).exp());
        let t = t_cond * 0.5 + t_adv * 0.5;
        profile.push((z, t));
    }
    profile
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boiling_point_atmospheric() {
        let t = boiling_point(P_ATM);
        // Should be close to 373.15 K (100°C)
        assert!((t - 373.15).abs() < 5.0, "Boiling point at 1 atm: {t}");
    }

    #[test]
    fn test_water_density_liquid() {
        let rho = water_density(300.0, P_ATM);
        // Should be around 995–1000 kg/m³
        assert!(rho > 990.0 && rho < 1010.0, "Liquid water density: {rho}");
    }

    #[test]
    fn test_water_density_steam() {
        // At 200°C and 1 atm, steam should be present
        let rho = water_density(473.15, P_ATM);
        // steam density ~ 0.6–0.8 kg/m³ at 1 atm, 200°C
        assert!(rho < 2.0, "Steam density should be low: {rho}");
    }

    #[test]
    fn test_water_density_supercritical() {
        let rho = water_density(T_CRIT_WATER + 10.0, P_CRIT_WATER + 1.0e6);
        assert!(rho > 50.0 && rho < 500.0, "Supercritical density: {rho}");
    }

    #[test]
    fn test_fluid_phase_liquid() {
        let phase = FluidPhase::from_tp(300.0, P_ATM);
        assert_eq!(phase, FluidPhase::LiquidWater);
    }

    #[test]
    fn test_fluid_phase_steam() {
        let phase = FluidPhase::from_tp(500.0, P_ATM);
        assert_eq!(phase, FluidPhase::Steam);
    }

    #[test]
    fn test_fluid_phase_supercritical() {
        let phase = FluidPhase::from_tp(660.0, 23.0e6);
        assert_eq!(phase, FluidPhase::Supercritical);
    }

    #[test]
    fn test_supercritical_fluid_new() {
        let sc = SupercriticalFluid::new(660.0, 23.0e6);
        assert!(sc.is_supercritical());
        assert!(sc.density > 0.0);
        assert!(sc.viscosity > 0.0);
        assert!(sc.cp > 0.0);
    }

    #[test]
    fn test_supercritical_reduced_properties() {
        let sc = SupercriticalFluid::new(700.0, 25.0e6);
        let tau = sc.reduced_temperature();
        assert!((tau - T_CRIT_WATER / 700.0).abs() < 1.0e-10);
        let pi = sc.reduced_pressure();
        assert!((pi - 25.0e6 / P_CRIT_WATER).abs() < 1.0e-6);
    }

    #[test]
    fn test_rock_type_permeability_ordering() {
        assert!(RockType::FracturedRock.permeability() > RockType::Sandstone.permeability());
        assert!(RockType::Sandstone.permeability() > RockType::Granite.permeability());
    }

    #[test]
    fn test_darcy_velocity() {
        let v = darcy_velocity(1.0e-12, 1.0e-3, [1000.0, 0.0, 0.0]);
        // Should be negative (opposite to grad_p)
        assert!(v[0] < 0.0);
        assert!(v[1].abs() < 1.0e-30);
    }

    #[test]
    fn test_buoyancy_force_hot() {
        let f = buoyancy_force(1000.0, 2.5e-4, 100.0);
        // Upward buoyancy for hot fluid (f[1] > 0)
        assert!(f[1] > 0.0, "Buoyancy should be positive for hot fluid");
    }

    #[test]
    fn test_geo_particle_creation() {
        let p = GeoParticle::new([0.0; 3], 1.0, 400.0, P_ATM);
        assert_eq!(p.phase, FluidPhase::Steam);
        assert!(p.density > 0.0);
        assert!(p.enthalpy > 0.0);
    }

    #[test]
    fn test_crustal_gradient_temperature() {
        let crust = CrustalThermalGradient::continental();
        let t5 = crust.temperature_at_depth(5000.0);
        assert!((t5 - (crust.t_surface + 0.03 * 5000.0)).abs() < 1.0e-9);
    }

    #[test]
    fn test_volcanic_gradient_higher() {
        let cont = CrustalThermalGradient::continental();
        let volc = CrustalThermalGradient::volcanic();
        assert!(volc.gradient > cont.gradient);
    }

    #[test]
    fn test_hydrothermal_vent_mass_flow() {
        let vent = HydrothermalVent::black_smoker([0.0; 3]);
        let mdot = vent.mass_flow_rate();
        assert!(mdot > 0.0);
    }

    #[test]
    fn test_hydrothermal_vent_heat_flux() {
        let vent = HydrothermalVent::hot_spring([0.0; 3], 350.0);
        let q = vent.heat_flux();
        assert!(q > 0.0);
    }

    #[test]
    fn test_silica_solubility_increases_with_temp() {
        let s_cold = silica_solubility(300.0);
        let s_hot = silica_solubility(500.0);
        assert!(s_hot > s_cold, "Silica more soluble at higher T");
    }

    #[test]
    fn test_mineral_reaction_rate_precipitation() {
        let rate = mineral_reaction_rate(0.01, 0.005, 1.0e-6, 10.0);
        assert!(rate > 0.0, "Oversaturated: should precipitate");
    }

    #[test]
    fn test_mineral_reaction_rate_dissolution() {
        let rate = mineral_reaction_rate(0.001, 0.01, 1.0e-6, 10.0);
        assert!(rate < 0.0, "Undersaturated: should dissolve");
    }

    #[test]
    fn test_heat_pipe_effective_conductivity() {
        let hp = HeatPipeModel::new(1000.0, 1.0e-13, 0.1, 373.15);
        let k_eff = hp.effective_conductivity();
        assert!(k_eff > 0.0);
    }

    #[test]
    fn test_egs_thermal_power() {
        let egs = EgsSystem::new(3000.0, RockType::Granite);
        let p = egs.thermal_power();
        assert!(p > 1.0e6, "EGS should produce at least 1 MW");
    }

    #[test]
    fn test_egs_carnot_efficiency() {
        let egs = EgsSystem::new(3000.0, RockType::Granite);
        let eta = egs.thermal_efficiency();
        assert!(eta > 0.0 && eta < 1.0);
    }

    #[test]
    fn test_alteration_tracker_step() {
        let mut alt = AlterationTracker::feldspar_kaolinite(350.0);
        let initial = alt.primary_fraction;
        alt.step(1.0e6); // 11.6 days
        assert!(alt.primary_fraction < initial, "Alteration should proceed");
        assert!((alt.primary_fraction + alt.secondary_fraction - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_flash_separator_quality() {
        let flash = FlashSeparator::new(473.15, 10.0e5, P_ATM, 1.0);
        let x = flash.flash_quality();
        assert!((0.0..=1.0).contains(&x));
    }

    #[test]
    fn test_flash_separator_flow_balance() {
        let flash = FlashSeparator::new(473.15, 10.0e5, P_ATM, 10.0);
        let total = flash.steam_flow() + flash.brine_flow();
        assert!((total - 10.0).abs() < 1.0e-10, "Mass balance: {total}");
    }

    #[test]
    fn test_cubic_kernel_unity() {
        // Test kernel is positive
        let w = cubic_kernel(0.0, 0.1);
        assert!(w > 0.0);
    }

    #[test]
    fn test_cubic_kernel_zero_outside() {
        let w = cubic_kernel(0.21, 0.1);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_sph_density_sum() {
        let pos_i = [0.0; 3];
        let neighbors = vec![(1.0, [0.01, 0.0, 0.0]), (1.0, [-0.01, 0.0, 0.0])];
        let rho = sph_density_sum(pos_i, &neighbors, 0.05);
        assert!(rho > 0.0);
    }

    #[test]
    fn test_nusselt_porous_conductive() {
        let nu = nusselt_porous(10.0);
        assert!((nu - 1.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_nusselt_porous_convective() {
        let nu = nusselt_porous(1000.0);
        assert!(nu > 1.0);
    }

    #[test]
    fn test_rayleigh_darcy_positive() {
        let ra = rayleigh_darcy(1.0e-13, 100.0, 2.5e-4, 200.0, 1.4e-7, 1.0e-4);
        assert!(ra > 0.0);
    }

    #[test]
    fn test_temperature_profile_length() {
        let profile = geothermal_temperature_profile(5000.0, 50, 0.001, 1.4e-7, 0.03);
        assert_eq!(profile.len(), 50);
    }

    #[test]
    fn test_temperature_profile_monotonic() {
        let profile = geothermal_temperature_profile(1000.0, 10, 0.0001, 1.4e-7, 0.03);
        for i in 1..profile.len() {
            assert!(
                profile[i].1 >= profile[i - 1].1,
                "Temperature should increase with depth"
            );
        }
    }

    #[test]
    fn test_simulation_step() {
        let config = GeothermalConfig {
            smoothing_length: 0.2,
            dt: 1.0e-3,
            ..Default::default()
        };
        let mut sim = GeothermalSimulation::new(config);
        for i in 0..3 {
            let p = GeoParticle::new([i as f64 * 0.1, 0.0, 0.0], 1.0, 350.0, P_ATM);
            sim.add_particle(p);
        }
        sim.step_forward();
        assert_eq!(sim.step, 1);
        assert!(sim.time > 0.0);
    }

    #[test]
    fn test_simulation_energy_positive() {
        let mut sim = GeothermalSimulation::new(GeothermalConfig::default());
        let p = GeoParticle::new([0.0; 3], 1.0, 400.0, P_ATM);
        sim.add_particle(p);
        let ke = sim.total_kinetic_energy();
        let te = sim.total_thermal_energy();
        assert!(ke >= 0.0);
        assert!(te > 0.0);
    }

    #[test]
    fn test_water_viscosity_decreases_with_temp() {
        let mu_cold = water_viscosity(280.0);
        let mu_hot = water_viscosity(360.0);
        assert!(
            mu_cold > mu_hot,
            "Viscosity should decrease with temperature"
        );
    }

    #[test]
    fn test_water_cp_reasonable() {
        let cp = water_cp(300.0);
        // Specific heat of water ~4180 J/(kg K)
        assert!(cp > 4000.0 && cp < 4500.0, "Water Cp: {cp}");
    }
}
