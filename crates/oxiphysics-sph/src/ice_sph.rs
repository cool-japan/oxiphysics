#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ice/glaciology SPH module for ice sheet dynamics, calving, and rheology.
//!
//! This module provides smoothed-particle hydrodynamics methods for simulating:
//! - Ice sheet flow using the Shallow Ice Approximation (SIA)
//! - Glen's flow law (power-law rheology)
//! - Calving via crevasse-depth and Von Mises stress criteria
//! - Thermal evolution in ice (strain heating, geothermal flux)
//! - Crystal orientation fabric and enhancement factors
//! - Glacial erosion (abrasion, plucking, sediment transport)
//! - Ice volume analysis and mass balance tracking

use rand::RngExt;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Ice density in kg/m^3.
const ICE_DENSITY: f64 = 917.0;

/// Gravitational acceleration in m/s^2.
const GRAVITY: f64 = 9.81;

/// Universal gas constant in J/(mol K).
const GAS_CONSTANT: f64 = 8.314;

/// Glen flow law exponent (dimensionless).
const GLEN_N: f64 = 3.0;

/// Thermal conductivity of ice in W/(m K).
const ICE_THERMAL_CONDUCTIVITY: f64 = 2.1;

/// Specific heat capacity of ice in J/(kg K).
const ICE_SPECIFIC_HEAT: f64 = 2090.0;

/// Latent heat of fusion for ice in J/kg.
#[allow(dead_code)]
const ICE_LATENT_HEAT: f64 = 334_000.0;

/// Melting point of ice at 1 atm in K.
const ICE_MELTING_POINT: f64 = 273.15;

/// Clausius-Clapeyron slope in K/Pa.
const CLAUSIUS_CLAPEYRON: f64 = 7.42e-8;

/// Seconds per year.
const SECONDS_PER_YEAR: f64 = 365.25 * 24.0 * 3600.0;

// ---------------------------------------------------------------------------
// IceParticle
// ---------------------------------------------------------------------------

/// A single SPH particle representing a parcel of ice.
#[derive(Debug, Clone)]
pub struct IceParticle {
    /// Position \[x, y, z\] in metres.
    pub pos: [f64; 3],
    /// Velocity \[vx, vy, vz\] in m/s.
    pub vel: [f64; 3],
    /// Deviatoric stress tensor (symmetric, stored as 6 components: xx yy zz xy xz yz).
    pub stress: [f64; 6],
    /// Temperature in Kelvin.
    pub temperature: f64,
    /// Effective viscosity in Pa s.
    pub viscosity: f64,
    /// Particle age in years.
    pub age: f64,
    /// Particle mass in kg.
    pub mass: f64,
    /// Smoothing length in metres.
    pub smoothing_length: f64,
    /// Density in kg/m^3.
    pub density: f64,
    /// Pressure in Pa.
    pub pressure: f64,
    /// Damage parameter \[0, 1\].
    pub damage: f64,
    /// Water content fraction \[0, 1\].
    pub water_content: f64,
}

impl IceParticle {
    /// Create a new ice particle at a given position and temperature.
    pub fn new(pos: [f64; 3], temperature: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            stress: [0.0; 6],
            temperature,
            viscosity: 1.0e13,
            age: 0.0,
            mass: ICE_DENSITY * 1.0, // unit volume default
            smoothing_length: 1.0,
            density: ICE_DENSITY,
            pressure: 0.0,
            damage: 0.0,
            water_content: 0.0,
        }
    }

    /// Compute overburden pressure at this particle's depth.
    pub fn overburden_pressure(&self, surface_elevation: f64) -> f64 {
        let depth = (surface_elevation - self.pos[2]).max(0.0);
        ICE_DENSITY * GRAVITY * depth
    }

    /// Pressure-melting temperature accounting for overburden.
    pub fn pressure_melting_point(&self, surface_elevation: f64) -> f64 {
        let p = self.overburden_pressure(surface_elevation);
        ICE_MELTING_POINT - CLAUSIUS_CLAPEYRON * p
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
}

// ---------------------------------------------------------------------------
// GlenFlowLaw
// ---------------------------------------------------------------------------

/// Glen's flow law parameters for ice rheology.
///
/// The constitutive relation is: strain_rate = A(T) * tau^n
/// where n = 3 (Glen exponent), A is a temperature-dependent rate factor,
/// and tau is the effective stress.
#[derive(Debug, Clone)]
pub struct GlenFlowLaw {
    /// Flow law exponent (typically 3).
    pub n: f64,
    /// Pre-exponential constant in Pa^{-n} s^{-1}.
    pub a0: f64,
    /// Activation energy for creep in J/mol (below -10C).
    pub activation_energy_cold: f64,
    /// Activation energy for creep in J/mol (above -10C).
    pub activation_energy_warm: f64,
    /// Enhancement factor (anisotropy, impurities, etc.).
    pub enhancement_factor: f64,
}

impl Default for GlenFlowLaw {
    fn default() -> Self {
        Self {
            n: GLEN_N,
            a0: 3.985e-13, // Pa^-3 s^-1 for T < 263K
            activation_energy_cold: 60_000.0,
            activation_energy_warm: 139_000.0,
            enhancement_factor: 1.0,
        }
    }
}

impl GlenFlowLaw {
    /// Create a new Glen flow law with default parameters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create with custom enhancement factor.
    pub fn with_enhancement(enhancement: f64) -> Self {
        Self {
            enhancement_factor: enhancement,
            ..Default::default()
        }
    }

    /// Temperature-dependent rate factor A(T) in Pa^{-n} s^{-1}.
    ///
    /// Uses Arrhenius relationship with different activation energies
    /// above and below 263 K (-10 C).
    pub fn rate_factor(&self, temperature: f64) -> f64 {
        let t_star = 263.15; // -10 C
        let q = if temperature < t_star {
            self.activation_energy_cold
        } else {
            self.activation_energy_warm
        };
        let a0 = if temperature < t_star {
            self.a0
        } else {
            1.916e3 // Pa^-3 s^-1 for T >= 263K
        };
        self.enhancement_factor * a0 * (-q / (GAS_CONSTANT * temperature)).exp()
    }

    /// Effective strain rate from effective stress and temperature.
    ///
    /// strain_rate_eff = A(T) * tau_eff^n
    pub fn strain_rate(&self, effective_stress: f64, temperature: f64) -> f64 {
        let a = self.rate_factor(temperature);
        a * effective_stress.powf(self.n)
    }

    /// Effective viscosity eta = tau / (2 * strain_rate).
    ///
    /// Regularised to avoid division by zero.
    pub fn effective_viscosity(&self, effective_stress: f64, temperature: f64) -> f64 {
        let sr = self.strain_rate(effective_stress, temperature);
        let sr_reg = sr.max(1.0e-30);
        effective_stress / (2.0 * sr_reg)
    }

    /// Inverse: effective stress from effective strain rate.
    pub fn effective_stress(&self, strain_rate_eff: f64, temperature: f64) -> f64 {
        let a = self.rate_factor(temperature);
        if a < 1.0e-50 {
            return 0.0;
        }
        (strain_rate_eff / a).powf(1.0 / self.n)
    }

    /// Compute the deviatoric stress tensor from a strain rate tensor.
    ///
    /// tau_ij = 2 * eta * eps_ij  where eta is the effective viscosity.
    pub fn deviatoric_stress(&self, strain_rate_tensor: &[f64; 6], temperature: f64) -> [f64; 6] {
        let eps_eff = second_invariant_strain_rate(strain_rate_tensor);
        let sigma_eff = self.effective_stress(eps_eff, temperature);
        let eta = if eps_eff > 1.0e-30 {
            sigma_eff / (2.0 * eps_eff)
        } else {
            1.0e15
        };
        let mut tau = [0.0; 6];
        for i in 0..6 {
            tau[i] = 2.0 * eta * strain_rate_tensor[i];
        }
        tau
    }
}

/// Compute the second invariant of a symmetric strain-rate tensor.
///
/// eps_eff = sqrt(0.5 * eps_ij * eps_ij)
pub fn second_invariant_strain_rate(eps: &[f64; 6]) -> f64 {
    let sum = eps[0] * eps[0]
        + eps[1] * eps[1]
        + eps[2] * eps[2]
        + 2.0 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]);
    (0.5 * sum).sqrt()
}

/// Compute the second invariant of a symmetric stress tensor.
pub fn second_invariant_stress(tau: &[f64; 6]) -> f64 {
    second_invariant_strain_rate(tau)
}

// ---------------------------------------------------------------------------
// IceSheetFlow
// ---------------------------------------------------------------------------

/// Ice sheet flow model using the Shallow Ice Approximation (SIA).
///
/// The SIA relates the depth-averaged velocity to the ice thickness
/// and surface slope via Glen's flow law.
#[derive(Debug, Clone)]
pub struct IceSheetFlow {
    /// Glen flow law parameters.
    pub glen: GlenFlowLaw,
    /// Basal sliding coefficient in m Pa^{-1} yr^{-1}.
    pub sliding_coefficient: f64,
    /// Basal sliding exponent (typically 1 or 3).
    pub sliding_exponent: f64,
    /// Bed elevation profile (x, z_bed) pairs.
    pub bed_profile: Vec<[f64; 2]>,
    /// Surface elevation profile (x, z_surface) pairs.
    pub surface_profile: Vec<[f64; 2]>,
}

impl IceSheetFlow {
    /// Create a new ice sheet flow model.
    pub fn new(glen: GlenFlowLaw) -> Self {
        Self {
            glen,
            sliding_coefficient: 0.0,
            sliding_exponent: 1.0,
            bed_profile: Vec::new(),
            surface_profile: Vec::new(),
        }
    }

    /// Set bed and surface profiles from arrays.
    pub fn set_profiles(&mut self, bed: Vec<[f64; 2]>, surface: Vec<[f64; 2]>) {
        self.bed_profile = bed;
        self.surface_profile = surface;
    }

    /// Enable basal sliding with given coefficient and exponent.
    pub fn set_sliding(&mut self, coeff: f64, exponent: f64) {
        self.sliding_coefficient = coeff;
        self.sliding_exponent = exponent;
    }

    /// Driving stress at a given surface slope and ice thickness.
    ///
    /// tau_d = rho * g * H * |dS/dx|
    pub fn driving_stress(thickness: f64, surface_slope: f64) -> f64 {
        ICE_DENSITY * GRAVITY * thickness * surface_slope.abs()
    }

    /// SIA depth-averaged velocity (deformational component only).
    ///
    /// u_d = 2A / (n+1) * (rho g |dS/dx|)^n * H^{n+1}
    pub fn sia_deformation_velocity(
        &self,
        thickness: f64,
        surface_slope: f64,
        temperature: f64,
    ) -> f64 {
        let a = self.glen.rate_factor(temperature);
        let n = self.glen.n;
        let tau_base = ICE_DENSITY * GRAVITY * surface_slope.abs();
        2.0 * a / (n + 1.0) * tau_base.powf(n) * thickness.powf(n + 1.0)
    }

    /// Basal sliding velocity.
    ///
    /// u_b = C * tau_b^m  where tau_b is the basal shear stress.
    pub fn basal_sliding_velocity(&self, thickness: f64, surface_slope: f64) -> f64 {
        let tau_b = Self::driving_stress(thickness, surface_slope);
        self.sliding_coefficient * tau_b.powf(self.sliding_exponent)
    }

    /// Total surface velocity (deformation + sliding).
    pub fn surface_velocity(&self, thickness: f64, surface_slope: f64, temperature: f64) -> f64 {
        let u_d = self.sia_deformation_velocity(thickness, surface_slope, temperature);
        let u_b = self.basal_sliding_velocity(thickness, surface_slope);
        u_d + u_b
    }

    /// Velocity profile with depth (z from bed).
    ///
    /// Returns velocity at height z above bed.
    pub fn velocity_at_depth(
        &self,
        z: f64,
        thickness: f64,
        surface_slope: f64,
        temperature: f64,
    ) -> f64 {
        let a = self.glen.rate_factor(temperature);
        let n = self.glen.n;
        let tau_factor = (ICE_DENSITY * GRAVITY * surface_slope.abs()).powf(n);
        let u_b = self.basal_sliding_velocity(thickness, surface_slope);
        let u_def = 2.0 * a * tau_factor / (n + 1.0)
            * (thickness.powf(n + 1.0) - (thickness - z).max(0.0).powf(n + 1.0));
        u_b + u_def
    }

    /// Compute ice flux per unit width: q = integral_0^H u(z) dz.
    pub fn ice_flux(&self, thickness: f64, surface_slope: f64, temperature: f64) -> f64 {
        let a = self.glen.rate_factor(temperature);
        let n = self.glen.n;
        let tau_factor = (ICE_DENSITY * GRAVITY * surface_slope.abs()).powf(n);
        let q_def = 2.0 * a * tau_factor / (n + 2.0) * thickness.powf(n + 2.0);
        let q_slide = self.basal_sliding_velocity(thickness, surface_slope) * thickness;
        q_def + q_slide
    }

    /// SIA ice thickness evolution: dH/dt = M - d(q)/dx.
    ///
    /// Returns thickness rate of change at a point given upstream and downstream fluxes.
    pub fn thickness_change_rate(
        flux_upstream: f64,
        flux_downstream: f64,
        dx: f64,
        surface_mass_balance: f64,
    ) -> f64 {
        surface_mass_balance - (flux_downstream - flux_upstream) / dx
    }
}

// ---------------------------------------------------------------------------
// CalvingModel
// ---------------------------------------------------------------------------

/// Calving criterion enumeration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CalvingCriterion {
    /// Benn crevasse-depth model.
    CrevasseDepth,
    /// Von Mises stress criterion.
    VonMises,
    /// Thickness threshold.
    ThicknessThreshold,
}

/// Model for iceberg calving at marine-terminating glaciers.
#[derive(Debug, Clone)]
pub struct CalvingModel {
    /// Active calving criterion.
    pub criterion: CalvingCriterion,
    /// Critical Von Mises stress in Pa (for Von Mises criterion).
    pub critical_stress: f64,
    /// Minimum ice thickness before calving (m) (for thickness threshold).
    pub min_thickness: f64,
    /// Fracture toughness in Pa m^{1/2}.
    pub fracture_toughness: f64,
    /// Water depth at terminus in m.
    pub water_depth: f64,
    /// Water density in kg/m^3.
    pub water_density: f64,
}

impl CalvingModel {
    /// Create a new calving model with the given criterion.
    pub fn new(criterion: CalvingCriterion) -> Self {
        Self {
            criterion,
            critical_stress: 1.0e6,
            min_thickness: 50.0,
            fracture_toughness: 1.0e5,
            water_depth: 0.0,
            water_density: 1028.0,
        }
    }

    /// Set water depth at the terminus.
    pub fn set_water_depth(&mut self, depth: f64) {
        self.water_depth = depth;
    }

    /// Surface crevasse depth via the Nye formulation.
    ///
    /// d_s = R_xx / (rho_i * g)  where R_xx is the tensile stress.
    pub fn surface_crevasse_depth(tensile_stress: f64) -> f64 {
        if tensile_stress <= 0.0 {
            return 0.0;
        }
        tensile_stress / (ICE_DENSITY * GRAVITY)
    }

    /// Basal crevasse depth accounting for water pressure.
    ///
    /// d_b = (rho_w * g * d_w - sigma_xx) / ((rho_w - rho_i) * g)
    pub fn basal_crevasse_depth(&self, compressive_stress: f64) -> f64 {
        let rho_diff = self.water_density - ICE_DENSITY;
        if rho_diff.abs() < 1.0e-6 {
            return 0.0;
        }
        let numerator = self.water_density * GRAVITY * self.water_depth + compressive_stress;
        let depth = numerator / (rho_diff * GRAVITY);
        depth.max(0.0)
    }

    /// Benn calving criterion: calving occurs when d_s + d_b >= H.
    pub fn benn_calving_check(
        &self,
        thickness: f64,
        tensile_stress: f64,
        compressive_stress: f64,
    ) -> bool {
        let d_s = Self::surface_crevasse_depth(tensile_stress);
        let d_b = self.basal_crevasse_depth(compressive_stress);
        d_s + d_b >= thickness
    }

    /// Von Mises effective stress from a stress tensor.
    pub fn von_mises_stress(stress: &[f64; 6]) -> f64 {
        let s = stress;
        let val = 0.5
            * ((s[0] - s[1]).powi(2)
                + (s[1] - s[2]).powi(2)
                + (s[2] - s[0]).powi(2)
                + 6.0 * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2)));
        val.sqrt()
    }

    /// Check whether calving occurs for a particle.
    pub fn should_calve(&self, particle: &IceParticle, thickness: f64) -> bool {
        match self.criterion {
            CalvingCriterion::CrevasseDepth => {
                let tensile = particle.stress[0].max(0.0);
                let compressive = (-particle.stress[0]).max(0.0);
                self.benn_calving_check(thickness, tensile, compressive)
            }
            CalvingCriterion::VonMises => {
                let vm = Self::von_mises_stress(&particle.stress);
                vm >= self.critical_stress
            }
            CalvingCriterion::ThicknessThreshold => thickness < self.min_thickness,
        }
    }

    /// Calving rate estimate based on Von Mises stress (Morlighem et al.).
    pub fn calving_rate_von_mises(&self, von_mises: f64) -> f64 {
        if von_mises <= 0.0 {
            return 0.0;
        }
        // Simplified calving rate: c = v * (sigma_vm / sigma_crit)^2
        let ratio = von_mises / self.critical_stress;
        ratio * ratio
    }
}

// ---------------------------------------------------------------------------
// ThermalIce
// ---------------------------------------------------------------------------

/// Thermal model for heat transport in an ice sheet.
///
/// Solves the 1-D vertical heat equation with strain heating,
/// geothermal flux, and surface temperature boundary conditions.
#[derive(Debug, Clone)]
pub struct ThermalIce {
    /// Thermal conductivity in W/(m K).
    pub conductivity: f64,
    /// Specific heat capacity in J/(kg K).
    pub specific_heat: f64,
    /// Geothermal heat flux at the bed in W/m^2.
    pub geothermal_flux: f64,
    /// Surface temperature in K.
    pub surface_temperature: f64,
    /// Number of vertical layers.
    pub n_layers: usize,
}

impl ThermalIce {
    /// Create a new thermal ice model.
    pub fn new(geothermal_flux: f64, surface_temperature: f64, n_layers: usize) -> Self {
        Self {
            conductivity: ICE_THERMAL_CONDUCTIVITY,
            specific_heat: ICE_SPECIFIC_HEAT,
            geothermal_flux,
            surface_temperature,
            n_layers: n_layers.max(3),
        }
    }

    /// Thermal diffusivity kappa = k / (rho * cp).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.conductivity / (ICE_DENSITY * self.specific_heat)
    }

    /// Steady-state temperature profile (linear with geothermal gradient).
    ///
    /// Returns temperatures from bed to surface for the given thickness.
    pub fn steady_state_profile(&self, thickness: f64) -> Vec<f64> {
        let n = self.n_layers;
        let dz = thickness / (n as f64 - 1.0);
        let grad = self.geothermal_flux / self.conductivity; // K/m
        let t_bed = self.surface_temperature + grad * thickness;
        let t_bed_clamped = t_bed.min(ICE_MELTING_POINT);
        (0..n)
            .map(|i| {
                let z = i as f64 * dz;
                let frac = z / thickness;
                // Linear interpolation bed -> surface
                t_bed_clamped * (1.0 - frac) + self.surface_temperature * frac
            })
            .collect()
    }

    /// Strain heating rate per unit volume in W/m^3.
    ///
    /// Q_strain = 2 * eta * eps_eff^2  (or equivalently tau_eff * eps_eff).
    pub fn strain_heating(effective_stress: f64, effective_strain_rate: f64) -> f64 {
        effective_stress * effective_strain_rate
    }

    /// Explicit Euler step for 1-D vertical heat equation.
    ///
    /// dT/dt = kappa * d^2T/dz^2 + Q/(rho*cp)
    pub fn evolve_temperature(
        &self,
        profile: &mut [f64],
        thickness: f64,
        strain_heat: &[f64],
        dt: f64,
    ) {
        let n = profile.len();
        if n < 3 {
            return;
        }
        let dz = thickness / (n as f64 - 1.0);
        let kappa = self.thermal_diffusivity();
        let coeff = kappa * dt / (dz * dz);

        let old = profile.to_vec();

        // Interior points
        for i in 1..n - 1 {
            let diffusion = coeff * (old[i + 1] - 2.0 * old[i] + old[i - 1]);
            let source = if i < strain_heat.len() {
                strain_heat[i] * dt / (ICE_DENSITY * self.specific_heat)
            } else {
                0.0
            };
            profile[i] = old[i] + diffusion + source;
            // Clamp to melting point
            profile[i] = profile[i].min(ICE_MELTING_POINT);
        }
        // Boundary conditions
        // Bed: geothermal flux (Neumann BC)
        profile[0] = profile[1] + self.geothermal_flux * dz / self.conductivity;
        profile[0] = profile[0].min(ICE_MELTING_POINT);
        // Surface: fixed temperature (Dirichlet BC)
        profile[n - 1] = self.surface_temperature;
    }

    /// Robin number (Peclet-like) for advection vs diffusion balance.
    pub fn peclet_number(&self, vertical_velocity: f64, thickness: f64) -> f64 {
        let kappa = self.thermal_diffusivity();
        if kappa < 1.0e-30 {
            return f64::INFINITY;
        }
        vertical_velocity.abs() * thickness / kappa
    }
}

// ---------------------------------------------------------------------------
// IceFabric
// ---------------------------------------------------------------------------

/// Crystal orientation fabric model for polycrystalline ice.
///
/// Tracks the c-axis distribution and computes enhancement factors
/// for Glen's flow law.
#[derive(Debug, Clone)]
pub struct IceFabric {
    /// c-axis orientations as (theta, phi) in radians — polar coordinates on the unit sphere.
    pub orientations: Vec<[f64; 2]>,
    /// Eigenvalues of the second-order orientation tensor (a1 <= a2 <= a3).
    pub eigenvalues: [f64; 3],
    /// Enhancement factor for shear deformation.
    pub enhancement_shear: f64,
    /// Enhancement factor for compression.
    pub enhancement_compression: f64,
}

impl IceFabric {
    /// Create an isotropic fabric with n crystals.
    pub fn isotropic(n: usize) -> Self {
        let mut rng = rand::rng();
        let orientations: Vec<[f64; 2]> = (0..n)
            .map(|_| {
                let theta: f64 = rng.random_range(0.0..PI);
                let phi: f64 = rng.random_range(0.0..2.0 * PI);
                [theta, phi]
            })
            .collect();
        let mut fab = Self {
            orientations,
            eigenvalues: [1.0 / 3.0; 3],
            enhancement_shear: 1.0,
            enhancement_compression: 1.0,
        };
        fab.update_eigenvalues();
        fab
    }

    /// Create a single-maximum fabric (all c-axes aligned vertically).
    pub fn single_maximum(n: usize) -> Self {
        let orientations = vec![[0.0, 0.0]; n]; // theta=0 => vertical
        let mut fab = Self {
            orientations,
            eigenvalues: [0.0, 0.0, 1.0],
            enhancement_shear: 1.0,
            enhancement_compression: 1.0,
        };
        fab.update_eigenvalues();
        fab.compute_enhancement();
        fab
    }

    /// Compute the second-order orientation tensor and its eigenvalues.
    pub fn update_eigenvalues(&mut self) {
        let n = self.orientations.len();
        if n == 0 {
            self.eigenvalues = [0.0; 3];
            return;
        }
        // Build orientation tensor a2 = <c_i c_j>
        let mut a = [[0.0f64; 3]; 3];
        for orient in &self.orientations {
            let (st, ct) = orient[0].sin_cos();
            let (sp, cp) = orient[1].sin_cos();
            let c = [st * cp, st * sp, ct];
            for ii in 0..3 {
                for jj in 0..3 {
                    a[ii][jj] += c[ii] * c[jj];
                }
            }
        }
        let inv_n = 1.0 / n as f64;
        for row in &mut a {
            for v in row.iter_mut() {
                *v *= inv_n;
            }
        }
        // Eigenvalues of 3x3 symmetric tensor (diagonal dominance approximation)
        let mut eigs = [a[0][0], a[1][1], a[2][2]];
        eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        self.eigenvalues = eigs;
    }

    /// Compute enhancement factors from eigenvalues.
    ///
    /// Single maximum fabric enhances shear by up to ~10x,
    /// while reducing compression.
    pub fn compute_enhancement(&mut self) {
        let a3 = self.eigenvalues[2];
        // Empirical formula (Thorsteinsson 2001 simplified)
        self.enhancement_shear = 1.0 + 9.0 * a3 * a3;
        self.enhancement_compression = 1.0 / (1.0 + 9.0 * a3 * a3);
    }

    /// Woodcock parameter: k = ln(a3/a2) / ln(a2/a1).
    ///
    /// k > 1 => single maximum (girdle if k < 1).
    pub fn woodcock_parameter(&self) -> f64 {
        let [a1, a2, a3] = self.eigenvalues;
        let a1c = a1.max(1.0e-10);
        let a2c = a2.max(1.0e-10);
        let num = (a3 / a2c).ln();
        let den = (a2c / a1c).ln();
        if den.abs() < 1.0e-12 {
            return 1.0;
        }
        num / den
    }

    /// Rotate fabric towards vertical under compressive strain.
    pub fn apply_compression(&mut self, strain_increment: f64) {
        for orient in &mut self.orientations {
            // Rotate theta towards 0 (vertical)
            orient[0] *= (1.0 - strain_increment.abs()).max(0.0);
        }
        self.update_eigenvalues();
        self.compute_enhancement();
    }
}

// ---------------------------------------------------------------------------
// GlacialErosion
// ---------------------------------------------------------------------------

/// Model for glacial erosion processes: abrasion and plucking.
#[derive(Debug, Clone)]
pub struct GlacialErosion {
    /// Abrasion coefficient (dimensionless).
    pub abrasion_coefficient: f64,
    /// Abrasion exponent on sliding velocity.
    pub abrasion_exponent: f64,
    /// Plucking coefficient.
    pub plucking_coefficient: f64,
    /// Plucking exponent on effective pressure.
    pub plucking_exponent: f64,
    /// Sediment transport efficiency \[0, 1\].
    pub transport_efficiency: f64,
    /// Accumulated erosion depth in metres.
    pub total_erosion: f64,
    /// Accumulated sediment volume per unit width in m^2.
    pub sediment_volume: f64,
}

impl GlacialErosion {
    /// Create a new glacial erosion model.
    pub fn new(abrasion_coeff: f64, plucking_coeff: f64) -> Self {
        Self {
            abrasion_coefficient: abrasion_coeff,
            abrasion_exponent: 2.0,
            plucking_coefficient: plucking_coeff,
            plucking_exponent: 1.0,
            transport_efficiency: 0.5,
            total_erosion: 0.0,
            sediment_volume: 0.0,
        }
    }

    /// Abrasion rate in m/yr.
    ///
    /// E_a = K_a * u_b^l  where u_b is basal sliding velocity.
    pub fn abrasion_rate(&self, sliding_velocity: f64) -> f64 {
        self.abrasion_coefficient * sliding_velocity.abs().powf(self.abrasion_exponent)
    }

    /// Plucking rate in m/yr.
    ///
    /// E_p = K_p * N^m  where N is effective pressure (ice overburden - water pressure).
    pub fn plucking_rate(&self, effective_pressure: f64) -> f64 {
        if effective_pressure <= 0.0 {
            return 0.0;
        }
        self.plucking_coefficient * effective_pressure.powf(self.plucking_exponent)
    }

    /// Total erosion rate in m/yr.
    pub fn total_erosion_rate(&self, sliding_velocity: f64, effective_pressure: f64) -> f64 {
        self.abrasion_rate(sliding_velocity) + self.plucking_rate(effective_pressure)
    }

    /// Update erosion and sediment for a time step (in years).
    pub fn step(&mut self, sliding_velocity: f64, effective_pressure: f64, dt_years: f64) {
        let rate = self.total_erosion_rate(sliding_velocity, effective_pressure);
        let eroded = rate * dt_years;
        self.total_erosion += eroded;
        self.sediment_volume += eroded * self.transport_efficiency;
    }

    /// Basal shear stress from ice thickness and surface slope.
    pub fn basal_shear_stress(thickness: f64, surface_slope: f64) -> f64 {
        ICE_DENSITY * GRAVITY * thickness * surface_slope.abs()
    }

    /// Sediment flux estimate (simplified diffusive transport).
    pub fn sediment_flux(&self, sliding_velocity: f64) -> f64 {
        self.sediment_volume * sliding_velocity.abs() * self.transport_efficiency
    }
}

// ---------------------------------------------------------------------------
// IceAnalysis
// ---------------------------------------------------------------------------

/// Analysis tools for ice sheet simulations.
#[derive(Debug, Clone)]
pub struct IceAnalysis {
    /// Recorded total ice volume over time (time, volume).
    pub volume_history: Vec<[f64; 2]>,
    /// Mass balance record (time, smb, calving, basal_melt).
    pub mass_balance_history: Vec<[f64; 4]>,
    /// Grounding line positions over time (time, x_gl).
    pub grounding_line_history: Vec<[f64; 2]>,
}

impl IceAnalysis {
    /// Create a new analysis tracker.
    pub fn new() -> Self {
        Self {
            volume_history: Vec::new(),
            mass_balance_history: Vec::new(),
            grounding_line_history: Vec::new(),
        }
    }

    /// Record ice volume at a given time.
    pub fn record_volume(&mut self, time: f64, volume: f64) {
        self.volume_history.push([time, volume]);
    }

    /// Record mass balance components.
    pub fn record_mass_balance(&mut self, time: f64, smb: f64, calving: f64, basal_melt: f64) {
        self.mass_balance_history
            .push([time, smb, calving, basal_melt]);
    }

    /// Record grounding line position.
    pub fn record_grounding_line(&mut self, time: f64, x_gl: f64) {
        self.grounding_line_history.push([time, x_gl]);
    }

    /// Compute total ice volume from particle data.
    pub fn compute_volume(particles: &[IceParticle]) -> f64 {
        particles.iter().map(|p| p.mass / p.density).sum()
    }

    /// Compute mean ice temperature.
    pub fn mean_temperature(particles: &[IceParticle]) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }
        let sum: f64 = particles.iter().map(|p| p.temperature).sum();
        sum / particles.len() as f64
    }

    /// Compute mean ice speed.
    pub fn mean_speed(particles: &[IceParticle]) -> f64 {
        if particles.is_empty() {
            return 0.0;
        }
        let sum: f64 = particles.iter().map(|p| p.speed()).sum();
        sum / particles.len() as f64
    }

    /// Find grounding line position from a flotation criterion.
    ///
    /// The grounding line is where ice thickness equals flotation thickness.
    pub fn grounding_line_position(
        thickness: &[f64],
        bed_elevation: &[f64],
        dx: f64,
        water_density: f64,
    ) -> Option<f64> {
        if thickness.len() != bed_elevation.len() || thickness.is_empty() {
            return None;
        }
        let rho_ratio = water_density / ICE_DENSITY;
        for i in 0..thickness.len() - 1 {
            let water_depth = (-bed_elevation[i]).max(0.0);
            let flotation_thickness = water_depth * rho_ratio;
            let water_depth_next = (-bed_elevation[i + 1]).max(0.0);
            let flotation_next = water_depth_next * rho_ratio;

            let diff_i = thickness[i] - flotation_thickness;
            let diff_next = thickness[i + 1] - flotation_next;
            if diff_i >= 0.0 && diff_next < 0.0 {
                // Linear interpolation
                let frac = diff_i / (diff_i - diff_next);
                return Some(i as f64 * dx + frac * dx);
            }
        }
        None
    }

    /// Volume change rate (m^3/yr) between the last two records.
    pub fn volume_change_rate(&self) -> Option<f64> {
        let n = self.volume_history.len();
        if n < 2 {
            return None;
        }
        let [t1, v1] = self.volume_history[n - 2];
        let [t2, v2] = self.volume_history[n - 1];
        let dt = t2 - t1;
        if dt.abs() < 1.0e-30 {
            return None;
        }
        Some((v2 - v1) / dt)
    }

    /// Net mass balance from the last record (SMB - calving - basal melt).
    pub fn net_mass_balance(&self) -> Option<f64> {
        self.mass_balance_history.last().map(|r| r[1] - r[2] - r[3])
    }

    /// Sea level equivalent contribution in metres.
    ///
    /// SLE = delta_V * rho_i / (rho_w * A_ocean)
    pub fn sea_level_equivalent(volume_change: f64, ocean_area: f64, water_density: f64) -> f64 {
        if ocean_area.abs() < 1.0e-10 {
            return 0.0;
        }
        volume_change * ICE_DENSITY / (water_density * ocean_area)
    }
}

impl Default for IceAnalysis {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// SPH kernel helper
// ---------------------------------------------------------------------------

/// Cubic spline kernel (3D) for ice SPH.
pub fn cubic_spline_kernel_3d(r: f64, h: f64) -> f64 {
    if h <= 0.0 || r < 0.0 {
        return 0.0;
    }
    let q = r / h;
    let norm = 1.0 / (PI * h * h * h);
    if q <= 1.0 {
        norm * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q <= 2.0 {
        norm * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}

/// Gradient magnitude of cubic spline kernel (3D).
pub fn cubic_spline_gradient_3d(r: f64, h: f64) -> f64 {
    if h <= 0.0 || r < 1.0e-30 {
        return 0.0;
    }
    let q = r / h;
    let norm = 1.0 / (PI * h * h * h * h);
    if q <= 1.0 {
        norm * (-3.0 * q + 2.25 * q * q)
    } else if q <= 2.0 {
        norm * (-0.75 * (2.0 - q).powi(2))
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Ice SPH step
// ---------------------------------------------------------------------------

/// Perform one SPH time step for ice particles.
///
/// Simplified leapfrog integration with viscous forces from Glen's flow law.
#[allow(clippy::too_many_arguments)]
pub fn ice_sph_step(
    particles: &mut [IceParticle],
    glen: &GlenFlowLaw,
    dt: f64,
    surface_slope: f64,
    thickness: f64,
    _gravity_vec: [f64; 3],
) {
    let n = particles.len();
    // Compute driving stress
    let tau_d = IceSheetFlow::driving_stress(thickness, surface_slope);

    for i in 0..n {
        let t = particles[i].temperature;
        let sr = glen.strain_rate(tau_d, t);
        let eta = glen.effective_viscosity(tau_d, t);
        particles[i].viscosity = eta;

        // Simple gravity-driven acceleration
        let acc_x = tau_d / (ICE_DENSITY * thickness.max(1.0));
        particles[i].vel[0] += acc_x * dt;
        // Viscous damping
        let damping = (-dt * ICE_DENSITY * GRAVITY / eta.max(1.0)).exp();
        particles[i].vel[0] *= damping;

        // Update position
        for d in 0..3 {
            particles[i].pos[d] += particles[i].vel[d] * dt;
        }

        // Update stress (simplified)
        particles[i].stress[0] = 2.0 * eta * sr;

        // Age
        particles[i].age += dt / SECONDS_PER_YEAR;
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
    fn test_ice_particle_creation() {
        let p = IceParticle::new([0.0, 0.0, 100.0], 260.0);
        assert!((p.temperature - 260.0).abs() < TOL);
        assert!((p.density - ICE_DENSITY).abs() < TOL);
        assert!(p.speed() < TOL);
    }

    #[test]
    fn test_overburden_pressure() {
        let p = IceParticle::new([0.0, 0.0, 500.0], 260.0);
        let op = p.overburden_pressure(1000.0);
        let expected = ICE_DENSITY * GRAVITY * 500.0;
        assert!((op - expected).abs() / expected < 1.0e-10);
    }

    #[test]
    fn test_pressure_melting_point() {
        let p = IceParticle::new([0.0, 0.0, 0.0], 260.0);
        let pmp = p.pressure_melting_point(1000.0);
        // Should be below 273.15
        assert!(pmp < ICE_MELTING_POINT);
        assert!(pmp > 260.0);
    }

    #[test]
    fn test_glen_flow_law_default() {
        let g = GlenFlowLaw::default();
        assert!((g.n - 3.0).abs() < TOL);
        assert!((g.enhancement_factor - 1.0).abs() < TOL);
    }

    #[test]
    fn test_rate_factor_cold() {
        let g = GlenFlowLaw::default();
        let a_cold = g.rate_factor(250.0);
        let a_warm = g.rate_factor(270.0);
        // Warmer ice should have larger rate factor
        assert!(a_warm > a_cold);
    }

    #[test]
    fn test_glen_strain_rate_increases_with_stress() {
        let g = GlenFlowLaw::default();
        let sr1 = g.strain_rate(1.0e5, 260.0);
        let sr2 = g.strain_rate(2.0e5, 260.0);
        // With n=3, doubling stress => 8x strain rate
        let ratio = sr2 / sr1;
        assert!((ratio - 8.0).abs() < 0.1);
    }

    #[test]
    fn test_glen_effective_viscosity_positive() {
        let g = GlenFlowLaw::default();
        let eta = g.effective_viscosity(1.0e5, 260.0);
        assert!(eta > 0.0);
        assert!(eta.is_finite());
    }

    #[test]
    fn test_glen_stress_strain_roundtrip() {
        let g = GlenFlowLaw::default();
        let tau = 1.0e5;
        let temp = 260.0;
        let sr = g.strain_rate(tau, temp);
        let tau_back = g.effective_stress(sr, temp);
        assert!((tau - tau_back).abs() / tau < 1.0e-8);
    }

    #[test]
    fn test_second_invariant_uniaxial() {
        // Pure shear: eps_xy = e, all others zero
        let eps = [0.0, 0.0, 0.0, 0.01, 0.0, 0.0];
        let inv = second_invariant_strain_rate(&eps);
        // sqrt(0.5 * 2 * 0.01^2) = 0.01
        assert!((inv - 0.01).abs() < 1.0e-10);
    }

    #[test]
    fn test_driving_stress() {
        let tau = IceSheetFlow::driving_stress(1000.0, 0.01);
        let expected = ICE_DENSITY * GRAVITY * 1000.0 * 0.01;
        assert!((tau - expected).abs() < 1.0e-6);
    }

    #[test]
    fn test_sia_velocity_profile_increases_with_height() {
        let flow = IceSheetFlow::new(GlenFlowLaw::default());
        let h = 1000.0;
        let slope = 0.01;
        let temp = 260.0;
        let v_bed = flow.velocity_at_depth(0.0, h, slope, temp);
        let v_mid = flow.velocity_at_depth(500.0, h, slope, temp);
        let v_top = flow.velocity_at_depth(1000.0, h, slope, temp);
        // Velocity should increase with height
        assert!(v_mid > v_bed || (v_bed.abs() < 1.0e-30 && v_mid >= 0.0));
        assert!(v_top >= v_mid);
    }

    #[test]
    fn test_sia_deformation_velocity_positive() {
        let flow = IceSheetFlow::new(GlenFlowLaw::default());
        let v = flow.sia_deformation_velocity(1000.0, 0.01, 260.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_basal_sliding() {
        let mut flow = IceSheetFlow::new(GlenFlowLaw::default());
        flow.set_sliding(1.0e-10, 1.0);
        let v_b = flow.basal_sliding_velocity(1000.0, 0.01);
        assert!(v_b > 0.0);
    }

    #[test]
    fn test_ice_flux_positive() {
        let flow = IceSheetFlow::new(GlenFlowLaw::default());
        let q = flow.ice_flux(1000.0, 0.01, 260.0);
        assert!(q > 0.0);
    }

    #[test]
    fn test_thickness_change_rate() {
        let dh = IceSheetFlow::thickness_change_rate(100.0, 110.0, 1000.0, 0.5);
        // dh = 0.5 - (110-100)/1000 = 0.5 - 0.01 = 0.49
        assert!((dh - 0.49).abs() < 1.0e-10);
    }

    #[test]
    fn test_surface_crevasse_depth() {
        let d = CalvingModel::surface_crevasse_depth(1.0e5);
        let expected = 1.0e5 / (ICE_DENSITY * GRAVITY);
        assert!((d - expected).abs() < 1.0e-6);
    }

    #[test]
    fn test_surface_crevasse_depth_no_tension() {
        let d = CalvingModel::surface_crevasse_depth(-1.0e5);
        assert!((d - 0.0).abs() < TOL);
    }

    #[test]
    fn test_von_mises_stress_uniaxial() {
        // Uniaxial tension: sigma_xx = S, rest zero
        let s = 1.0e6;
        let stress = [s, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = CalvingModel::von_mises_stress(&stress);
        assert!((vm - s).abs() / s < 1.0e-10);
    }

    #[test]
    fn test_calving_von_mises_criterion() {
        let model = CalvingModel::new(CalvingCriterion::VonMises);
        let mut p = IceParticle::new([0.0, 0.0, 0.0], 260.0);
        p.stress = [2.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(model.should_calve(&p, 100.0));
    }

    #[test]
    fn test_calving_thickness_threshold() {
        let model = CalvingModel::new(CalvingCriterion::ThicknessThreshold);
        let p = IceParticle::new([0.0, 0.0, 0.0], 260.0);
        assert!(model.should_calve(&p, 30.0)); // below 50m threshold
        assert!(!model.should_calve(&p, 100.0));
    }

    #[test]
    fn test_thermal_steady_state() {
        let thermal = ThermalIce::new(0.05, 243.15, 11);
        let profile = thermal.steady_state_profile(1000.0);
        // Surface should be at surface_temperature
        assert!((profile[10] - 243.15).abs() < TOL);
        // Bed should be warmer
        assert!(profile[0] > profile[10]);
    }

    #[test]
    fn test_thermal_diffusivity_positive() {
        let thermal = ThermalIce::new(0.05, 243.15, 11);
        assert!(thermal.thermal_diffusivity() > 0.0);
    }

    #[test]
    fn test_strain_heating_positive() {
        let q = ThermalIce::strain_heating(1.0e5, 1.0e-10);
        assert!(q > 0.0);
        assert!((q - 1.0e-5).abs() < 1.0e-15);
    }

    #[test]
    fn test_fabric_isotropic_eigenvalues() {
        let fab = IceFabric::isotropic(10000);
        // For isotropic fabric, eigenvalues should sum to ~1 and be broadly similar
        let sum: f64 = fab.eigenvalues.iter().sum();
        assert!((sum - 1.0).abs() < 0.15);
        // The largest eigenvalue should not dominate
        assert!(fab.eigenvalues[2] < 0.6);
    }

    #[test]
    fn test_fabric_single_maximum() {
        let fab = IceFabric::single_maximum(100);
        // a3 should be close to 1
        assert!(fab.eigenvalues[2] > 0.9);
    }

    #[test]
    fn test_fabric_compression_aligns() {
        let mut fab = IceFabric::isotropic(1000);
        let a3_before = fab.eigenvalues[2];
        for _ in 0..20 {
            fab.apply_compression(0.05);
        }
        assert!(fab.eigenvalues[2] > a3_before);
    }

    #[test]
    fn test_erosion_abrasion_rate() {
        let e = GlacialErosion::new(5.0e-4, 1.0e-7);
        let rate = e.abrasion_rate(100.0); // 100 m/yr sliding
        assert!(rate > 0.0);
    }

    #[test]
    fn test_erosion_plucking_zero_pressure() {
        let e = GlacialErosion::new(5.0e-4, 1.0e-7);
        let rate = e.plucking_rate(0.0);
        assert!((rate - 0.0).abs() < TOL);
    }

    #[test]
    fn test_erosion_step_accumulates() {
        let mut e = GlacialErosion::new(1.0e-3, 0.0);
        e.step(10.0, 0.0, 1.0);
        assert!(e.total_erosion > 0.0);
        assert!(e.sediment_volume > 0.0);
    }

    #[test]
    fn test_analysis_volume() {
        let particles = vec![
            IceParticle::new([0.0, 0.0, 0.0], 260.0),
            IceParticle::new([1.0, 0.0, 0.0], 260.0),
        ];
        let vol = IceAnalysis::compute_volume(&particles);
        assert!(vol > 0.0);
    }

    #[test]
    fn test_analysis_grounding_line() {
        let thickness = vec![500.0, 400.0, 300.0, 200.0, 100.0];
        let bed = vec![0.0, -100.0, -200.0, -300.0, -400.0];
        let pos = IceAnalysis::grounding_line_position(&thickness, &bed, 1000.0, 1028.0);
        // Should find a grounding line somewhere
        assert!(pos.is_some());
    }

    #[test]
    fn test_analysis_volume_change_rate() {
        let mut a = IceAnalysis::new();
        a.record_volume(0.0, 1000.0);
        a.record_volume(1.0, 900.0);
        let rate = a.volume_change_rate().unwrap();
        assert!((rate - (-100.0)).abs() < TOL);
    }

    #[test]
    fn test_cubic_spline_kernel() {
        let w = cubic_spline_kernel_3d(0.0, 1.0);
        assert!(w > 0.0);
        let w2 = cubic_spline_kernel_3d(3.0, 1.0);
        assert!((w2 - 0.0).abs() < TOL);
    }

    #[test]
    fn test_ice_sph_step_moves_particles() {
        let mut particles = vec![IceParticle::new([0.0, 0.0, 500.0], 260.0)];
        let glen = GlenFlowLaw::default();
        ice_sph_step(
            &mut particles,
            &glen,
            1.0,
            0.01,
            1000.0,
            [0.0, 0.0, -GRAVITY],
        );
        // Particle should have moved in x
        assert!(particles[0].pos[0].abs() > 0.0 || particles[0].vel[0].abs() > 0.0);
    }

    #[test]
    fn test_sea_level_equivalent() {
        let sle = IceAnalysis::sea_level_equivalent(1.0e12, 3.625e14, 1028.0);
        assert!(sle > 0.0);
        assert!(sle.is_finite());
    }
}
