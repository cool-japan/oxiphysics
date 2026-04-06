// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Nanoscale SPH module with surface and molecular effects.
//!
//! Models molecular-continuum coupling, surface forces, and thermal fluctuations:
//!
//! - [`NanoscaleSph`]: Knudsen number effects, slip boundary, molecular-continuum coupling
//! - [`SurfaceTensionNano`]: van der Waals energy, Young–Dupré equation, disjoining pressure
//! - [`ThermalFluctuationSph`]: Landau-Lifshitz fluctuating hydrodynamics, Brownian SPH
//! - [`ElectricDoubleLayer`]: Debye length, EDL formation, electro-osmotic flow, streaming potential
//! - [`MolecularConfinement`]: density oscillations near walls, flow enhancement
//! - [`NanobubbleSph`]: nanobubble stability, surface and bulk nanobubbles, dissolution kinetics
//! - [`NanofluidThermal`]: Maxwell model, Bruggeman model, nanoparticle thermal enhancement
//! - [`SlipFlowModel`]: Navier slip, Maxwell slip coefficient, accommodation coefficient
//! - [`DissipativeParticleNano`]: DPD-SPH hybrid with conservative/dissipative/random forces
//! - [`NonequilibriumNano`]: NEMD-SPH, heat flux, shear stress, Green-Kubo relations

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use rand::RngExt;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
const K_B: f64 = 1.380_649e-23;

/// Avogadro's number (mol⁻¹).
const N_AV: f64 = 6.022_140_76e23;

/// Elementary charge (C).
const E_CHARGE: f64 = 1.602_176_634e-19;

/// Permittivity of free space (F m⁻¹).
const EPS_0: f64 = 8.854_187_817e-12;

/// Reference temperature (K).
const T_REF: f64 = 298.15;

/// Dielectric constant of water (dimensionless).
const EPS_WATER: f64 = 78.5;

/// Universal gas constant (J mol⁻¹ K⁻¹).
const GAS_CONSTANT: f64 = 8.314_462_618;

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
pub fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Normalise a 3-vector.
#[inline]
pub fn normalise3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1.0e-300 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
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
// NanoscaleSph
// ---------------------------------------------------------------------------

/// Knudsen number flow regime classification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KnudsenRegime {
    /// Kn < 0.001 — continuum flow.
    Continuum,
    /// 0.001 ≤ Kn < 0.1 — slip flow.
    SlipFlow,
    /// 0.1 ≤ Kn < 10 — transitional flow.
    Transitional,
    /// Kn ≥ 10 — free molecular flow.
    FreeMolecular,
}

/// Nanoscale SPH particle with molecular-continuum coupling state.
#[derive(Debug, Clone)]
pub struct NanoscaleParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Local density (kg m⁻³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Particle mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub h: f64,
}

impl NanoscaleParticle {
    /// Create a new nanoscale SPH particle.
    pub fn new(position: [f64; 3], mass: f64, h: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            density: 1000.0,
            pressure: 0.0,
            temperature: T_REF,
            mass,
            h,
        }
    }

    /// Kinetic energy of this particle (J).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }

    /// Thermal de Broglie wavelength (m) for a particle of molecular mass `m_mol` (kg).
    pub fn thermal_de_broglie(&self, m_mol: f64) -> f64 {
        let h_planck = 6.626_070_15e-34;
        h_planck / (2.0 * PI * m_mol * K_B * self.temperature).sqrt()
    }
}

/// Main nanoscale SPH solver coupling molecular and continuum descriptions.
#[derive(Debug, Clone)]
pub struct NanoscaleSph {
    /// Particle data.
    pub particles: Vec<NanoscaleParticle>,
    /// Mean free path of the fluid molecules (m).
    pub mean_free_path: f64,
    /// Characteristic length scale of the channel/domain (m).
    pub characteristic_length: f64,
    /// Dynamic viscosity (Pa s).
    pub viscosity: f64,
    /// Fluid thermal conductivity (W m⁻¹ K⁻¹).
    pub conductivity: f64,
    /// Fluid density (kg m⁻³).
    pub fluid_density: f64,
}

impl NanoscaleSph {
    /// Create a new `NanoscaleSph` solver.
    pub fn new(
        particles: Vec<NanoscaleParticle>,
        mean_free_path: f64,
        characteristic_length: f64,
        viscosity: f64,
        conductivity: f64,
        fluid_density: f64,
    ) -> Self {
        Self {
            particles,
            mean_free_path,
            characteristic_length,
            viscosity,
            conductivity,
            fluid_density,
        }
    }

    /// Knudsen number Kn = λ / L.
    pub fn knudsen_number(&self) -> f64 {
        self.mean_free_path / self.characteristic_length
    }

    /// Classify the current flow regime based on the Knudsen number.
    pub fn flow_regime(&self) -> KnudsenRegime {
        let kn = self.knudsen_number();
        if kn < 0.001 {
            KnudsenRegime::Continuum
        } else if kn < 0.1 {
            KnudsenRegime::SlipFlow
        } else if kn < 10.0 {
            KnudsenRegime::Transitional
        } else {
            KnudsenRegime::FreeMolecular
        }
    }

    /// Effective viscosity with Knudsen correction (Beskok-Karniadakis model).
    ///
    /// μ_eff = μ₀ / (1 − b·Kn)  where b ≈ −1 for slip regime.
    pub fn effective_viscosity(&self) -> f64 {
        let kn = self.knudsen_number();
        let b = -1.0;
        self.viscosity / (1.0 - b * kn)
    }

    /// Velocity slip length (m) = μ / (density * (π/8) * v_th) approximately.
    ///
    /// Slip length β ≈ 2 λ (1 − σ_v) / σ_v   where σ_v = tangential momentum accommodation.
    pub fn slip_length(&self, accommodation_coefficient: f64) -> f64 {
        2.0 * self.mean_free_path * (1.0 - accommodation_coefficient) / accommodation_coefficient
    }

    /// Number of particles.
    pub fn n_particles(&self) -> usize {
        self.particles.len()
    }

    /// Mean kinetic temperature (K) from particle velocities.
    pub fn mean_temperature(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        self.particles.iter().map(|p| p.temperature).sum::<f64>() / self.particles.len() as f64
    }
}

// ---------------------------------------------------------------------------
// SurfaceTensionNano
// ---------------------------------------------------------------------------

/// van der Waals / nanoscale surface tension model.
#[derive(Debug, Clone)]
pub struct SurfaceTensionNano {
    /// Liquid–vapour surface energy γ_LV (J m⁻²).
    pub gamma_lv: f64,
    /// Solid–liquid surface energy γ_SL (J m⁻²).
    pub gamma_sl: f64,
    /// Solid–vapour surface energy γ_SV (J m⁻²).
    pub gamma_sv: f64,
    /// Hamaker constant for van der Waals interaction (J).
    pub hamaker_constant: f64,
    /// Reference film thickness for disjoining pressure (m).
    pub reference_film_thickness: f64,
}

impl SurfaceTensionNano {
    /// Create a new nanoscale surface tension model.
    pub fn new(
        gamma_lv: f64,
        gamma_sl: f64,
        gamma_sv: f64,
        hamaker_constant: f64,
        reference_film_thickness: f64,
    ) -> Self {
        Self {
            gamma_lv,
            gamma_sl,
            gamma_sv,
            hamaker_constant,
            reference_film_thickness,
        }
    }

    /// Young contact angle (radians) from Young's equation:
    ///
    /// cos θ = (γ_SV − γ_SL) / γ_LV
    pub fn contact_angle(&self) -> f64 {
        let cos_theta = (self.gamma_sv - self.gamma_sl) / self.gamma_lv;
        cos_theta.clamp(-1.0, 1.0).acos()
    }

    /// Work of adhesion (J m⁻²) via Young–Dupré equation:
    ///
    /// W_adh = γ_LV * (1 + cos θ)
    pub fn work_of_adhesion(&self) -> f64 {
        self.gamma_lv * (1.0 + self.contact_angle().cos())
    }

    /// Disjoining pressure (Pa) for a thin film of thickness h (m) from van der Waals forces.
    ///
    /// Π_vdW(h) = −A_H / (6π h³)
    pub fn disjoining_pressure_vdw(&self, film_thickness: f64) -> f64 {
        -self.hamaker_constant / (6.0 * PI * film_thickness.powi(3))
    }

    /// Electrostatic component of disjoining pressure (Pa) using DLVO theory.
    ///
    /// Π_el(h) = 64 c₀ k_B T tanh²(ze ψ₀/(4 k_B T)) exp(−h/λ_D)
    pub fn disjoining_pressure_electrostatic(
        &self,
        film_thickness: f64,
        debye_length: f64,
        concentration: f64,
        temperature: f64,
        surface_potential: f64,
    ) -> f64 {
        let kbt = K_B * temperature;
        let tanh_term = (E_CHARGE * surface_potential / (4.0 * kbt)).tanh();
        64.0 * concentration * kbt * tanh_term * tanh_term * (-film_thickness / debye_length).exp()
    }

    /// Total DLVO disjoining pressure (Pa).
    pub fn disjoining_pressure_total(
        &self,
        film_thickness: f64,
        debye_length: f64,
        concentration: f64,
        temperature: f64,
        surface_potential: f64,
    ) -> f64 {
        self.disjoining_pressure_vdw(film_thickness)
            + self.disjoining_pressure_electrostatic(
                film_thickness,
                debye_length,
                concentration,
                temperature,
                surface_potential,
            )
    }

    /// Laplace pressure across a spherical nanodroplet of radius r (Pa).
    pub fn laplace_pressure(&self, radius: f64) -> f64 {
        2.0 * self.gamma_lv / radius
    }

    /// Critical disjoining pressure film thickness (m) — balance between electrostatic
    /// repulsion and van der Waals attraction.
    pub fn critical_film_thickness(
        &self,
        debye_length: f64,
        concentration: f64,
        temperature: f64,
        surface_potential: f64,
    ) -> f64 {
        // Simple bisection between 0.1 nm and 100 nm
        let mut lo = 1.0e-10_f64;
        let mut hi = 1.0e-7_f64;
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            let pi = self.disjoining_pressure_total(
                mid,
                debye_length,
                concentration,
                temperature,
                surface_potential,
            );
            if pi > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

// ---------------------------------------------------------------------------
// ThermalFluctuationSph
// ---------------------------------------------------------------------------

/// Landau-Lifshitz fluctuating hydrodynamics SPH model.
///
/// Adds a stochastic stress tensor to the SPH momentum equation
/// according to the fluctuation-dissipation theorem.
#[derive(Debug, Clone)]
pub struct ThermalFluctuationSph {
    /// Number of particles.
    pub n_particles: usize,
    /// Particle mass (kg).
    pub particle_mass: f64,
    /// Dynamic viscosity (Pa s).
    pub viscosity: f64,
    /// Thermal conductivity (W m⁻¹ K⁻¹).
    pub conductivity: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Particle volume (m³).
    pub particle_volume: f64,
    /// Stored random force realisations for each particle (N).
    pub random_forces: Vec<[f64; 3]>,
}

impl ThermalFluctuationSph {
    /// Create a new fluctuating hydrodynamics SPH model.
    pub fn new(
        n_particles: usize,
        particle_mass: f64,
        viscosity: f64,
        conductivity: f64,
        temperature: f64,
        particle_volume: f64,
    ) -> Self {
        Self {
            n_particles,
            particle_mass,
            viscosity,
            conductivity,
            temperature,
            particle_volume,
            random_forces: vec![[0.0; 3]; n_particles],
        }
    }

    /// Amplitude of thermal noise force per particle (N).
    ///
    /// σ_F = sqrt(2 k_B T μ V / dt)
    pub fn noise_amplitude(&self, dt: f64) -> f64 {
        (2.0 * K_B * self.temperature * self.viscosity * self.particle_volume / dt).sqrt()
    }

    /// Generate random thermal forces for all particles using the Landau-Lifshitz prescription.
    pub fn generate_random_forces(&mut self, dt: f64) {
        let mut rng = rand::rng();
        let sigma = self.noise_amplitude(dt);
        for force in self.random_forces.iter_mut() {
            // Box-Muller for Gaussian random numbers
            let u1: f64 = rng.random_range(1.0e-15..1.0);
            let u2: f64 = rng.random_range(0.0..2.0 * PI);
            let u3: f64 = rng.random_range(1.0e-15..1.0);
            let u4: f64 = rng.random_range(0.0..2.0 * PI);
            let z0 = (-2.0 * u1.ln()).sqrt() * u2.cos();
            let z1 = (-2.0 * u1.ln()).sqrt() * u2.sin();
            let z2 = (-2.0 * u3.ln()).sqrt() * u4.cos();
            force[0] = sigma * z0;
            force[1] = sigma * z1;
            force[2] = sigma * z2;
        }
    }

    /// Thermal diffusivity (m² s⁻¹).
    pub fn thermal_diffusivity(&self, density: f64, specific_heat: f64) -> f64 {
        self.conductivity / (density * specific_heat)
    }

    /// Mean square displacement from thermal fluctuations over time dt (m²).
    ///
    /// ⟨r²⟩ = 6 D dt   where D = k_B T / (3π μ d_p)
    pub fn mean_square_displacement(&self, dt: f64, particle_diameter: f64) -> f64 {
        let d_coeff = K_B * self.temperature / (3.0 * PI * self.viscosity * particle_diameter);
        6.0 * d_coeff * dt
    }

    /// Stokes-Einstein diffusion coefficient (m² s⁻¹) for particle of diameter `d`.
    pub fn stokes_einstein_diffusion(&self, particle_diameter: f64) -> f64 {
        K_B * self.temperature / (3.0 * PI * self.viscosity * particle_diameter)
    }
}

// ---------------------------------------------------------------------------
// ElectricDoubleLayer
// ---------------------------------------------------------------------------

/// Electric double layer (EDL) model for nanoscale electro-osmotic flows.
#[derive(Debug, Clone)]
pub struct ElectricDoubleLayer {
    /// Electrolyte concentration (mol m⁻³).
    pub concentration: f64,
    /// Ion valence z.
    pub valence: i32,
    /// Surface potential ψ₀ (V).
    pub surface_potential: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Relative permittivity of the electrolyte.
    pub relative_permittivity: f64,
}

impl ElectricDoubleLayer {
    /// Create a new EDL model.
    pub fn new(
        concentration: f64,
        valence: i32,
        surface_potential: f64,
        temperature: f64,
        relative_permittivity: f64,
    ) -> Self {
        Self {
            concentration,
            valence,
            surface_potential,
            temperature,
            relative_permittivity,
        }
    }

    /// Debye length λ_D (m).
    ///
    /// λ_D = sqrt(ε₀ ε_r k_B T / (2 N_A z² e² c))
    pub fn debye_length(&self) -> f64 {
        let c_si = self.concentration; // mol m⁻³
        let z2 = (self.valence * self.valence) as f64;
        let numerator = EPS_0 * self.relative_permittivity * K_B * self.temperature;
        let denominator = 2.0 * N_AV * z2 * E_CHARGE * E_CHARGE * c_si;
        (numerator / denominator).sqrt()
    }

    /// Electric potential at distance y from wall (V) using linearised Poisson-Boltzmann.
    ///
    /// ψ(y) = ψ₀ exp(−y / λ_D)
    pub fn potential_profile(&self, y: f64) -> f64 {
        self.surface_potential * (-y / self.debye_length()).exp()
    }

    /// Net charge density at distance y (C m⁻³).
    ///
    /// ρ_e(y) = −ε₀ε_r * d²ψ/dy² = ε₀ε_r ψ₀ / λ_D² * exp(−y/λ_D)
    pub fn charge_density(&self, y: f64) -> f64 {
        let lambda = self.debye_length();
        EPS_0 * self.relative_permittivity * self.surface_potential / (lambda * lambda)
            * (-y / lambda).exp()
    }

    /// Electro-osmotic mobility μ_EO (m² V⁻¹ s⁻¹) via Helmholtz-Smoluchowski.
    ///
    /// μ_EO = −ε₀ε_r ζ / μ   where ζ ≈ surface_potential
    pub fn electro_osmotic_mobility(&self, viscosity: f64) -> f64 {
        -EPS_0 * self.relative_permittivity * self.surface_potential / viscosity
    }

    /// Electro-osmotic velocity (m s⁻¹) for applied electric field E_app (V m⁻¹).
    pub fn electro_osmotic_velocity(&self, e_field: f64, viscosity: f64) -> f64 {
        self.electro_osmotic_mobility(viscosity) * e_field
    }

    /// Streaming potential (V) for a pressure-driven flow through a channel of half-height h.
    ///
    /// V_str = −(ε₀ε_r ζ / μ σ_bulk) * ΔP   (Onsager reciprocal)
    pub fn streaming_potential(&self, delta_p: f64, viscosity: f64, conductivity_bulk: f64) -> f64 {
        let mob = self.electro_osmotic_mobility(viscosity);
        -mob / conductivity_bulk * delta_p
    }

    /// Dimensionless EDL overlap parameter κh = h / λ_D.
    pub fn edl_overlap_parameter(&self, half_channel_height: f64) -> f64 {
        half_channel_height / self.debye_length()
    }
}

// ---------------------------------------------------------------------------
// MolecularConfinement
// ---------------------------------------------------------------------------

/// Molecular confinement model for nanoscale flow between parallel walls.
#[derive(Debug, Clone)]
pub struct MolecularConfinement {
    /// Channel half-width (m).
    pub channel_half_width: f64,
    /// Fluid molecule diameter (m).
    pub molecule_diameter: f64,
    /// Bulk fluid density (kg m⁻³).
    pub bulk_density: f64,
    /// Lennard-Jones energy parameter ε_LJ (J).
    pub lj_epsilon: f64,
    /// Lennard-Jones length parameter σ_LJ (m).
    pub lj_sigma: f64,
}

impl MolecularConfinement {
    /// Create a new molecular confinement model.
    pub fn new(
        channel_half_width: f64,
        molecule_diameter: f64,
        bulk_density: f64,
        lj_epsilon: f64,
        lj_sigma: f64,
    ) -> Self {
        Self {
            channel_half_width,
            molecule_diameter,
            bulk_density,
            lj_epsilon,
            lj_sigma,
        }
    }

    /// Number of molecular layers that fit in the channel.
    pub fn n_molecular_layers(&self) -> f64 {
        2.0 * self.channel_half_width / self.molecule_diameter
    }

    /// Density oscillation amplitude at distance y from the wall.
    ///
    /// ρ(y) = ρ_bulk * (1 + A * exp(−y/d_0) * cos(2π y / d_mol + φ))
    /// A ≈ 0.7, d_0 ≈ 2 * σ, d_mol ≈ σ, φ = 0
    pub fn density_profile(&self, y: f64) -> f64 {
        let amplitude = 0.7;
        let decay = 2.0 * self.lj_sigma;
        let period = self.lj_sigma;
        self.bulk_density * (1.0 + amplitude * (-y / decay).exp() * (2.0 * PI * y / period).cos())
    }

    /// Flow enhancement factor for water in carbon nanotubes (dimensionless).
    ///
    /// EF = 1 + (A / (1 + (d/d_c)^n)) where A ≈ 1000, d_c ≈ 1 nm, n ≈ 2.
    pub fn flow_enhancement_factor(&self) -> f64 {
        let d_nm = 2.0 * self.channel_half_width * 1.0e9;
        1.0 + 1000.0 / (1.0 + (d_nm / 1.0).powi(2))
    }

    /// Effective slip length (m) from flow enhancement.
    ///
    /// β_eff = r/4 * (EF - 1) for a cylindrical nanotube of radius r.
    pub fn effective_slip_length(&self) -> f64 {
        let r = self.channel_half_width;
        let ef = self.flow_enhancement_factor();
        r / 4.0 * (ef - 1.0)
    }

    /// Lennard-Jones fluid-wall interaction energy at distance r (J).
    pub fn lj_wall_energy(&self, r: f64) -> f64 {
        let sr6 = (self.lj_sigma / r).powi(6);
        4.0 * self.lj_epsilon * (sr6 * sr6 - sr6)
    }
}

// ---------------------------------------------------------------------------
// NanobubbleSph
// ---------------------------------------------------------------------------

/// Nanobubble type classification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NanobubbleType {
    /// Nanobubble pinned to a solid surface.
    Surface,
    /// Nanobubble freely suspended in the bulk liquid.
    Bulk,
}

/// Model for nanobubble stability and dissolution kinetics.
#[derive(Debug, Clone)]
pub struct NanobubbleSph {
    /// Nanobubble type.
    pub bubble_type: NanobubbleType,
    /// Bubble radius (m).
    pub radius: f64,
    /// Gas type (0 = air, 1 = O2, 2 = N2).
    pub gas_type: u8,
    /// Surface tension of liquid–gas interface (N m⁻¹).
    pub surface_tension: f64,
    /// Contact angle of surface nanobubble (radians, from liquid side).
    pub contact_angle: f64,
    /// Henry's law constant for the gas (mol m⁻³ Pa⁻¹).
    pub henry_constant: f64,
    /// Dissolved gas concentration in bulk liquid (mol m⁻³).
    pub bulk_concentration: f64,
}

impl NanobubbleSph {
    /// Create a new nanobubble model.
    pub fn new(
        bubble_type: NanobubbleType,
        radius: f64,
        gas_type: u8,
        surface_tension: f64,
        contact_angle: f64,
        henry_constant: f64,
        bulk_concentration: f64,
    ) -> Self {
        Self {
            bubble_type,
            radius,
            gas_type,
            surface_tension,
            contact_angle,
            henry_constant,
            bulk_concentration,
        }
    }

    /// Internal pressure of the nanobubble (Pa) above atmospheric pressure.
    ///
    /// ΔP = 2γ / R  (spherical) or  ΔP = 2γ sin(θ) / R  (surface bubble, cap geometry)
    pub fn internal_pressure(&self, p_atm: f64) -> f64 {
        match self.bubble_type {
            NanobubbleType::Bulk => p_atm + 2.0 * self.surface_tension / self.radius,
            NanobubbleType::Surface => {
                p_atm + 2.0 * self.surface_tension * self.contact_angle.sin() / self.radius
            }
        }
    }

    /// Saturated gas concentration at bubble interface (mol m⁻³).
    pub fn interface_concentration(&self, p_atm: f64) -> f64 {
        self.henry_constant * self.internal_pressure(p_atm)
    }

    /// Dissolution rate dr/dt (m s⁻¹) via diffusion-limited model.
    ///
    /// dr/dt = −D_g * (c_s - c_inf) * V_m / R
    /// where D_g is gas diffusivity in liquid and V_m is molar volume of gas.
    pub fn dissolution_rate(&self, diffusivity: f64, molar_volume: f64, p_atm: f64) -> f64 {
        let c_s = self.interface_concentration(p_atm);
        let driving_force = c_s - self.bulk_concentration;
        -diffusivity * driving_force * molar_volume / self.radius
    }

    /// Lifetime estimate (s) for a freely dissolving bulk nanobubble.
    ///
    /// t_life ≈ R² / (2 D V_m (c_s − c_inf))
    pub fn lifetime_estimate(&self, diffusivity: f64, molar_volume: f64, p_atm: f64) -> f64 {
        let c_s = self.interface_concentration(p_atm);
        let driving_force = (c_s - self.bulk_concentration).abs().max(1.0e-20);
        self.radius * self.radius / (2.0 * diffusivity * molar_volume * driving_force)
    }

    /// Radius of curvature for a surface nanobubble (spherical cap) (m).
    pub fn radius_of_curvature(&self) -> f64 {
        match self.bubble_type {
            NanobubbleType::Bulk => self.radius,
            NanobubbleType::Surface => self.radius / self.contact_angle.sin(),
        }
    }
}

// ---------------------------------------------------------------------------
// NanofluidThermal
// ---------------------------------------------------------------------------

/// Nanofluid thermal conductivity enhancement model.
#[derive(Debug, Clone)]
pub struct NanofluidThermal {
    /// Base fluid thermal conductivity (W m⁻¹ K⁻¹).
    pub k_fluid: f64,
    /// Nanoparticle thermal conductivity (W m⁻¹ K⁻¹).
    pub k_particle: f64,
    /// Nanoparticle volume fraction φ (0–1).
    pub volume_fraction: f64,
    /// Nanoparticle radius (m).
    pub particle_radius: f64,
    /// Interfacial thermal resistance (Kapitza) (m² K W⁻¹).
    pub kapitza_resistance: f64,
}

impl NanofluidThermal {
    /// Create a new nanofluid thermal model.
    pub fn new(
        k_fluid: f64,
        k_particle: f64,
        volume_fraction: f64,
        particle_radius: f64,
        kapitza_resistance: f64,
    ) -> Self {
        Self {
            k_fluid,
            k_particle,
            volume_fraction,
            particle_radius,
            kapitza_resistance,
        }
    }

    /// Effective particle conductivity accounting for Kapitza resistance.
    ///
    /// k_p_eff = k_p / (1 + R_k k_p / r)
    pub fn effective_particle_conductivity(&self) -> f64 {
        self.k_particle / (1.0 + self.kapitza_resistance * self.k_particle / self.particle_radius)
    }

    /// Maxwell model effective thermal conductivity (W m⁻¹ K⁻¹).
    ///
    /// k_eff / k_f = (k_p + 2 k_f + 2 φ(k_p − k_f)) / (k_p + 2 k_f − φ(k_p − k_f))
    pub fn maxwell_model(&self) -> f64 {
        let kp = self.effective_particle_conductivity();
        let kf = self.k_fluid;
        let phi = self.volume_fraction;
        kf * (kp + 2.0 * kf + 2.0 * phi * (kp - kf)) / (kp + 2.0 * kf - phi * (kp - kf))
    }

    /// Bruggeman model effective thermal conductivity (W m⁻¹ K⁻¹).
    ///
    /// Solves: (1−φ)(k_f − k_eff)/(k_f + 2 k_eff) + φ(k_p − k_eff)/(k_p + 2 k_eff) = 0
    pub fn bruggeman_model(&self) -> f64 {
        let kp = self.effective_particle_conductivity();
        let kf = self.k_fluid;
        let phi = self.volume_fraction;
        // Quadratic form: a k_eff² + b k_eff + c = 0
        let _a = 3.0 * phi - 1.0;
        let _b = -(3.0 * phi * kp + (2.0 - 3.0 * phi) * kf);
        let _c = 0.5 * kf * kp;
        // Wait — use the closed-form expression:
        // k_eff = 1/4 * [(3φ-1)k_p + (2-3φ)k_f +
        //          sqrt(((3φ-1)k_p + (2-3φ)k_f)² + 8 k_p k_f)]
        let disc = (3.0 * phi - 1.0) * kp + (2.0 - 3.0 * phi) * kf;
        0.25 * (disc + (disc * disc + 8.0 * kp * kf).sqrt())
    }

    /// Thermal enhancement ratio k_eff / k_fluid (dimensionless).
    pub fn enhancement_ratio_maxwell(&self) -> f64 {
        self.maxwell_model() / self.k_fluid
    }

    /// Bruggeman enhancement ratio.
    pub fn enhancement_ratio_bruggeman(&self) -> f64 {
        self.bruggeman_model() / self.k_fluid
    }

    /// Brownian motion contribution to thermal conductivity (W m⁻¹ K⁻¹).
    ///
    /// Uses Jang-Choi model: Δk_Brown = β k_f Re_p² Pr φ
    pub fn brownian_contribution(
        &self,
        temperature: f64,
        viscosity: f64,
        specific_heat_fluid: f64,
    ) -> f64 {
        let d_b = K_B * temperature / (3.0 * PI * viscosity * 2.0 * self.particle_radius);
        let re_sq = (d_b / (2.0 * self.particle_radius)).powi(2);
        let pr = viscosity * specific_heat_fluid / self.k_fluid;
        let beta = 7.647e4; // empirical
        beta * self.k_fluid * re_sq * pr * self.volume_fraction
    }
}

// ---------------------------------------------------------------------------
// SlipFlowModel
// ---------------------------------------------------------------------------

/// Navier-slip and Maxwell slip boundary condition model.
#[derive(Debug, Clone)]
pub struct SlipFlowModel {
    /// Mean free path (m).
    pub mean_free_path: f64,
    /// Tangential momentum accommodation coefficient (TMAC) σ_v (0–1).
    pub accommodation_coefficient: f64,
    /// Temperature jump accommodation coefficient σ_T.
    pub temperature_accommodation: f64,
    /// Specific heat ratio γ.
    pub specific_heat_ratio: f64,
    /// Prandtl number.
    pub prandtl_number: f64,
}

impl SlipFlowModel {
    /// Create a new slip flow model.
    pub fn new(
        mean_free_path: f64,
        accommodation_coefficient: f64,
        temperature_accommodation: f64,
        specific_heat_ratio: f64,
        prandtl_number: f64,
    ) -> Self {
        Self {
            mean_free_path,
            accommodation_coefficient,
            temperature_accommodation,
            specific_heat_ratio,
            prandtl_number,
        }
    }

    /// Maxwell slip velocity at the wall (m s⁻¹).
    ///
    /// u_slip = ((2 − σ_v)/σ_v) λ ∂u/∂y|_wall
    ///        + (3/4) μ/(ρ T) ∂T/∂x|_wall  (thermal creep, second term)
    pub fn velocity_slip(&self, velocity_gradient: f64) -> f64 {
        (2.0 - self.accommodation_coefficient) / self.accommodation_coefficient
            * self.mean_free_path
            * velocity_gradient
    }

    /// Temperature jump at the wall (K).
    ///
    /// T_jump = ((2 − σ_T)/σ_T) * 2γ/(γ+1) * λ/Pr * ∂T/∂y|_wall
    pub fn temperature_jump(&self, temperature_gradient: f64) -> f64 {
        let gamma = self.specific_heat_ratio;
        (2.0 - self.temperature_accommodation) / self.temperature_accommodation * 2.0 * gamma
            / (gamma + 1.0)
            * self.mean_free_path
            / self.prandtl_number
            * temperature_gradient
    }

    /// Slip length β (m) = μ_eff / viscosity = ((2-σ_v)/σ_v) λ.
    pub fn slip_length(&self) -> f64 {
        (2.0 - self.accommodation_coefficient) / self.accommodation_coefficient
            * self.mean_free_path
    }

    /// Poiseuille flow correction factor for slip in a channel of half-width h.
    ///
    /// Q_slip/Q_no_slip = 1 + 6 β / h
    pub fn poiseuille_correction(&self, half_width: f64) -> f64 {
        1.0 + 6.0 * self.slip_length() / half_width
    }

    /// Knudsen number for a given characteristic length (dimensionless).
    pub fn knudsen_number(&self, characteristic_length: f64) -> f64 {
        self.mean_free_path / characteristic_length
    }
}

// ---------------------------------------------------------------------------
// DissipativeParticleNano
// ---------------------------------------------------------------------------

/// DPD (dissipative particle dynamics) force parameter set.
#[derive(Debug, Clone)]
pub struct DpdParameters {
    /// Conservative force amplitude A_ij (N).
    pub conservative_amplitude: f64,
    /// Dissipative force coefficient γ_dpd (kg s⁻¹).
    pub dissipative_coefficient: f64,
    /// Random force amplitude σ_dpd (related to γ via FDT).
    pub random_amplitude: f64,
    /// Cut-off radius r_c (m).
    pub cutoff_radius: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl DpdParameters {
    /// Create DPD parameters; random amplitude enforced by FDT: σ² = 2 γ k_B T.
    pub fn new(
        conservative_amplitude: f64,
        dissipative_coefficient: f64,
        cutoff_radius: f64,
        temperature: f64,
    ) -> Self {
        let sigma = (2.0 * dissipative_coefficient * K_B * temperature).sqrt();
        Self {
            conservative_amplitude,
            dissipative_coefficient,
            random_amplitude: sigma,
            cutoff_radius,
            temperature,
        }
    }

    /// DPD weight function w_R(r) = 1 − r/r_c for r < r_c, else 0.
    pub fn weight(&self, r: f64) -> f64 {
        if r < self.cutoff_radius {
            1.0 - r / self.cutoff_radius
        } else {
            0.0
        }
    }

    /// Conservative force (scalar, along e_ij) between two particles at distance r.
    pub fn conservative_force(&self, r: f64) -> f64 {
        self.conservative_amplitude * self.weight(r)
    }

    /// Dissipative force coefficient at distance r.
    pub fn dissipative_force_coeff(&self, r: f64) -> f64 {
        let w = self.weight(r);
        -self.dissipative_coefficient * w * w
    }

    /// Random force amplitude at distance r per sqrt(dt).
    pub fn random_force_amplitude(&self, r: f64) -> f64 {
        let w = self.weight(r);
        self.random_amplitude * w
    }
}

/// DPD-SPH hybrid particle for nanoscale simulations.
#[derive(Debug, Clone)]
pub struct DissipativeParticleNano {
    /// Particle positions (m).
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities (m s⁻¹).
    pub velocities: Vec<[f64; 3]>,
    /// Particle masses (kg).
    pub masses: Vec<f64>,
    /// DPD force parameters.
    pub dpd_params: DpdParameters,
    /// SPH smoothing length (m).
    pub h_sph: f64,
    /// Total force on each particle (N).
    pub forces: Vec<[f64; 3]>,
}

impl DissipativeParticleNano {
    /// Create a new DPD-SPH hybrid solver.
    pub fn new(
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
        dpd_params: DpdParameters,
        h_sph: f64,
    ) -> Self {
        let n = positions.len();
        Self {
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            dpd_params,
            h_sph,
            forces: vec![[0.0; 3]; n],
        }
    }

    /// Number of particles.
    pub fn n_particles(&self) -> usize {
        self.positions.len()
    }

    /// Compute all DPD pairwise forces and accumulate into `forces`.
    pub fn compute_dpd_forces(&mut self, dt: f64) {
        let n = self.n_particles();
        self.forces = vec![[0.0; 3]; n];
        let mut rng = rand::rng();

        for i in 0..n {
            for j in (i + 1)..n {
                let r_vec = sub3(self.positions[j], self.positions[i]);
                let r = norm3(r_vec);
                if r < 1.0e-15 || r >= self.dpd_params.cutoff_radius {
                    continue;
                }
                let e_ij = scale3(r_vec, 1.0 / r);

                // Conservative
                let f_c = self.dpd_params.conservative_force(r);

                // Dissipative
                let v_ij = sub3(self.velocities[j], self.velocities[i]);
                let vdot = dot3(v_ij, e_ij);
                let gamma_w2 = self.dpd_params.dissipative_force_coeff(r);
                let f_d = gamma_w2 * vdot;

                // Random (Gaussian)
                let u1: f64 = rng.random_range(1.0e-15..1.0);
                let u2: f64 = rng.random_range(0.0..2.0 * PI);
                let xi = (-2.0 * u1.ln()).sqrt() * u2.cos();
                let sigma_r = self.dpd_params.random_force_amplitude(r);
                let f_r = sigma_r * xi / dt.sqrt();

                let f_total = f_c + f_d + f_r;
                let f_vec = scale3(e_ij, f_total);
                self.forces[i] = add3(self.forces[i], f_vec);
                self.forces[j] = sub3(self.forces[j], f_vec);
            }
        }
    }

    /// Advance particle positions and velocities by dt using velocity-Verlet.
    pub fn integrate(&mut self, dt: f64) {
        let n = self.n_particles();
        for i in 0..n {
            let inv_m = 1.0 / self.masses[i];
            let acc = scale3(self.forces[i], inv_m);
            self.velocities[i] = add3(self.velocities[i], scale3(acc, 0.5 * dt));
            self.positions[i] = add3(self.positions[i], scale3(self.velocities[i], dt));
        }
    }
}

// ---------------------------------------------------------------------------
// NonequilibriumNano
// ---------------------------------------------------------------------------

/// NEMD-SPH (non-equilibrium molecular dynamics coupled to SPH) model.
///
/// Provides Green-Kubo estimators and non-equilibrium transport coefficients.
#[derive(Debug, Clone)]
pub struct NonequilibriumNano {
    /// Number of particles.
    pub n_particles: usize,
    /// Particle velocities (m s⁻¹).
    pub velocities: Vec<[f64; 3]>,
    /// Particle temperatures (K).
    pub temperatures: Vec<f64>,
    /// Particle masses (kg).
    pub masses: Vec<f64>,
    /// Collected velocity autocorrelation time series.
    pub vacf: Vec<f64>,
    /// Collected heat flux autocorrelation time series.
    pub hfacf: Vec<f64>,
    /// Heat flux vector at current step (W m⁻²).
    pub heat_flux: [f64; 3],
    /// Shear stress at current step (Pa).
    pub shear_stress: f64,
}

impl NonequilibriumNano {
    /// Create a new NEMD-SPH model.
    pub fn new(n_particles: usize, particle_mass: f64) -> Self {
        Self {
            n_particles,
            velocities: vec![[0.0; 3]; n_particles],
            temperatures: vec![T_REF; n_particles],
            masses: vec![particle_mass; n_particles],
            vacf: Vec::new(),
            hfacf: Vec::new(),
            heat_flux: [0.0; 3],
            shear_stress: 0.0,
        }
    }

    /// Compute the velocity autocorrelation function (VACF) value at lag 0.
    pub fn vacf_zero(&self) -> f64 {
        if self.velocities.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.velocities.iter().map(|&v| dot3(v, v)).sum();
        sum / self.n_particles as f64
    }

    /// Green-Kubo diffusion coefficient (m² s⁻¹) from VACF integral.
    ///
    /// D = (1/3) ∫₀^∞ ⟨v(0)·v(t)⟩ dt  ≈ (1/3) Σ VACF(t_i) dt
    pub fn diffusion_coefficient_green_kubo(&self, dt: f64) -> f64 {
        let integral: f64 = self.vacf.iter().sum::<f64>() * dt;
        integral / 3.0
    }

    /// Green-Kubo thermal conductivity (W m⁻¹ K⁻¹) from heat flux ACF.
    ///
    /// k = V/(k_B T²) ∫₀^∞ ⟨J(0)·J(t)⟩ dt
    pub fn thermal_conductivity_green_kubo(&self, dt: f64, volume: f64, temperature: f64) -> f64 {
        let integral: f64 = self.hfacf.iter().sum::<f64>() * dt;
        volume / (K_B * temperature * temperature) * integral
    }

    /// Compute instantaneous heat flux vector (W m⁻²) from kinetic energies and velocities.
    pub fn compute_heat_flux(&mut self) {
        let mut jx = 0.0;
        let mut jy = 0.0;
        let mut jz = 0.0;
        for i in 0..self.n_particles {
            let v = self.velocities[i];
            let ke = 0.5 * self.masses[i] * dot3(v, v);
            jx += ke * v[0];
            jy += ke * v[1];
            jz += ke * v[2];
        }
        self.heat_flux = [jx, jy, jz];
        self.hfacf.push(dot3(self.heat_flux, self.heat_flux));
    }

    /// Compute shear stress from velocity gradient (Pa).
    ///
    /// τ = −μ * (∂u_x/∂y)
    pub fn compute_shear_stress(&mut self, viscosity: f64, velocity_gradient: f64) {
        self.shear_stress = -viscosity * velocity_gradient;
    }

    /// Non-equilibrium viscosity from imposed shear stress and velocity gradient (Pa s).
    pub fn nemd_viscosity(&self, velocity_gradient: f64) -> f64 {
        if velocity_gradient.abs() < 1.0e-30 {
            return 0.0;
        }
        (-self.shear_stress / velocity_gradient).abs()
    }

    /// Mean temperature (K) of all particles.
    pub fn mean_temperature(&self) -> f64 {
        self.temperatures.iter().sum::<f64>() / self.n_particles as f64
    }

    /// Total kinetic energy (J).
    pub fn total_kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(&v, &m)| 0.5 * m * dot3(v, v))
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- NanoscaleSph tests ---

    #[test]
    fn test_knudsen_number() {
        let p = NanoscaleParticle::new([0.0; 3], 1.0e-20, 10.0e-9);
        let sph = NanoscaleSph::new(vec![p], 66.0e-9, 1.0e-6, 1.0e-3, 0.6, 1000.0);
        let kn = sph.knudsen_number();
        assert!((kn - 0.066).abs() < 1.0e-10);
    }

    #[test]
    fn test_flow_regime_slip() {
        let p = NanoscaleParticle::new([0.0; 3], 1.0e-20, 10.0e-9);
        let sph = NanoscaleSph::new(vec![p], 66.0e-9, 1.0e-6, 1.0e-3, 0.6, 1000.0);
        assert_eq!(sph.flow_regime(), KnudsenRegime::SlipFlow);
    }

    #[test]
    fn test_flow_regime_continuum() {
        let p = NanoscaleParticle::new([0.0; 3], 1.0e-20, 10.0e-9);
        let sph = NanoscaleSph::new(vec![p], 1.0e-9, 1.0e-3, 1.0e-3, 0.6, 1000.0);
        assert_eq!(sph.flow_regime(), KnudsenRegime::Continuum);
    }

    #[test]
    fn test_effective_viscosity_nonzero_positive() {
        // Beskok-Karniadakis with b = -1: μ_eff = μ₀/(1 + Kn), which is less than μ₀ for Kn > 0.
        // The formula is correct and the result must be positive.
        let p = NanoscaleParticle::new([0.0; 3], 1.0e-20, 10.0e-9);
        let sph = NanoscaleSph::new(vec![p], 66.0e-9, 1.0e-6, 1.0e-3, 0.6, 1000.0);
        let mu_eff = sph.effective_viscosity();
        assert!(mu_eff > 0.0 && mu_eff < sph.viscosity * 10.0);
    }

    #[test]
    fn test_slip_length_positive() {
        let p = NanoscaleParticle::new([0.0; 3], 1.0e-20, 10.0e-9);
        let sph = NanoscaleSph::new(vec![p], 66.0e-9, 1.0e-6, 1.0e-3, 0.6, 1000.0);
        let beta = sph.slip_length(0.9);
        assert!(beta > 0.0);
    }

    #[test]
    fn test_particle_kinetic_energy() {
        let mut p = NanoscaleParticle::new([0.0; 3], 1.0e-20, 10.0e-9);
        p.velocity = [1.0, 0.0, 0.0];
        let ke = p.kinetic_energy();
        assert!((ke - 0.5e-20).abs() < 1.0e-35);
    }

    // --- SurfaceTensionNano tests ---

    #[test]
    fn test_contact_angle_hydrophilic() {
        // γ_SV > γ_SL → cos θ > 0 → θ < 90°
        let st = SurfaceTensionNano::new(0.072, 0.02, 0.08, 1.0e-20, 1.0e-9);
        let theta = st.contact_angle();
        assert!(theta < PI / 2.0);
    }

    #[test]
    fn test_work_of_adhesion_positive() {
        let st = SurfaceTensionNano::new(0.072, 0.02, 0.08, 1.0e-20, 1.0e-9);
        assert!(st.work_of_adhesion() > 0.0);
    }

    #[test]
    fn test_disjoining_pressure_vdw_negative() {
        let st = SurfaceTensionNano::new(0.072, 0.02, 0.08, 1.0e-20, 5.0e-9);
        let pi = st.disjoining_pressure_vdw(5.0e-9);
        assert!(
            pi < 0.0,
            "vdW disjoining pressure should be attractive (negative) for A>0"
        );
    }

    #[test]
    fn test_laplace_pressure_nanodroplet() {
        let st = SurfaceTensionNano::new(0.072, 0.02, 0.08, 1.0e-20, 1.0e-9);
        let r = 50.0e-9;
        let lp = st.laplace_pressure(r);
        let expected = 2.0 * 0.072 / r;
        assert!((lp - expected).abs() < 1.0e-3);
    }

    // --- ThermalFluctuationSph tests ---

    #[test]
    fn test_thermal_fluctuation_noise_amplitude_positive() {
        let tf = ThermalFluctuationSph::new(100, 1.0e-20, 1.0e-3, 0.6, 298.15, 1.0e-24);
        let sigma = tf.noise_amplitude(1.0e-12);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_stokes_einstein_diffusion() {
        let tf = ThermalFluctuationSph::new(100, 1.0e-20, 1.0e-3, 0.6, 298.15, 1.0e-24);
        let d = tf.stokes_einstein_diffusion(10.0e-9);
        let expected = K_B * 298.15 / (3.0 * PI * 1.0e-3 * 10.0e-9);
        assert!((d - expected).abs() / expected < 1.0e-10);
    }

    #[test]
    fn test_thermal_fluctuation_generate_forces() {
        let mut tf = ThermalFluctuationSph::new(10, 1.0e-20, 1.0e-3, 0.6, 298.15, 1.0e-24);
        tf.generate_random_forces(1.0e-12);
        assert_eq!(tf.random_forces.len(), 10);
    }

    #[test]
    fn test_mean_square_displacement_positive() {
        let tf = ThermalFluctuationSph::new(100, 1.0e-20, 1.0e-3, 0.6, 298.15, 1.0e-24);
        let msd = tf.mean_square_displacement(1.0e-9, 10.0e-9);
        assert!(msd > 0.0);
    }

    // --- ElectricDoubleLayer tests ---

    #[test]
    fn test_debye_length_order_of_magnitude() {
        // 1 mM NaCl → λ_D ≈ 9.6 nm
        let edl = ElectricDoubleLayer::new(1.0, 1, -0.05, 298.15, EPS_WATER);
        let lambda = edl.debye_length();
        // Should be in the nm range (1e-9 to 100e-9)
        assert!(lambda > 1.0e-9 && lambda < 1.0e-7, "λ_D = {:.3e} m", lambda);
    }

    #[test]
    fn test_potential_profile_decays() {
        let edl = ElectricDoubleLayer::new(1.0, 1, -0.05, 298.15, EPS_WATER);
        let lambda = edl.debye_length();
        let psi0 = edl.potential_profile(0.0);
        let psi1 = edl.potential_profile(lambda);
        assert!(psi1.abs() < psi0.abs());
    }

    #[test]
    fn test_electro_osmotic_velocity_sign() {
        let edl = ElectricDoubleLayer::new(1.0, 1, -0.05, 298.15, EPS_WATER);
        let v = edl.electro_osmotic_velocity(1000.0, 1.0e-3);
        // Negative surface potential → positive EOF velocity for positive E field
        assert!(v > 0.0, "Expected positive EOF velocity, got {:.4e}", v);
    }

    #[test]
    fn test_edl_overlap_parameter() {
        let edl = ElectricDoubleLayer::new(1.0, 1, -0.05, 298.15, EPS_WATER);
        let lambda = edl.debye_length();
        let kappa_h = edl.edl_overlap_parameter(lambda);
        assert!((kappa_h - 1.0).abs() < 1.0e-10);
    }

    // --- MolecularConfinement tests ---

    #[test]
    fn test_molecular_confinement_layers() {
        let mc = MolecularConfinement::new(5.0e-9, 0.3e-9, 1000.0, 1.0e-21, 0.3e-9);
        let n = mc.n_molecular_layers();
        assert!((n - 10.0 / 0.3).abs() < 1.0e-3, "n_layers = {:.3}", n);
    }

    #[test]
    fn test_flow_enhancement_factor_greater_than_one() {
        let mc = MolecularConfinement::new(1.0e-9, 0.3e-9, 1000.0, 1.0e-21, 0.3e-9);
        let ef = mc.flow_enhancement_factor();
        assert!(ef > 1.0);
    }

    #[test]
    fn test_density_profile_bulk_far_from_wall() {
        let mc = MolecularConfinement::new(5.0e-9, 0.3e-9, 1000.0, 1.0e-21, 0.3e-9);
        // Very far from wall → should approach bulk density
        let rho = mc.density_profile(50.0e-9);
        assert!((rho - 1000.0).abs() < 1.0, "rho = {:.4}", rho);
    }

    // --- NanobubbleSph tests ---

    #[test]
    fn test_nanobubble_internal_pressure_exceeds_atm() {
        let nb = NanobubbleSph::new(NanobubbleType::Bulk, 100.0e-9, 0, 0.072, 0.0, 7.0e-9, 0.0);
        let p_int = nb.internal_pressure(101_325.0);
        assert!(p_int > 101_325.0);
    }

    #[test]
    fn test_nanobubble_laplace_pressure_large_r() {
        // Larger bubble → smaller Laplace pressure
        let nb1 = NanobubbleSph::new(NanobubbleType::Bulk, 50.0e-9, 0, 0.072, 0.0, 7.0e-9, 0.0);
        let nb2 = NanobubbleSph::new(NanobubbleType::Bulk, 500.0e-9, 0, 0.072, 0.0, 7.0e-9, 0.0);
        assert!(nb1.internal_pressure(0.0) > nb2.internal_pressure(0.0));
    }

    #[test]
    fn test_nanobubble_dissolution_rate_sign() {
        // If c_s > c_inf, bubble dissolves → dr/dt < 0
        let nb = NanobubbleSph::new(NanobubbleType::Bulk, 100.0e-9, 0, 0.072, 0.0, 7.0e-9, 0.0);
        let rate = nb.dissolution_rate(2.0e-9, 0.024, 101_325.0);
        assert!(rate < 0.0);
    }

    // --- NanofluidThermal tests ---

    #[test]
    fn test_maxwell_model_greater_than_base() {
        let nf = NanofluidThermal::new(0.6, 400.0, 0.01, 20.0e-9, 1.0e-8);
        let k_eff = nf.maxwell_model();
        assert!(k_eff > nf.k_fluid);
    }

    #[test]
    fn test_bruggeman_model_greater_than_base() {
        let nf = NanofluidThermal::new(0.6, 400.0, 0.01, 20.0e-9, 1.0e-8);
        let k_eff = nf.bruggeman_model();
        assert!(k_eff > nf.k_fluid);
    }

    #[test]
    fn test_maxwell_bruggeman_agree_low_phi() {
        // For small φ, Maxwell and Bruggeman should give similar results
        let nf = NanofluidThermal::new(0.6, 400.0, 0.001, 20.0e-9, 1.0e-8);
        let km = nf.maxwell_model();
        let kb = nf.bruggeman_model();
        let rel_diff = (km - kb).abs() / km;
        assert!(rel_diff < 0.01, "rel diff = {:.4}", rel_diff);
    }

    #[test]
    fn test_kapitza_resistance_reduces_conductivity() {
        let nf_no_kapitza = NanofluidThermal::new(0.6, 400.0, 0.01, 20.0e-9, 0.0);
        let nf_kapitza = NanofluidThermal::new(0.6, 400.0, 0.01, 20.0e-9, 1.0e-8);
        assert!(
            nf_no_kapitza.maxwell_model() >= nf_kapitza.maxwell_model(),
            "Kapitza resistance should reduce effective conductivity"
        );
    }

    // --- SlipFlowModel tests ---

    #[test]
    fn test_slip_length_formula() {
        let sfm = SlipFlowModel::new(66.0e-9, 0.9, 0.9, 1.4, 0.71);
        let beta = sfm.slip_length();
        let expected = (2.0 - 0.9) / 0.9 * 66.0e-9;
        assert!((beta - expected).abs() < 1.0e-20);
    }

    #[test]
    fn test_poiseuille_correction_greater_than_one() {
        let sfm = SlipFlowModel::new(66.0e-9, 0.9, 0.9, 1.4, 0.71);
        let corr = sfm.poiseuille_correction(500.0e-9);
        assert!(corr > 1.0);
    }

    #[test]
    fn test_velocity_slip_proportional_to_gradient() {
        let sfm = SlipFlowModel::new(66.0e-9, 0.9, 0.9, 1.4, 0.71);
        let u1 = sfm.velocity_slip(1.0e6);
        let u2 = sfm.velocity_slip(2.0e6);
        assert!((u2 - 2.0 * u1).abs() < 1.0e-25);
    }

    // --- DissipativeParticleNano tests ---

    #[test]
    fn test_dpd_weight_at_cutoff() {
        let params = DpdParameters::new(25.0, 4.5, 1.0e-9, 298.15);
        assert!((params.weight(1.0e-9) - 0.0).abs() < 1.0e-15);
    }

    #[test]
    fn test_dpd_weight_at_zero() {
        let params = DpdParameters::new(25.0, 4.5, 1.0e-9, 298.15);
        assert!((params.weight(0.0) - 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn test_dpd_fdt_relation() {
        let params = DpdParameters::new(25.0, 4.5, 1.0e-9, 298.15);
        let sigma_sq = params.random_amplitude * params.random_amplitude;
        let expected = 2.0 * params.dissipative_coefficient * K_B * params.temperature;
        assert!((sigma_sq - expected).abs() / expected < 1.0e-10);
    }

    #[test]
    fn test_dpd_sph_n_particles() {
        let positions = vec![[0.0; 3], [1.0e-9, 0.0, 0.0], [2.0e-9, 0.0, 0.0]];
        let masses = vec![1.0e-24; 3];
        let params = DpdParameters::new(25.0, 4.5, 3.0e-9, 298.15);
        let dpd = DissipativeParticleNano::new(positions, masses, params, 1.5e-9);
        assert_eq!(dpd.n_particles(), 3);
    }

    // --- NonequilibriumNano tests ---

    #[test]
    fn test_vacf_zero_from_velocities() {
        let mut nemd = NonequilibriumNano::new(3, 1.0e-24);
        nemd.velocities = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let vacf0 = nemd.vacf_zero();
        assert!((vacf0 - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_total_kinetic_energy_positive() {
        let mut nemd = NonequilibriumNano::new(2, 1.0e-24);
        nemd.velocities = vec![[100.0, 0.0, 0.0], [0.0, 100.0, 0.0]];
        let ke = nemd.total_kinetic_energy();
        assert!(ke > 0.0);
    }

    #[test]
    fn test_nemd_viscosity_from_shear_stress() {
        let mut nemd = NonequilibriumNano::new(10, 1.0e-24);
        nemd.compute_shear_stress(1.0e-3, 1.0e6);
        let mu = nemd.nemd_viscosity(1.0e6);
        assert!((mu - 1.0e-3).abs() < 1.0e-15);
    }

    #[test]
    fn test_nemd_mean_temperature() {
        let mut nemd = NonequilibriumNano::new(4, 1.0e-24);
        nemd.temperatures = vec![300.0, 310.0, 290.0, 300.0];
        let t_mean = nemd.mean_temperature();
        assert!((t_mean - 300.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_compute_heat_flux_stores_hfacf() {
        let mut nemd = NonequilibriumNano::new(3, 1.0e-24);
        nemd.velocities = vec![[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        nemd.compute_heat_flux();
        assert_eq!(nemd.hfacf.len(), 1);
    }
}
