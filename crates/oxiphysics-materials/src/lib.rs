// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Material properties and material library.
//!
//! Defines physical material properties such as density, friction,
//! and restitution, along with a library for managing named materials.
#![warn(missing_docs)]

mod error;
pub use error::*;

pub mod acoustics;
pub mod additive_manufacturing;
pub mod aerospace;
pub mod anisotropic;
pub mod battery_materials;
pub mod biological_materials;
pub mod biomaterials;
pub mod biomechanics;
pub mod ceramic_materials;
pub mod combination;
pub mod composite;
pub mod composite_failure;
pub mod composite_materials;
pub mod constitutive;
pub mod construction;
pub mod creep;
pub mod crystal_plasticity;
pub mod damage;
pub mod dielectric_materials;
pub mod elastic;
pub mod electrochemistry;
pub mod energy_materials;
pub mod eos;
pub mod fatigue;
pub mod fiber_composites;
pub mod fracture;
pub mod geological;
pub mod geomechanics;
pub mod homogenization;
pub mod hydrogen_storage;
pub mod hyperelastic;
pub mod metamaterials;
pub mod multiphysics;
pub mod nano;
pub mod nano_materials;
pub mod nanocomposites;
pub mod nanomaterials;
pub mod nuclear_materials;
pub mod optical;
pub mod optical_materials;
pub mod phase_transform;
pub mod plasticity;
pub mod polymer_mechanics;
pub mod polymer_physics;
pub mod porous_media;
pub mod presets;
pub mod quantum_materials;
pub mod radiation;
pub mod radiation_shielding;
pub mod semiconductor;
pub mod shape_memory;
pub mod smart_materials;
pub mod superconductor;
pub mod thermal;
pub mod thermoelectrics;
pub mod tribology;
pub mod tribology_ext;
pub mod viscoelastic;

pub mod aerospace_materials;

pub mod corrosion;

pub mod additive_manufacturing_materials;
pub mod foam_materials;
pub mod magnetocaloric_materials;

pub mod materials_bench;
pub mod simd_paths;
pub use simd_paths::{
    MaterialBatchStats, SoaMaterialPoints, elastic_stress_batch, elastic_stress_batch_x4, lame,
    miner_damage_batch, neo_hookean_stress_batch, return_mapping_batch,
    thermal_expansion_stress_batch, viscoplastic_rate_batch, von_mises_yield_batch,
};

pub use combination::{
    ContactMaterialPair, FrictionCombineRule, ModulusCombineRule, RestitutionCombineRule,
    ThermalContactResistance, combine_friction, combine_modulus, combine_restitution,
    hertz_contact_force, hertz_effective_modulus, hertz_effective_radius, maxwell_diffusivity,
};
pub use composite::{
    HalpinTsai, HashinResult, IsotropicConstituent, Laminate, MoriTanaka, OrthotropicPly,
    PlyFailureResult, PlyStrength, PlyStress, ProgressiveFailureResult, PuckResult, hashin_failure,
    hashin_shtrikman_lower, hashin_shtrikman_upper, interlaminar_shear_estimate,
    max_stress_failure, ply_by_ply_stress, progressive_failure_analysis, reuss_modulus,
    thermal_residual_stresses, tsai_wu_failure, voigt_density, voigt_modulus, voigt_poisson,
};
pub use constitutive::{ConstitutiveModel, ConstitutiveResponse};
pub use creep::functions::*;
pub use creep::types::*;
pub use elastic::{LinearElastic, NeoHookean};
pub use eos::{
    EosWithEnergy, EquationOfState, IdealGasEos, MieGruneisenEos as MieGruneisenEosShock,
    PolynomialEos, StiffenedGasEos, TaitEos, TillotsonEos, VanDerWaalsEos,
};
pub use hyperelastic::{DruckerPrager, J2Plasticity, JwlEos, MieGruneisenEos, MooneyRivlin, Ogden};
pub use phase_transform::functions::*;
pub use phase_transform::types::*;
pub use viscoelastic::{KelvinVoigt, Maxwell, StandardLinearSolid};

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Physical material properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Material {
    /// Material name
    pub name: String,
    /// Density in kg/m^3
    pub density: f64,
    /// Coefficient of friction (0..1)
    pub friction: f64,
    /// Coefficient of restitution (bounciness, 0..1)
    pub restitution: f64,
}

impl Material {
    /// Create a new material with the given properties.
    pub fn new(name: impl Into<String>, density: f64, friction: f64, restitution: f64) -> Self {
        Self {
            name: name.into(),
            density,
            friction,
            restitution,
        }
    }
}

/// A library of named materials.
#[derive(Debug, Default)]
pub struct MaterialLibrary {
    materials: HashMap<String, Material>,
}

impl MaterialLibrary {
    /// Create a new empty material library.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a material to the library.
    pub fn add(&mut self, material: Material) {
        self.materials.insert(material.name.clone(), material);
    }

    /// Look up a material by name.
    pub fn get(&self, name: &str) -> Option<&Material> {
        self.materials.get(name)
    }

    /// Return the number of materials in the library.
    pub fn len(&self) -> usize {
        self.materials.len()
    }

    /// Check whether the library is empty.
    pub fn is_empty(&self) -> bool {
        self.materials.is_empty()
    }

    /// Remove a material by name. Returns the removed material if it existed.
    pub fn remove(&mut self, name: &str) -> Option<Material> {
        self.materials.remove(name)
    }

    /// Return an iterator over all materials.
    pub fn iter(&self) -> impl Iterator<Item = (&String, &Material)> {
        self.materials.iter()
    }

    /// Check if a material with the given name exists.
    pub fn contains(&self, name: &str) -> bool {
        self.materials.contains_key(name)
    }

    /// Get all material names.
    pub fn names(&self) -> Vec<String> {
        self.materials.keys().cloned().collect()
    }
}

impl Material {
    /// Compute the impedance (density * speed_of_sound).
    ///
    /// `speed_of_sound` - Speed of sound in the material (m/s).
    pub fn impedance(&self, speed_of_sound: f64) -> f64 {
        self.density * speed_of_sound
    }

    /// Check whether this material is heavier than another.
    pub fn is_heavier_than(&self, other: &Material) -> bool {
        self.density > other.density
    }

    /// Check whether this material is bouncier than another.
    pub fn is_bouncier_than(&self, other: &Material) -> bool {
        self.restitution > other.restitution
    }

    /// Combine this material with another using average rules.
    pub fn combine_average(&self, other: &Material) -> Material {
        Material {
            name: format!("{}+{}", self.name, other.name),
            density: (self.density + other.density) * 0.5,
            friction: (self.friction + other.friction) * 0.5,
            restitution: (self.restitution + other.restitution) * 0.5,
        }
    }

    /// Compute the weight of a given volume of this material (density * volume * g).
    pub fn weight(&self, volume: f64, gravity: f64) -> f64 {
        self.density * volume * gravity
    }

    /// Compute the mass of a given volume of this material.
    pub fn mass_from_volume(&self, volume: f64) -> f64 {
        self.density * volume
    }
}

/// Linearly interpolate between two material property values.
pub fn lerp_property(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Compute the effective coefficient of restitution for a collision
/// using the velocity-dependent model:
/// e_eff = e * max(0, 1 - k * |v_rel|)
///
/// where `k` controls velocity-dependence (typically small, e.g. 0.01).
pub fn velocity_dependent_restitution(e: f64, v_rel: f64, k: f64) -> f64 {
    (e * (1.0 - k * v_rel.abs())).max(0.0)
}

/// Compute the critical damping coefficient for a spring-mass system.
///
/// c_crit = 2 * sqrt(k * m)
pub fn critical_damping(stiffness: f64, mass: f64) -> f64 {
    2.0 * (stiffness * mass).sqrt()
}

/// Compute the damping ratio: zeta = c / c_crit.
pub fn damping_ratio(damping: f64, stiffness: f64, mass: f64) -> f64 {
    let c_crit = critical_damping(stiffness, mass);
    if c_crit.abs() < 1e-15 {
        0.0
    } else {
        damping / c_crit
    }
}

/// Compute the natural frequency of a spring-mass system: omega_n = sqrt(k/m).
pub fn natural_frequency(stiffness: f64, mass: f64) -> f64 {
    if mass.abs() < 1e-15 {
        0.0
    } else {
        (stiffness / mass).sqrt()
    }
}

// ─── Extended Material with Full Physical Properties ──────────────────────────

/// All scalar properties required to construct an [`ExtendedMaterial`].
///
/// Grouping the 12 physical quantities into a single struct lets callers use
/// named-field syntax instead of a long positional argument list.
#[derive(Debug, Clone)]
pub struct ExtendedMaterialProps {
    /// Density \[kg/m³\]
    pub density: f64,
    /// Young's modulus \[Pa\]
    pub young_modulus: f64,
    /// Poisson's ratio \[-\]
    pub poisson_ratio: f64,
    /// Yield strength \[Pa\]
    pub yield_strength: f64,
    /// Ultimate tensile strength \[Pa\]
    pub ultimate_strength: f64,
    /// Thermal conductivity \[W/(m·K)\]
    pub thermal_conductivity: f64,
    /// Specific heat capacity at constant pressure \[J/(kg·K)\]
    pub specific_heat: f64,
    /// Thermal expansion coefficient \[1/K\]
    pub thermal_expansion: f64,
    /// Melting point \[K\]
    pub melting_point: f64,
    /// Speed of sound \[m/s\]
    pub speed_of_sound: f64,
    /// Electrical resistivity \[Ω·m\]
    pub electrical_resistivity: f64,
    /// Emissivity (0 = perfect reflector, 1 = blackbody)
    pub emissivity: f64,
}

/// Extended material with full physical, mechanical, thermal, acoustic, and
/// optical properties for simulation of complex multi-physics problems.
#[derive(Debug, Clone)]
pub struct ExtendedMaterial {
    /// Material name
    pub name: String,
    /// Density \[kg/m³\]
    pub density: f64,
    /// Young's modulus \[Pa\]
    pub young_modulus: f64,
    /// Poisson's ratio \[-\]
    pub poisson_ratio: f64,
    /// Yield strength \[Pa\]
    pub yield_strength: f64,
    /// Ultimate tensile strength \[Pa\]
    pub ultimate_strength: f64,
    /// Thermal conductivity \[W/(m·K)\]
    pub thermal_conductivity: f64,
    /// Specific heat capacity at constant pressure \[J/(kg·K)\]
    pub specific_heat: f64,
    /// Thermal expansion coefficient \[1/K\]
    pub thermal_expansion: f64,
    /// Melting point \[K\]
    pub melting_point: f64,
    /// Speed of sound \[m/s\]
    pub speed_of_sound: f64,
    /// Electrical resistivity \[Ω·m\]
    pub electrical_resistivity: f64,
    /// Emissivity (0 = perfect reflector, 1 = blackbody)
    pub emissivity: f64,
}

impl ExtendedMaterial {
    /// Create a new extended material from a name and an [`ExtendedMaterialProps`]
    /// parameter struct.
    pub fn new(name: impl Into<String>, props: ExtendedMaterialProps) -> Self {
        Self {
            name: name.into(),
            density: props.density,
            young_modulus: props.young_modulus,
            poisson_ratio: props.poisson_ratio,
            yield_strength: props.yield_strength,
            ultimate_strength: props.ultimate_strength,
            thermal_conductivity: props.thermal_conductivity,
            specific_heat: props.specific_heat,
            thermal_expansion: props.thermal_expansion,
            melting_point: props.melting_point,
            speed_of_sound: props.speed_of_sound,
            electrical_resistivity: props.electrical_resistivity,
            emissivity: props.emissivity,
        }
    }

    /// Shear modulus G = E / (2(1+ν)) \[Pa\].
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Bulk modulus K = E / (3(1-2ν)) \[Pa\].
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Thermal diffusivity α = λ / (ρ·cp) \[m²/s\].
    pub fn thermal_diffusivity(&self) -> f64 {
        self.thermal_conductivity / (self.density * self.specific_heat)
    }

    /// Acoustic impedance Z = ρ·c \[kg/(m²·s) = Pa·s/m\].
    pub fn acoustic_impedance(&self) -> f64 {
        self.density * self.speed_of_sound
    }

    /// Specific strength (yield strength / density) \[Pa·m³/kg = N·m/kg\].
    pub fn specific_strength(&self) -> f64 {
        self.yield_strength / self.density
    }

    /// Specific stiffness (Young's modulus / density) \[Pa·m³/kg\].
    pub fn specific_stiffness(&self) -> f64 {
        self.young_modulus / self.density
    }

    /// Thermal shock resistance (simplified Biot-type):
    /// R = σ_y · λ / (E · α) \[W/m\].
    pub fn thermal_shock_resistance(&self) -> f64 {
        self.yield_strength * self.thermal_conductivity
            / (self.young_modulus * self.thermal_expansion)
    }

    /// Reflection coefficient at normal incidence from medium 1 (air, Z₁ = 415)
    /// to this material: R = ((Z₂ - Z₁) / (Z₂ + Z₁))².
    pub fn acoustic_reflection_coefficient_from_air(&self) -> f64 {
        const Z_AIR: f64 = 415.0; // kg/(m²·s) at 20°C
        let z2 = self.acoustic_impedance();
        let r = (z2 - Z_AIR) / (z2 + Z_AIR);
        r * r
    }

    /// Skin depth at frequency f \[Hz\] for electromagnetic waves:
    /// δ = sqrt(ρ_e / (π · f · μ₀)) \[m\].
    pub fn skin_depth(&self, frequency_hz: f64) -> f64 {
        const MU_0: f64 = 1.256_637_061_4e-6; // H/m
        (self.electrical_resistivity / (std::f64::consts::PI * frequency_hz * MU_0)).sqrt()
    }

    /// Thermal boundary layer thickness for flow over a flat plate at x:
    /// δ_T = 5x / sqrt(Re_x) · Pr^{-1/3},  Re_x = ρ·v·x/η.
    ///
    /// Returns δ_T \[m\].
    pub fn thermal_boundary_layer(
        &self,
        x: f64,
        velocity: f64,
        dynamic_viscosity: f64,
        prandtl: f64,
    ) -> f64 {
        let re_x = self.density * velocity * x / dynamic_viscosity;
        if re_x < 1e-15 {
            return 0.0;
        }
        5.0 * x / re_x.sqrt() * prandtl.powf(-1.0 / 3.0)
    }

    /// Return `true` if this material has higher specific stiffness than another.
    pub fn is_stiffer_specific(&self, other: &ExtendedMaterial) -> bool {
        self.specific_stiffness() > other.specific_stiffness()
    }

    /// Preset: structural steel (S355).
    pub fn steel() -> Self {
        Self::new(
            "steel_s355",
            ExtendedMaterialProps {
                density: 7850.0,
                young_modulus: 200.0e9,
                poisson_ratio: 0.30,
                yield_strength: 355.0e6,
                ultimate_strength: 490.0e6,
                thermal_conductivity: 50.0,
                specific_heat: 486.0,
                thermal_expansion: 12.0e-6,
                melting_point: 1808.0,
                speed_of_sound: 5960.0,
                electrical_resistivity: 1.7e-7,
                emissivity: 0.28,
            },
        )
    }

    /// Preset: aluminium alloy (6061-T6).
    pub fn aluminium_6061() -> Self {
        Self::new(
            "aluminium_6061",
            ExtendedMaterialProps {
                density: 2700.0,
                young_modulus: 68.9e9,
                poisson_ratio: 0.33,
                yield_strength: 276.0e6,
                ultimate_strength: 310.0e6,
                thermal_conductivity: 167.0,
                specific_heat: 896.0,
                thermal_expansion: 23.6e-6,
                melting_point: 933.0,
                speed_of_sound: 6320.0,
                electrical_resistivity: 4.0e-8,
                emissivity: 0.09,
            },
        )
    }

    /// Preset: Ti-6Al-4V titanium alloy.
    pub fn titanium_6al4v() -> Self {
        Self::new(
            "titanium_6al4v",
            ExtendedMaterialProps {
                density: 4430.0,
                young_modulus: 114.0e9,
                poisson_ratio: 0.34,
                yield_strength: 880.0e6,
                ultimate_strength: 950.0e6,
                thermal_conductivity: 6.7,
                specific_heat: 560.0,
                thermal_expansion: 8.6e-6,
                melting_point: 1878.0,
                speed_of_sound: 6070.0,
                electrical_resistivity: 1.7e-6,
                emissivity: 0.35,
            },
        )
    }

    /// Preset: Carbon Fibre Reinforced Polymer (quasi-isotropic laminate).
    pub fn cfrp_quasi_isotropic() -> Self {
        Self::new(
            "cfrp_quasi_iso",
            ExtendedMaterialProps {
                density: 1550.0,
                young_modulus: 70.0e9,
                poisson_ratio: 0.30,
                yield_strength: 600.0e6,
                ultimate_strength: 700.0e6,
                thermal_conductivity: 5.0,
                specific_heat: 800.0,
                thermal_expansion: 2.0e-6,
                melting_point: 3800.0,
                speed_of_sound: 3000.0,
                electrical_resistivity: 1e4,
                emissivity: 0.95,
            },
        )
    }

    /// Preset: borosilicate glass (Pyrex).
    pub fn borosilicate_glass() -> Self {
        Self::new(
            "borosilicate_glass",
            ExtendedMaterialProps {
                density: 2230.0,
                young_modulus: 64.0e9,
                poisson_ratio: 0.20,
                yield_strength: 40.0e6,
                ultimate_strength: 40.0e6,
                thermal_conductivity: 1.2,
                specific_heat: 830.0,
                thermal_expansion: 3.3e-6,
                melting_point: 1100.0,
                speed_of_sound: 5640.0,
                electrical_resistivity: 1e12,
                emissivity: 0.92,
            },
        )
    }
}

// ─── Material Mixture and Interpolation ───────────────────────────────────────

/// Material mixing rules for composite and alloy properties.
pub struct MaterialMixture;

impl MaterialMixture {
    /// Linear (isostrain / Voigt) rule for Young's modulus.
    ///
    /// E_V = φ·E₂ + (1-φ)·E₁
    pub fn voigt_modulus(e1: f64, e2: f64, phi: f64) -> f64 {
        (1.0 - phi) * e1 + phi * e2
    }

    /// Isostress (Reuss) rule for Young's modulus.
    ///
    /// 1/E_R = φ/E₂ + (1-φ)/E₁
    pub fn reuss_modulus(e1: f64, e2: f64, phi: f64) -> f64 {
        let inv = (1.0 - phi) / e1 + phi / e2;
        1.0 / inv
    }

    /// Hill average: E_H = (E_V + E_R) / 2.
    pub fn hill_modulus(e1: f64, e2: f64, phi: f64) -> f64 {
        0.5 * (Self::voigt_modulus(e1, e2, phi) + Self::reuss_modulus(e1, e2, phi))
    }

    /// Voigt (linear) rule for thermal conductivity.
    pub fn voigt_conductivity(k1: f64, k2: f64, phi: f64) -> f64 {
        (1.0 - phi) * k1 + phi * k2
    }

    /// Hashin-Shtrikman upper bound for bulk modulus.
    ///
    /// # Arguments
    /// * `k1`, `k2` — bulk moduli of phases 1 and 2 \[Pa\]
    /// * `g1`, `g2` — shear moduli of phases 1 and 2 \[Pa\]
    /// * `phi`      — volume fraction of phase 2
    pub fn hashin_shtrikman_upper_bulk(k1: f64, k2: f64, g1: f64, _g2: f64, phi: f64) -> f64 {
        // HS+ uses phase 1 (stiffer) as reference if k1 > k2
        let (k_ref, k_inc, _g_ref, phi_inc) = if k1 >= k2 {
            (k1, k2, g1, phi)
        } else {
            (k2, k1, _g2, 1.0 - phi)
        };
        let denom = 3.0 * k_ref + 4.0 * _g_ref;
        k_ref + phi_inc / ((1.0 / (k_inc - k_ref)) + (1.0 - phi_inc) * 3.0 / denom)
    }

    /// Hashin-Shtrikman lower bound for bulk modulus.
    pub fn hashin_shtrikman_lower_bulk(k1: f64, k2: f64, g1: f64, g2: f64, phi: f64) -> f64 {
        // HS- uses the softer phase as reference
        Self::hashin_shtrikman_upper_bulk(k2, k1, g2, g1, 1.0 - phi)
    }

    /// Mixture density: ρ = φ·ρ₂ + (1-φ)·ρ₁ \[kg/m³\].
    pub fn mixture_density(rho1: f64, rho2: f64, phi: f64) -> f64 {
        (1.0 - phi) * rho1 + phi * rho2
    }

    /// Effective thermal expansion coefficient (Turner model):
    ///
    /// α_eff = (φ·α₂·K₂ + (1-φ)·α₁·K₁) / (φ·K₂ + (1-φ)·K₁)
    ///
    /// where K_i are bulk moduli of the individual phases.
    pub fn turner_thermal_expansion(alpha1: f64, k1: f64, alpha2: f64, k2: f64, phi: f64) -> f64 {
        let num = (1.0 - phi) * alpha1 * k1 + phi * alpha2 * k2;
        let den = (1.0 - phi) * k1 + phi * k2;
        if den.abs() < f64::EPSILON {
            0.0
        } else {
            num / den
        }
    }
}

// ─── Material Selection Utilities ─────────────────────────────────────────────

/// Multi-criterion material selection index (Ashby-type).
pub struct MaterialIndex;

impl MaterialIndex {
    /// Minimum mass beam (stiffness-limited): E^{1/2} / ρ.
    ///
    /// Higher is better for minimum weight beams under bending.
    pub fn beam_stiffness_index(young_modulus: f64, density: f64) -> f64 {
        young_modulus.sqrt() / density
    }

    /// Minimum mass column (buckling-limited): E^{1/3} / ρ.
    pub fn column_buckling_index(young_modulus: f64, density: f64) -> f64 {
        young_modulus.powf(1.0 / 3.0) / density
    }

    /// Minimum mass panel (stiffness-limited): E^{1/3} / ρ.
    pub fn panel_stiffness_index(young_modulus: f64, density: f64) -> f64 {
        young_modulus.powf(1.0 / 3.0) / density
    }

    /// Strength-limited beam: σ_y^{2/3} / ρ.
    pub fn beam_strength_index(yield_strength: f64, density: f64) -> f64 {
        yield_strength.powf(2.0 / 3.0) / density
    }

    /// Thermal insulation index (low conductivity, high specific heat):
    /// λ / (ρ · cp).  Lower thermal diffusivity = better insulator.
    pub fn thermal_insulation_index(conductivity: f64, density: f64, specific_heat: f64) -> f64 {
        conductivity / (density * specific_heat)
    }

    /// Pressure vessel (fracture-limited): K_Ic² / (π · σ_y²)
    ///
    /// Proportional to critical crack length — higher is safer.
    pub fn pressure_vessel_index(k_ic: f64, yield_strength: f64) -> f64 {
        k_ic * k_ic / (std::f64::consts::PI * yield_strength * yield_strength)
    }

    /// Wear resistance index: H / (E · λ) where H ≈ 3σ_y (Vickers).
    ///
    /// Higher yield → better wear resistance relative to contact mechanics.
    pub fn wear_resistance_index(yield_strength: f64, young_modulus: f64) -> f64 {
        3.0 * yield_strength / young_modulus
    }
}

// ─── Elastic Moduli Conversions ────────────────────────────────────────────────

/// Convert (E, ν) → (K, G, λ) — the standard elastic moduli family.
///
/// Returns `(K, G, lambda)` where:
/// - K = bulk modulus \[Pa\]
/// - G = shear modulus \[Pa\]
/// - λ = Lamé's first parameter \[Pa\]
pub fn elastic_moduli_from_e_nu(e: f64, nu: f64) -> (f64, f64, f64) {
    let k = e / (3.0 * (1.0 - 2.0 * nu));
    let g = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    (k, g, lambda)
}

/// Convert (K, G) → (E, ν, λ).
///
/// Returns `(E, nu, lambda)`.
pub fn elastic_moduli_from_kg(k: f64, g: f64) -> (f64, f64, f64) {
    let e = 9.0 * k * g / (3.0 * k + g);
    let nu = (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g));
    let lambda = k - 2.0 * g / 3.0;
    (e, nu, lambda)
}

/// Convert (λ, G) → (E, ν, K).
///
/// Returns `(E, nu, K)`.
pub fn elastic_moduli_from_lame(lambda: f64, g: f64) -> (f64, f64, f64) {
    let e = g * (3.0 * lambda + 2.0 * g) / (lambda + g);
    let nu = lambda / (2.0 * (lambda + g));
    let k = lambda + 2.0 * g / 3.0;
    (e, nu, k)
}

/// Check physical admissibility of (E, ν):
/// - E > 0
/// - -1 < ν < 0.5 (for isotropic solids; auxetic materials allow ν < 0)
pub fn elastic_moduli_admissible(e: f64, nu: f64) -> bool {
    e > 0.0 && nu > -1.0 && nu < 0.5
}

// ─── Wave Speed Utilities ──────────────────────────────────────────────────────

/// Longitudinal (P-wave) speed: c_p = sqrt((K + 4G/3) / ρ) \[m/s\].
///
/// This is the speed of dilatational waves in an infinite elastic solid.
pub fn p_wave_speed(bulk_modulus: f64, shear_modulus: f64, density: f64) -> f64 {
    let m = bulk_modulus + 4.0 * shear_modulus / 3.0;
    (m / density).sqrt()
}

/// Shear (S-wave) speed: c_s = sqrt(G / ρ) \[m/s\].
pub fn s_wave_speed(shear_modulus: f64, density: f64) -> f64 {
    (shear_modulus / density).sqrt()
}

/// Rayleigh surface wave speed (Viktorov approximation):
/// c_R ≈ c_s · (0.862 + 1.14·ν) / (1 + ν).
pub fn rayleigh_wave_speed(shear_modulus: f64, density: f64, poisson_ratio: f64) -> f64 {
    let cs = s_wave_speed(shear_modulus, density);
    cs * (0.862 + 1.14 * poisson_ratio) / (1.0 + poisson_ratio)
}

/// Bar (longitudinal in thin rod) wave speed: c_bar = sqrt(E / ρ) \[m/s\].
pub fn bar_wave_speed(young_modulus: f64, density: f64) -> f64 {
    (young_modulus / density).sqrt()
}

/// Reflection coefficient at normal incidence between two media:
/// R = ((Z₂ - Z₁) / (Z₂ + Z₁))² (intensity).
pub fn acoustic_reflection_coefficient(z1: f64, z2: f64) -> f64 {
    if (z1 + z2).abs() < f64::EPSILON {
        return 0.0;
    }
    let r = (z2 - z1) / (z2 + z1);
    r * r
}

/// Transmission coefficient at normal incidence:
/// T = 1 - R = 4·Z₁·Z₂ / (Z₁ + Z₂)².
pub fn acoustic_transmission_coefficient(z1: f64, z2: f64) -> f64 {
    1.0 - acoustic_reflection_coefficient(z1, z2)
}

// ─── Thermal Expansion Stress ──────────────────────────────────────────────────

/// Thermal stress in a fully constrained bar:
/// σ = -E · α · ΔT  \[Pa\].
///
/// Negative = compression when ΔT > 0 (material wants to expand but is constrained).
pub fn thermal_stress_constrained(young_modulus: f64, alpha: f64, delta_t: f64) -> f64 {
    -young_modulus * alpha * delta_t
}

/// Free thermal strain: ε_th = α · ΔT (dimensionless).
pub fn thermal_strain(alpha: f64, delta_t: f64) -> f64 {
    alpha * delta_t
}

/// Mismatch stress at a bi-material interface (Suhir model, simplified):
/// σ_mismatch = E_eff · (α₂ - α₁) · ΔT / (1 - ν_eff)
///
/// where E_eff and ν_eff are harmonic averages.
pub fn bimaterial_mismatch_stress(
    e1: f64,
    nu1: f64,
    alpha1: f64,
    e2: f64,
    nu2: f64,
    alpha2: f64,
    delta_t: f64,
) -> f64 {
    let e_eff = 2.0 * e1 * e2 / (e1 + e2);
    let nu_eff = (nu1 + nu2) / 2.0;
    e_eff * (alpha2 - alpha1) * delta_t / (1.0 - nu_eff)
}

// ─── Material Search / Filter Utilities ───────────────────────────────────────

/// Filter a slice of `ExtendedMaterial` by minimum Young's modulus.
pub fn filter_by_min_stiffness(
    materials: &[ExtendedMaterial],
    min_e: f64,
) -> Vec<&ExtendedMaterial> {
    materials
        .iter()
        .filter(|m| m.young_modulus >= min_e)
        .collect()
}

/// Filter by maximum density (lightweight materials).
pub fn filter_by_max_density(
    materials: &[ExtendedMaterial],
    max_rho: f64,
) -> Vec<&ExtendedMaterial> {
    materials.iter().filter(|m| m.density <= max_rho).collect()
}

/// Find the material with the highest specific stiffness (E/ρ).
pub fn highest_specific_stiffness(materials: &[ExtendedMaterial]) -> Option<&ExtendedMaterial> {
    materials.iter().max_by(|a, b| {
        a.specific_stiffness()
            .partial_cmp(&b.specific_stiffness())
            .expect("operation should succeed")
    })
}

/// Find the material with the best thermal shock resistance.
pub fn best_thermal_shock_resistance(materials: &[ExtendedMaterial]) -> Option<&ExtendedMaterial> {
    materials.iter().max_by(|a, b| {
        a.thermal_shock_resistance()
            .partial_cmp(&b.thermal_shock_resistance())
            .expect("operation should succeed")
    })
}

/// Rank materials by an Ashby beam-stiffness index E^(1/2)/ρ.
///
/// Returns indices into the input slice sorted best → worst.
pub fn rank_by_beam_stiffness_index(materials: &[ExtendedMaterial]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..materials.len()).collect();
    indices.sort_by(|&a, &b| {
        let ia =
            MaterialIndex::beam_stiffness_index(materials[a].young_modulus, materials[a].density);
        let ib =
            MaterialIndex::beam_stiffness_index(materials[b].young_modulus, materials[b].density);
        ib.partial_cmp(&ia).unwrap_or(std::cmp::Ordering::Equal)
    });
    indices
}

// ─── Simple Material Catalogue ────────────────────────────────────────────────

/// Return a small catalogue of extended materials for testing and demonstration.
pub fn standard_material_catalogue() -> Vec<ExtendedMaterial> {
    vec![
        ExtendedMaterial::steel(),
        ExtendedMaterial::aluminium_6061(),
        ExtendedMaterial::titanium_6al4v(),
        ExtendedMaterial::cfrp_quasi_isotropic(),
        ExtendedMaterial::borosilicate_glass(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_material_library() {
        let mut lib = MaterialLibrary::new();
        assert!(lib.is_empty());
        lib.add(Material::new("steel", 7800.0, 0.6, 0.3));
        assert_eq!(lib.len(), 1);
        let steel = lib.get("steel").unwrap();
        assert_eq!(steel.density, 7800.0);
    }

    #[test]
    fn test_linear_elastic_bulk_shear_modulus_steel() {
        // Steel: E = 200 GPa, nu = 0.3
        let mat = LinearElastic::new(200.0e9, 0.3);
        let k = mat.bulk_modulus();
        let g = mat.shear_modulus();
        // K = 200e9 / (3*(1-0.6)) = 200e9 / 1.2 ≈ 166.67e9
        assert!((k - 166.667e9).abs() < 1.0e8);
        // G = 200e9 / (2*1.3) ≈ 76.923e9
        assert!((g - 76.923e9).abs() < 1.0e8);
    }

    #[test]
    fn test_linear_elastic_stress_strain_symmetry() {
        let mat = LinearElastic::new(200.0e9, 0.3);
        let c = mat.stress_strain_matrix_3d();
        for (i, row) in c.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - c[j][i]).abs() < 1.0e-6, "C[{i}][{j}] != C[{j}][{i}]");
            }
        }
    }

    #[test]
    fn test_neo_hookean_identity_zero_stress() {
        let mat = NeoHookean::new(1.0e6, 1.0e9);
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let p = mat.first_piola_kirchhoff_stress(&identity);
        for (i, row) in p.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(val.abs() < 1.0e-6, "P[{i}][{j}] = {val} should be ~0");
            }
        }
    }

    #[test]
    fn test_tait_eos_pressure_at_reference() {
        let eos = TaitEos::new(1000.0, 0.0, 1500.0, 7.0);
        let p = eos.pressure(1000.0);
        // At reference density, (rho/rho0)^gamma - 1 = 0, so pressure = reference_pressure = 0
        assert!(
            p.abs() < 1.0e-6,
            "Pressure at ref density should be ~0, got {p}"
        );
    }

    #[test]
    fn test_tait_eos_sound_speed() {
        let eos = TaitEos::new(1000.0, 0.0, 1500.0, 7.0);
        let cs = eos.sound_speed(1000.0);
        // At reference density, sound speed should equal c0
        assert!(
            (cs - 1500.0).abs() < 1.0e-6,
            "Sound speed at ref density should be 1500, got {cs}"
        );
    }

    #[test]
    fn test_friction_combine_average() {
        let r = combine_friction(0.4, 0.6, FrictionCombineRule::Average);
        assert!((r - 0.5).abs() < 1.0e-10);
    }

    #[test]
    fn test_friction_combine_min() {
        let r = combine_friction(0.4, 0.6, FrictionCombineRule::Min);
        assert!((r - 0.4).abs() < 1.0e-10);
    }

    #[test]
    fn test_friction_combine_max() {
        let r = combine_friction(0.4, 0.6, FrictionCombineRule::Max);
        assert!((r - 0.6).abs() < 1.0e-10);
    }

    #[test]
    fn test_friction_combine_multiply() {
        let r = combine_friction(0.4, 0.6, FrictionCombineRule::Multiply);
        assert!((r - 0.24).abs() < 1.0e-10);
    }

    #[test]
    fn test_restitution_combine_rules() {
        assert!(
            (combine_restitution(0.3, 0.7, RestitutionCombineRule::Average) - 0.5).abs() < 1.0e-10
        );
        assert!((combine_restitution(0.3, 0.7, RestitutionCombineRule::Min) - 0.3).abs() < 1.0e-10);
        assert!((combine_restitution(0.3, 0.7, RestitutionCombineRule::Max) - 0.7).abs() < 1.0e-10);
        assert!(
            (combine_restitution(0.3, 0.7, RestitutionCombineRule::Multiply) - 0.21).abs()
                < 1.0e-10
        );
    }

    #[test]
    fn test_material_remove() {
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("steel", 7800.0, 0.6, 0.3));
        lib.add(Material::new("rubber", 1100.0, 0.8, 0.9));
        assert_eq!(lib.len(), 2);
        let removed = lib.remove("steel");
        assert!(removed.is_some());
        assert_eq!(lib.len(), 1);
        assert!(!lib.contains("steel"));
        assert!(lib.contains("rubber"));
    }

    #[test]
    fn test_material_contains() {
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("wood", 600.0, 0.5, 0.4));
        assert!(lib.contains("wood"));
        assert!(!lib.contains("gold"));
    }

    #[test]
    fn test_material_names() {
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("a", 1.0, 0.1, 0.1));
        lib.add(Material::new("b", 2.0, 0.2, 0.2));
        let names = lib.names();
        assert_eq!(names.len(), 2);
        assert!(names.contains(&"a".to_string()));
        assert!(names.contains(&"b".to_string()));
    }

    #[test]
    fn test_material_iter() {
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("x", 1.0, 0.1, 0.1));
        lib.add(Material::new("y", 2.0, 0.2, 0.2));
        let count = lib.iter().count();
        assert_eq!(count, 2);
    }

    #[test]
    fn test_material_impedance() {
        let m = Material::new("steel", 7800.0, 0.6, 0.3);
        let z = m.impedance(5000.0);
        assert!((z - 39_000_000.0).abs() < 1.0);
    }

    #[test]
    fn test_material_is_heavier_than() {
        let steel = Material::new("steel", 7800.0, 0.6, 0.3);
        let wood = Material::new("wood", 600.0, 0.5, 0.4);
        assert!(steel.is_heavier_than(&wood));
        assert!(!wood.is_heavier_than(&steel));
    }

    #[test]
    fn test_material_is_bouncier_than() {
        let rubber = Material::new("rubber", 1100.0, 0.8, 0.9);
        let steel = Material::new("steel", 7800.0, 0.6, 0.3);
        assert!(rubber.is_bouncier_than(&steel));
    }

    #[test]
    fn test_material_combine_average() {
        let a = Material::new("a", 1000.0, 0.4, 0.2);
        let b = Material::new("b", 2000.0, 0.6, 0.8);
        let c = a.combine_average(&b);
        assert!((c.density - 1500.0).abs() < 1e-10);
        assert!((c.friction - 0.5).abs() < 1e-10);
        assert!((c.restitution - 0.5).abs() < 1e-10);
        assert_eq!(c.name, "a+b");
    }

    #[test]
    fn test_material_weight() {
        let m = Material::new("water", 1000.0, 0.0, 0.0);
        let w = m.weight(1.0, 9.81);
        assert!((w - 9810.0).abs() < 1e-6);
    }

    #[test]
    fn test_material_mass_from_volume() {
        let m = Material::new("iron", 7874.0, 0.5, 0.3);
        let mass = m.mass_from_volume(2.0);
        assert!((mass - 15748.0).abs() < 1e-6);
    }

    #[test]
    fn test_lerp_property() {
        assert!((lerp_property(0.0, 10.0, 0.5) - 5.0).abs() < 1e-10);
        assert!((lerp_property(0.0, 10.0, 0.0)).abs() < 1e-10);
        assert!((lerp_property(0.0, 10.0, 1.0) - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_dependent_restitution() {
        let e = velocity_dependent_restitution(0.8, 5.0, 0.01);
        // 0.8 * (1 - 0.01 * 5) = 0.8 * 0.95 = 0.76
        assert!((e - 0.76).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_dependent_restitution_clamped() {
        let e = velocity_dependent_restitution(0.8, 200.0, 0.01);
        // 0.8 * (1 - 0.01 * 200) = 0.8 * (-1.0) = -0.8 => clamped to 0
        assert!(e.abs() < 1e-10);
    }

    #[test]
    fn test_critical_damping() {
        let cc = critical_damping(100.0, 1.0);
        // 2 * sqrt(100) = 20
        assert!((cc - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_damping_ratio() {
        let zeta = damping_ratio(20.0, 100.0, 1.0);
        // c_crit = 20, so zeta = 20/20 = 1.0 (critically damped)
        assert!((zeta - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_natural_frequency() {
        let omega = natural_frequency(100.0, 1.0);
        assert!((omega - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_natural_frequency_zero_mass() {
        let omega = natural_frequency(100.0, 0.0);
        assert!(omega.abs() < 1e-10);
    }

    #[test]
    fn test_preset_materials_reasonable_values() {
        let materials = [
            presets::steel(),
            presets::aluminum(),
            presets::rubber(),
            presets::wood(),
            presets::glass(),
            presets::water(),
            presets::concrete(),
        ];
        for m in &materials {
            assert!(m.density > 0.0, "{} density should be positive", m.name);
            assert!(
                m.friction >= 0.0,
                "{} friction should be non-negative",
                m.name
            );
            assert!(
                m.restitution >= 0.0 && m.restitution <= 1.0,
                "{} restitution should be in [0, 1]",
                m.name
            );
        }
        // Steel should be denser than water
        assert!(presets::steel().density > presets::water().density);
        // Rubber should have highest restitution
        assert!(presets::rubber().restitution > presets::steel().restitution);
    }

    // ── ExtendedMaterial tests ─────────────────────────────────────────────────

    #[test]
    fn test_extended_material_steel_shear_modulus() {
        let s = ExtendedMaterial::steel();
        let g = s.shear_modulus();
        // For E=200e9, nu=0.3: G ≈ 76.9 GPa
        assert!((g - 76.9e9).abs() / 76.9e9 < 0.01, "G = {g}");
    }

    #[test]
    fn test_extended_material_steel_bulk_modulus() {
        let s = ExtendedMaterial::steel();
        let k = s.bulk_modulus();
        // For E=200e9, nu=0.3: K ≈ 166.7 GPa
        assert!((k - 166.7e9).abs() / 166.7e9 < 0.01, "K = {k}");
    }

    #[test]
    fn test_extended_material_thermal_diffusivity_steel() {
        let s = ExtendedMaterial::steel();
        let alpha = s.thermal_diffusivity();
        // λ=50, ρ=7850, cp=486 → α ≈ 1.31e-5 m²/s
        assert!(
            alpha > 1.0e-5 && alpha < 2.0e-5,
            "thermal diffusivity = {alpha}"
        );
    }

    #[test]
    fn test_extended_material_acoustic_impedance() {
        let s = ExtendedMaterial::steel();
        let z = s.acoustic_impedance();
        // Z = 7850 * 5960 ≈ 4.68e7
        assert!(z > 4.0e7, "steel acoustic impedance = {z}");
    }

    #[test]
    fn test_extended_material_specific_strength_ti_gt_steel() {
        let ti = ExtendedMaterial::titanium_6al4v();
        let fe = ExtendedMaterial::steel();
        assert!(
            ti.specific_strength() > fe.specific_strength(),
            "Ti specific strength ({}) should exceed steel ({})",
            ti.specific_strength(),
            fe.specific_strength()
        );
    }

    #[test]
    fn test_extended_material_thermal_shock_resistance_positive() {
        let al = ExtendedMaterial::aluminium_6061();
        let r = al.thermal_shock_resistance();
        assert!(r > 0.0, "thermal shock resistance = {r}");
    }

    #[test]
    fn test_extended_material_acoustic_reflection_coefficient() {
        let s = ExtendedMaterial::steel();
        let r = s.acoustic_reflection_coefficient_from_air();
        // Steel has very high impedance compared to air → R ≈ 1
        assert!(
            r > 0.99,
            "reflection from air-steel interface should be ~1, got {r}"
        );
    }

    #[test]
    fn test_extended_material_skin_depth_steel() {
        let s = ExtendedMaterial::steel();
        // At 1 MHz: δ = sqrt(1.7e-7 / (π * 1e6 * 4π*1e-7)) ≈ 6.6 μm for steel
        let delta = s.skin_depth(1.0e6);
        assert!(
            delta > 1.0e-6 && delta < 1.0e-3,
            "skin depth at 1 MHz = {delta}"
        );
    }

    #[test]
    fn test_extended_material_cfrp_lighter_than_steel() {
        let cfrp = ExtendedMaterial::cfrp_quasi_isotropic();
        let steel = ExtendedMaterial::steel();
        assert!(cfrp.density < steel.density);
    }

    #[test]
    fn test_extended_material_specific_stiffness_comparison() {
        let al = ExtendedMaterial::aluminium_6061();
        let cfrp = ExtendedMaterial::cfrp_quasi_isotropic();
        // Both should have reasonable specific stiffness
        assert!(al.specific_stiffness() > 0.0);
        assert!(cfrp.specific_stiffness() > 0.0);
    }

    // ── MaterialMixture tests ──────────────────────────────────────────────────

    #[test]
    fn test_voigt_modulus_at_zero_phi() {
        let e = MaterialMixture::voigt_modulus(100.0e9, 300.0e9, 0.0);
        assert!((e - 100.0e9).abs() < 1.0, "Voigt at φ=0 should be E1");
    }

    #[test]
    fn test_voigt_modulus_at_full_phi() {
        let e = MaterialMixture::voigt_modulus(100.0e9, 300.0e9, 1.0);
        assert!((e - 300.0e9).abs() < 1.0, "Voigt at φ=1 should be E2");
    }

    #[test]
    fn test_reuss_modulus_bounded_below_voigt() {
        let e_voigt = MaterialMixture::voigt_modulus(100.0e9, 300.0e9, 0.4);
        let e_reuss = MaterialMixture::reuss_modulus(100.0e9, 300.0e9, 0.4);
        assert!(
            e_reuss <= e_voigt + 1.0,
            "Reuss ({e_reuss}) should ≤ Voigt ({e_voigt})"
        );
    }

    #[test]
    fn test_hill_modulus_between_bounds() {
        let e1 = 70.0e9_f64;
        let e2 = 400.0e9_f64;
        let phi = 0.3_f64;
        let hill = MaterialMixture::hill_modulus(e1, e2, phi);
        let voigt = MaterialMixture::voigt_modulus(e1, e2, phi);
        let reuss = MaterialMixture::reuss_modulus(e1, e2, phi);
        assert!(
            hill >= reuss - 1.0 && hill <= voigt + 1.0,
            "Hill ({hill}) not between Reuss ({reuss}) and Voigt ({voigt})"
        );
    }

    #[test]
    fn test_mixture_density_linear() {
        let rho = MaterialMixture::mixture_density(1000.0, 3000.0, 0.5);
        assert!((rho - 2000.0).abs() < 1e-10);
    }

    #[test]
    fn test_turner_thermal_expansion_equal_phases() {
        // If both phases have the same α and K, the mixture α = that value
        let alpha_eff =
            MaterialMixture::turner_thermal_expansion(10e-6, 1.0e11, 10e-6, 1.0e11, 0.4);
        assert!((alpha_eff - 10e-6).abs() < 1e-14, "α_eff = {alpha_eff}");
    }

    #[test]
    fn test_hashin_shtrikman_upper_at_zero_phi() {
        // At φ=0 the mixture should approach phase 1
        let k_hs =
            MaterialMixture::hashin_shtrikman_upper_bulk(100.0e9, 300.0e9, 40.0e9, 100.0e9, 0.0);
        // Very close to k1 = 100e9 or k2 = 300e9 depending on ordering
        assert!((100.0e9 - 1.0..=300.0e9 + 1.0).contains(&k_hs));
    }

    // ── MaterialIndex tests ────────────────────────────────────────────────────

    #[test]
    fn test_beam_stiffness_index_positive() {
        let idx = MaterialIndex::beam_stiffness_index(200.0e9, 7850.0);
        assert!(idx > 0.0, "index = {idx}");
    }

    #[test]
    fn test_column_buckling_index_positive() {
        let idx = MaterialIndex::column_buckling_index(200.0e9, 7850.0);
        assert!(idx > 0.0, "index = {idx}");
    }

    #[test]
    fn test_beam_strength_index_positive() {
        let idx = MaterialIndex::beam_strength_index(355.0e6, 7850.0);
        assert!(idx > 0.0);
    }

    #[test]
    fn test_pressure_vessel_index_positive() {
        let idx = MaterialIndex::pressure_vessel_index(50.0e6, 355.0e6);
        assert!(idx > 0.0, "index = {idx}");
    }

    #[test]
    fn test_wear_resistance_index_positive() {
        let idx = MaterialIndex::wear_resistance_index(355.0e6, 200.0e9);
        assert!(idx > 0.0);
    }

    #[test]
    fn test_thermal_insulation_index_lower_for_higher_conductivity() {
        // Lower conductivity → lower thermal diffusivity → better insulator
        let i_steel = MaterialIndex::thermal_insulation_index(50.0, 7850.0, 486.0);
        let i_glass = MaterialIndex::thermal_insulation_index(1.0, 2500.0, 750.0);
        assert!(
            i_glass < i_steel,
            "glass insulation index should be lower: glass={i_glass}, steel={i_steel}"
        );
    }

    // ── Elastic moduli conversion tests ───────────────────────────────────────

    #[test]
    fn test_elastic_moduli_round_trip() {
        let e = 200.0e9_f64;
        let nu = 0.3_f64;
        let (k, g, _lambda) = elastic_moduli_from_e_nu(e, nu);
        let (e2, nu2, _) = elastic_moduli_from_kg(k, g);
        assert!((e2 - e).abs() / e < 1e-10, "E round-trip: {e2} vs {e}");
        assert!((nu2 - nu).abs() < 1e-10, "ν round-trip: {nu2} vs {nu}");
    }

    #[test]
    fn test_elastic_moduli_from_lame_round_trip() {
        let lambda = 115.38e9_f64; // steel Lamé
        let g = 76.92e9_f64;
        let (e, nu, k) = elastic_moduli_from_lame(lambda, g);
        assert!(e > 0.0);
        assert!(nu > 0.0 && nu < 0.5);
        assert!(k > 0.0);
    }

    #[test]
    fn test_elastic_admissible_steel() {
        assert!(elastic_moduli_admissible(200.0e9, 0.30));
    }

    #[test]
    fn test_elastic_admissible_negative_nu_auxetic() {
        // Auxetic material: ν in (-1, 0) is valid
        assert!(elastic_moduli_admissible(50.0e9, -0.5));
    }

    #[test]
    fn test_elastic_inadmissible_incompressible() {
        // ν = 0.5 → K → ∞: not admissible
        assert!(!elastic_moduli_admissible(200.0e9, 0.5));
    }

    // ── Wave speed tests ──────────────────────────────────────────────────────

    #[test]
    fn test_p_wave_speed_steel() {
        // Steel: K ≈ 166.7e9, G ≈ 76.9e9, ρ = 7850
        let (k, g, _) = elastic_moduli_from_e_nu(200.0e9, 0.3);
        let cp = p_wave_speed(k, g, 7850.0);
        // Expected ~5960 m/s
        assert!(cp > 5000.0 && cp < 7000.0, "P-wave speed = {cp}");
    }

    #[test]
    fn test_s_wave_speed_steel() {
        let g = 76.9e9_f64;
        let cs = s_wave_speed(g, 7850.0);
        // ≈ 3130 m/s
        assert!(cs > 2500.0 && cs < 4000.0, "S-wave speed = {cs}");
    }

    #[test]
    fn test_rayleigh_wave_speed_less_than_s() {
        let g = 76.9e9_f64;
        let cs = s_wave_speed(g, 7850.0);
        let cr = rayleigh_wave_speed(g, 7850.0, 0.3);
        assert!(
            cr < cs,
            "Rayleigh speed should be < S-wave speed: {cr} vs {cs}"
        );
    }

    #[test]
    fn test_bar_wave_speed_positive() {
        let cb = bar_wave_speed(200.0e9, 7850.0);
        assert!(cb > 0.0);
    }

    #[test]
    fn test_acoustic_reflection_same_medium() {
        // Z1 = Z2 → R = 0 (no reflection)
        let r = acoustic_reflection_coefficient(1.0e7, 1.0e7);
        assert!(r.abs() < 1e-15);
    }

    #[test]
    fn test_acoustic_transmission_plus_reflection_unity() {
        let z1 = 415.0_f64;
        let z2 = 46.8e6_f64;
        let r = acoustic_reflection_coefficient(z1, z2);
        let t = acoustic_transmission_coefficient(z1, z2);
        assert!(
            (r + t - 1.0).abs() < 1e-12,
            "R + T should = 1, got {}",
            r + t
        );
    }

    // ── Thermal stress tests ───────────────────────────────────────────────────

    #[test]
    fn test_thermal_stress_constrained_compressive() {
        // ΔT > 0 in constrained bar → compressive (negative) stress
        let sigma = thermal_stress_constrained(200.0e9, 12e-6, 100.0);
        assert!(
            sigma < 0.0,
            "Thermal stress should be compressive, got {sigma}"
        );
    }

    #[test]
    fn test_thermal_strain_positive() {
        let eps = thermal_strain(12e-6, 50.0);
        assert!((eps - 6e-4).abs() < 1e-15);
    }

    #[test]
    fn test_bimaterial_mismatch_stress_direction() {
        // α₂ > α₁ and ΔT > 0 → σ > 0 (phase 2 wants to expand more)
        let sigma = bimaterial_mismatch_stress(200.0e9, 0.30, 12e-6, 70.0e9, 0.33, 23e-6, 100.0);
        assert!(sigma > 0.0, "mismatch stress should be positive: {sigma}");
    }

    #[test]
    fn test_bimaterial_zero_mismatch_equal_alphas() {
        let sigma = bimaterial_mismatch_stress(200.0e9, 0.30, 12e-6, 70.0e9, 0.33, 12e-6, 200.0);
        assert!(sigma.abs() < 1e-3, "zero mismatch for equal α: {sigma}");
    }

    // ── Filter / search tests ──────────────────────────────────────────────────

    #[test]
    fn test_filter_by_min_stiffness() {
        let catalogue = standard_material_catalogue();
        let stiff = filter_by_min_stiffness(&catalogue, 100.0e9);
        // Steel, Al, Ti, glass all above 100 GPa; CFRP ≈ 70 GPa might be excluded
        assert!(!stiff.is_empty());
        for m in &stiff {
            assert!(
                m.young_modulus >= 100.0e9,
                "{} E={}",
                m.name,
                m.young_modulus
            );
        }
    }

    #[test]
    fn test_filter_by_max_density() {
        let catalogue = standard_material_catalogue();
        let light = filter_by_max_density(&catalogue, 3000.0);
        for m in &light {
            assert!(m.density <= 3000.0);
        }
    }

    #[test]
    fn test_highest_specific_stiffness_found() {
        let catalogue = standard_material_catalogue();
        let best = highest_specific_stiffness(&catalogue);
        assert!(best.is_some());
        let best_ss = best.unwrap().specific_stiffness();
        for m in &catalogue {
            assert!(m.specific_stiffness() <= best_ss + 1.0);
        }
    }

    #[test]
    fn test_best_thermal_shock_resistance_found() {
        let catalogue = standard_material_catalogue();
        let best = best_thermal_shock_resistance(&catalogue);
        assert!(best.is_some());
    }

    #[test]
    fn test_rank_by_beam_stiffness_index() {
        let catalogue = standard_material_catalogue();
        let ranked = rank_by_beam_stiffness_index(&catalogue);
        assert_eq!(ranked.len(), catalogue.len());
        // First-ranked should have highest E^0.5/ρ
        let best_idx = ranked[0];
        let best_val = MaterialIndex::beam_stiffness_index(
            catalogue[best_idx].young_modulus,
            catalogue[best_idx].density,
        );
        for &i in &ranked[1..] {
            let val = MaterialIndex::beam_stiffness_index(
                catalogue[i].young_modulus,
                catalogue[i].density,
            );
            assert!(val <= best_val + 1.0);
        }
    }

    #[test]
    fn test_standard_catalogue_count() {
        let catalogue = standard_material_catalogue();
        assert_eq!(catalogue.len(), 5);
    }
}

pub mod electromagnetic;

pub mod geomaterials;

pub mod nuclear_materials_advanced;

pub mod alloy_materials;
pub mod biomedical_materials;
pub mod composites_advanced;

pub mod biomechanical_materials;
pub mod ceramics_materials;
pub mod geomaterial_models;
pub mod polymers_materials;
pub mod semiconductor_materials;
pub mod shape_memory_alloy;
pub mod tribology_materials;
