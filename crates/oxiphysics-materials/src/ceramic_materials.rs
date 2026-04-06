// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ceramic material models: property database, fracture mechanics, thermal shock,
//! creep, sintering densification, hardness conversion, Weibull failure statistics,
//! piezoelectric ceramics (PZT), ferroelectric hysteresis, zirconia phase
//! transformation toughening, alumina properties, thermal conductivity models,
//! Hasselman parameters, and grain growth kinetics.
//!
//! All functions use SI units unless otherwise stated.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

/// Universal gas constant (J/(mol K)).
const R_GAS: f64 = 8.314;

// ---------------------------------------------------------------------------
// CeramicMaterial
// ---------------------------------------------------------------------------

/// Structural and thermal properties of a ceramic material.
#[derive(Debug, Clone)]
pub struct CeramicMaterial {
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poissons_ratio: f64,
    /// Mode-I fracture toughness KIC (Pa sqrt(m)).
    pub fracture_toughness_kic: f64,
    /// Vickers hardness (HV).
    pub hardness_vickers: f64,
    /// Thermal conductivity (W/(m K)).
    pub thermal_conductivity: f64,
    /// Coefficient of thermal expansion (1/K).
    pub thermal_expansion: f64,
    /// Melting / decomposition point (K).
    pub melting_point: f64,
    /// Average grain size (um).
    pub grain_size_um: f64,
}

// ---------------------------------------------------------------------------
// CeramicType
// ---------------------------------------------------------------------------

/// Enumeration of common engineering ceramics supported by the preset database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CeramicType {
    /// Aluminium oxide (Al2O3) -- general-purpose structural ceramic.
    Alumina,
    /// Yttria-stabilised zirconia (3Y-TZP) -- high-toughness, biomedical use.
    ZirconiaYSZ,
    /// Silicon carbide (SiC) -- extreme hardness and thermal stability.
    SiliconCarbide,
    /// Silicon nitride (Si3N4) -- excellent fatigue and wear resistance.
    SiliconNitride,
    /// Boron carbide (B4C) -- one of the hardest known materials; armour ceramic.
    BoronCarbide,
    /// Hydroxyapatite (Ca10(PO4)6(OH)2) -- bone-like biocompatible ceramic.
    HydroxyApatite,
}

// ---------------------------------------------------------------------------
// Preset database
// ---------------------------------------------------------------------------

/// Return a [`CeramicMaterial`] with representative published properties for
/// the requested [`CeramicType`].
///
/// Sources: Engineered Materials Handbook Vol. 4 (ASM); CRC Handbook; Munz & Fett.
pub fn ceramic_preset(ct: CeramicType) -> CeramicMaterial {
    match ct {
        CeramicType::Alumina => CeramicMaterial {
            youngs_modulus: 380.0e9,
            poissons_ratio: 0.22,
            fracture_toughness_kic: 3.5e6,
            hardness_vickers: 1600.0,
            thermal_conductivity: 25.0,
            thermal_expansion: 8.1e-6,
            melting_point: 2323.0,
            grain_size_um: 3.0,
        },
        CeramicType::ZirconiaYSZ => CeramicMaterial {
            youngs_modulus: 210.0e9,
            poissons_ratio: 0.31,
            fracture_toughness_kic: 8.0e6,
            hardness_vickers: 1200.0,
            thermal_conductivity: 2.2,
            thermal_expansion: 10.5e-6,
            melting_point: 2988.0,
            grain_size_um: 0.5,
        },
        CeramicType::SiliconCarbide => CeramicMaterial {
            youngs_modulus: 410.0e9,
            poissons_ratio: 0.14,
            fracture_toughness_kic: 3.0e6,
            hardness_vickers: 2500.0,
            thermal_conductivity: 120.0,
            thermal_expansion: 4.5e-6,
            melting_point: 3003.0,
            grain_size_um: 5.0,
        },
        CeramicType::SiliconNitride => CeramicMaterial {
            youngs_modulus: 310.0e9,
            poissons_ratio: 0.27,
            fracture_toughness_kic: 7.0e6,
            hardness_vickers: 1600.0,
            thermal_conductivity: 30.0,
            thermal_expansion: 3.2e-6,
            melting_point: 2173.0,
            grain_size_um: 1.0,
        },
        CeramicType::BoronCarbide => CeramicMaterial {
            youngs_modulus: 460.0e9,
            poissons_ratio: 0.17,
            fracture_toughness_kic: 2.5e6,
            hardness_vickers: 3000.0,
            thermal_conductivity: 35.0,
            thermal_expansion: 5.0e-6,
            melting_point: 2618.0,
            grain_size_um: 10.0,
        },
        CeramicType::HydroxyApatite => CeramicMaterial {
            youngs_modulus: 80.0e9,
            poissons_ratio: 0.27,
            fracture_toughness_kic: 0.9e6,
            hardness_vickers: 500.0,
            thermal_conductivity: 1.3,
            thermal_expansion: 14.0e-6,
            melting_point: 1943.0,
            grain_size_um: 2.0,
        },
    }
}

// ---------------------------------------------------------------------------
// Hall-Petch strength
// ---------------------------------------------------------------------------

/// Hall-Petch grain-boundary strengthening: sigma = sigma_0 + k_hp / sqrt(d).
///
/// # Arguments
/// * `k_hp` -- Hall-Petch coefficient (Pa m^0.5).
/// * `sigma_0` -- friction stress / lattice resistance (Pa).
/// * `grain_size` -- grain size d (m).
///
/// Returns yield/fracture stress (Pa).
pub fn hall_petch_strength(k_hp: f64, sigma_0: f64, grain_size: f64) -> f64 {
    sigma_0 + k_hp / grain_size.sqrt()
}

// ---------------------------------------------------------------------------
// Griffith critical flaw size
// ---------------------------------------------------------------------------

/// Griffith / Irwin criterion -- critical flaw radius *a* that drives unstable fracture.
///
/// `a_c = (1/pi) * (K_IC / (Y * sigma))^2`
///
/// # Arguments
/// * `kic` -- fracture toughness (Pa sqrt(m)).
/// * `stress` -- applied tensile stress (Pa).
/// * `geometry_factor` -- dimensionless stress-intensity shape factor Y.
///
/// Returns the critical half-crack length (m).
pub fn critical_flaw_size(kic: f64, stress: f64, geometry_factor: f64) -> f64 {
    let ratio = kic / (geometry_factor * stress);
    ratio * ratio / PI
}

// ---------------------------------------------------------------------------
// Thermal shock resistance
// ---------------------------------------------------------------------------

/// Kingery first thermal-shock resistance parameter R.
///
/// `R = sigma_f * (1 - nu) / (E * alpha)`
///
/// A higher R indicates better thermal-shock resistance.
///
/// # Arguments
/// * `ceramic` -- material properties.
/// * `delta_t` -- temperature difference (K).
///
/// Returns the R parameter (K).
pub fn thermal_shock_resistance(ceramic: &CeramicMaterial, delta_t: f64) -> f64 {
    let grain_m = ceramic.grain_size_um * 1e-6 / 2.0;
    let sigma_f = ceramic.fracture_toughness_kic / (PI * grain_m).sqrt();
    let r = sigma_f * (1.0 - ceramic.poissons_ratio)
        / (ceramic.youngs_modulus * ceramic.thermal_expansion);
    r * delta_t
}

// ---------------------------------------------------------------------------
// Hasselman thermal shock parameters
// ---------------------------------------------------------------------------

/// Hasselman thermal shock parameters for ceramic design.
///
/// Contains the four Hasselman figures of merit: R, R', R''', R''''.
#[derive(Debug, Clone)]
pub struct HasselmanParameters {
    /// First parameter R (K) -- resistance to fracture initiation.
    pub r_first: f64,
    /// Second parameter R' (W/m) -- R scaled by thermal conductivity.
    pub r_second: f64,
    /// Third parameter R''' (m^-1) -- resistance to crack propagation.
    pub r_third: f64,
    /// Fourth parameter R'''' (m Pa) -- damage resistance.
    pub r_fourth: f64,
}

/// Compute all four Hasselman thermal shock resistance parameters.
///
/// # Arguments
/// * `strength` -- tensile/flexural strength sigma_f (Pa).
/// * `youngs_modulus` -- E (Pa).
/// * `poissons_ratio` -- nu (dimensionless).
/// * `cte` -- coefficient of thermal expansion alpha (1/K).
/// * `thermal_conductivity` -- k (W/(m K)).
/// * `fracture_toughness` -- K_IC (Pa sqrt(m)).
/// * `surface_energy` -- gamma_f (J/m^2).
///
/// Returns [`HasselmanParameters`].
pub fn hasselman_parameters(
    strength: f64,
    youngs_modulus: f64,
    poissons_ratio: f64,
    cte: f64,
    thermal_conductivity: f64,
    fracture_toughness: f64,
    surface_energy: f64,
) -> HasselmanParameters {
    // R = sigma_f (1 - nu) / (E alpha)
    let r_first = strength * (1.0 - poissons_ratio) / (youngs_modulus * cte);
    // R' = k * R
    let r_second = thermal_conductivity * r_first;
    // R''' = E / (K_IC^2 (1 - nu^2))  -- resistance to crack propagation
    let r_third = youngs_modulus
        / (fracture_toughness * fracture_toughness * (1.0 - poissons_ratio * poissons_ratio));
    // R'''' = E * gamma_f / (sigma_f^2 (1 - nu^2))  -- damage resistance
    let r_fourth = youngs_modulus * surface_energy
        / (strength * strength * (1.0 - poissons_ratio * poissons_ratio));
    HasselmanParameters {
        r_first,
        r_second,
        r_third,
        r_fourth,
    }
}

// ---------------------------------------------------------------------------
// Nabarro-Herring creep
// ---------------------------------------------------------------------------

/// Nabarro-Herring diffusional creep rate for fine-grained ceramics.
///
/// # Arguments
/// * `temp` -- absolute temperature (K).
/// * `stress` -- applied stress (Pa).
/// * `grain_size` -- grain diameter (m).
/// * `activation_energy` -- activation energy for grain-boundary diffusion (J/mol).
///
/// Returns creep strain rate (1/s, scaled).
pub fn creep_rate_ceramic(temp: f64, stress: f64, grain_size: f64, activation_energy: f64) -> f64 {
    let exponent = -activation_energy / (R_GAS * temp);
    let a = 1e-20_f64;
    a * stress * exponent.exp() / (grain_size * grain_size)
}

// ---------------------------------------------------------------------------
// Sintering densification
// ---------------------------------------------------------------------------

/// Simplified sintering densification kinetics (Coble-type).
///
/// # Arguments
/// * `temp` -- sintering temperature (K).
/// * `time` -- sintering time (s).
/// * `initial_density` -- green density as fraction of theoretical (0-1).
///
/// Returns instantaneous relative density (fraction, clamped to 1.0).
pub fn sintering_densification(temp: f64, time: f64, initial_density: f64) -> f64 {
    const Q: f64 = 300_000.0;
    let tau = (-Q / (R_GAS * temp)).exp();
    let rho_inf = 1.0_f64;
    let rho = rho_inf - (rho_inf - initial_density) * (-time * tau).exp();
    rho.min(1.0)
}

// ---------------------------------------------------------------------------
// Grain growth kinetics
// ---------------------------------------------------------------------------

/// Normal grain growth model: d^n - d0^n = K * t * exp(-Q/(R*T)).
///
/// # Arguments
/// * `d0` -- initial grain size (m).
/// * `exponent_n` -- grain growth exponent (typically 2-4).
/// * `prefactor_k` -- pre-exponential growth constant (m^n / s).
/// * `activation_energy` -- activation energy Q (J/mol).
/// * `temp` -- absolute temperature (K).
/// * `time` -- annealing time (s).
///
/// Returns final grain size (m).
pub fn grain_growth_size(
    d0: f64,
    exponent_n: f64,
    prefactor_k: f64,
    activation_energy: f64,
    temp: f64,
    time: f64,
) -> f64 {
    let d0n = d0.powf(exponent_n);
    let kt = prefactor_k * time * (-activation_energy / (R_GAS * temp)).exp();
    (d0n + kt).powf(1.0 / exponent_n)
}

/// Grain growth rate dd/dt at current grain size d and temperature.
///
/// # Arguments
/// * `d` -- current grain size (m).
/// * `exponent_n` -- grain growth exponent.
/// * `prefactor_k` -- pre-exponential growth constant (m^n / s).
/// * `activation_energy` -- activation energy Q (J/mol).
/// * `temp` -- absolute temperature (K).
///
/// Returns growth rate (m/s).
pub fn grain_growth_rate(
    d: f64,
    exponent_n: f64,
    prefactor_k: f64,
    activation_energy: f64,
    temp: f64,
) -> f64 {
    let k_eff = prefactor_k * (-activation_energy / (R_GAS * temp)).exp();
    k_eff / (exponent_n * d.powf(exponent_n - 1.0))
}

// ---------------------------------------------------------------------------
// Vickers hardness conversion
// ---------------------------------------------------------------------------

/// Convert Vickers hardness number to equivalent pressure in MPa.
///
/// `HV [MPa] = HV * 9.807`
pub fn vickers_hardness_to_mpa(hv: f64) -> f64 {
    hv * 9.807
}

// ---------------------------------------------------------------------------
// Weibull fracture probability
// ---------------------------------------------------------------------------

/// Two-parameter Weibull failure probability under uniaxial tension.
///
/// `P_f = 1 - exp(-(sigma / sigma_0)^m)`
///
/// # Arguments
/// * `stress` -- applied stress (Pa).
/// * `modulus_of_rupture` -- Weibull scale parameter sigma_0 (Pa).
/// * `weibull_m` -- Weibull shape / modulus m (dimensionless).
///
/// Returns failure probability (0-1).
pub fn ceramic_fracture_probability(stress: f64, modulus_of_rupture: f64, weibull_m: f64) -> f64 {
    let ratio = stress / modulus_of_rupture;
    1.0 - (-(ratio.powf(weibull_m))).exp()
}

/// Weibull strength at a specified failure probability.
///
/// `sigma = sigma_0 * (-ln(1 - P_f))^(1/m)`
///
/// # Arguments
/// * `prob_failure` -- target failure probability (0 < P_f < 1).
/// * `modulus_of_rupture` -- scale parameter sigma_0 (Pa).
/// * `weibull_m` -- shape parameter m.
///
/// Returns stress (Pa).
pub fn weibull_strength_at_probability(
    prob_failure: f64,
    modulus_of_rupture: f64,
    weibull_m: f64,
) -> f64 {
    let ln_term = -(1.0 - prob_failure).ln();
    modulus_of_rupture * ln_term.powf(1.0 / weibull_m)
}

/// Weibull effective volume scaling: sigma_2 = sigma_1 * (V1/V2)^(1/m).
///
/// Scales strength from volume V1 to volume V2 via the weakest-link model.
///
/// # Arguments
/// * `sigma_1` -- strength at reference volume V1 (Pa).
/// * `v1` -- reference volume (m^3).
/// * `v2` -- target volume (m^3).
/// * `weibull_m` -- Weibull modulus m.
///
/// Returns scaled strength (Pa).
pub fn weibull_volume_scaling(sigma_1: f64, v1: f64, v2: f64, weibull_m: f64) -> f64 {
    sigma_1 * (v1 / v2).powf(1.0 / weibull_m)
}

// ---------------------------------------------------------------------------
// Thermal conductivity models
// ---------------------------------------------------------------------------

/// Debye thermal conductivity model for crystalline ceramics.
///
/// `k = (1/3) * C_v * v_s * l_mfp`
///
/// # Arguments
/// * `heat_capacity_vol` -- volumetric heat capacity C_v (J/(m^3 K)).
/// * `sound_velocity` -- average phonon velocity v_s (m/s).
/// * `mean_free_path` -- phonon mean free path l (m).
///
/// Returns thermal conductivity (W/(m K)).
pub fn debye_thermal_conductivity(
    heat_capacity_vol: f64,
    sound_velocity: f64,
    mean_free_path: f64,
) -> f64 {
    heat_capacity_vol * sound_velocity * mean_free_path / 3.0
}

/// Klemens point-defect phonon scattering model.
///
/// Gives the reciprocal mean free path due to point defects:
/// `1/l_pd = A * omega^4 / v_s`
///
/// where A is a scattering parameter and omega is phonon frequency.
///
/// # Arguments
/// * `scatter_param` -- point-defect scattering parameter A (s^3/m).
/// * `omega` -- phonon angular frequency (rad/s).
/// * `sound_velocity` -- v_s (m/s).
///
/// Returns reciprocal mean free path (1/m).
pub fn klemens_point_defect_scattering(scatter_param: f64, omega: f64, sound_velocity: f64) -> f64 {
    scatter_param * omega.powi(4) / sound_velocity
}

/// Effective medium thermal conductivity for porous ceramics (Maxwell-Eucken).
///
/// `k_eff = k_solid * (1 - porosity) / (1 + porosity / 2)`
///
/// # Arguments
/// * `k_solid` -- thermal conductivity of dense ceramic (W/(m K)).
/// * `porosity` -- volume fraction of pores (0-1).
///
/// Returns effective thermal conductivity (W/(m K)).
pub fn porous_thermal_conductivity(k_solid: f64, porosity: f64) -> f64 {
    k_solid * (1.0 - porosity) / (1.0 + porosity / 2.0)
}

/// Temperature-dependent thermal conductivity via 1/T Umklapp model.
///
/// `k(T) = k_ref * (T_ref / T)`
///
/// Valid above the Debye temperature where Umklapp scattering dominates.
///
/// # Arguments
/// * `k_ref` -- conductivity at reference temperature (W/(m K)).
/// * `t_ref` -- reference temperature (K).
/// * `t` -- target temperature (K).
///
/// Returns thermal conductivity at temperature T (W/(m K)).
pub fn umklapp_conductivity(k_ref: f64, t_ref: f64, t: f64) -> f64 {
    k_ref * t_ref / t
}

// ---------------------------------------------------------------------------
// Piezoelectric ceramics (PZT)
// ---------------------------------------------------------------------------

/// Properties of a piezoelectric ceramic (e.g. PZT).
#[derive(Debug, Clone)]
pub struct PiezoelectricCeramic {
    /// Piezoelectric charge coefficient d33 (C/N or m/V).
    pub d33: f64,
    /// Piezoelectric charge coefficient d31 (C/N).
    pub d31: f64,
    /// Piezoelectric voltage coefficient g33 (V m/N).
    pub g33: f64,
    /// Relative permittivity (dielectric constant) epsilon_r.
    pub epsilon_r: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Electromechanical coupling factor k33 (dimensionless).
    pub k33: f64,
    /// Curie temperature (K).
    pub curie_temperature: f64,
    /// Mechanical quality factor Q_m (dimensionless).
    pub quality_factor: f64,
}

/// Vacuum permittivity (F/m).
const EPS_0: f64 = 8.854_187_817e-12;

/// PZT-4 (Navy Type I) hard PZT preset.
pub fn pzt4_preset() -> PiezoelectricCeramic {
    PiezoelectricCeramic {
        d33: 289.0e-12,
        d31: -123.0e-12,
        g33: 26.1e-3,
        epsilon_r: 1300.0,
        youngs_modulus: 66.0e9,
        k33: 0.70,
        curie_temperature: 601.0,
        quality_factor: 500.0,
    }
}

/// PZT-5A (Navy Type II) soft PZT preset.
pub fn pzt5a_preset() -> PiezoelectricCeramic {
    PiezoelectricCeramic {
        d33: 374.0e-12,
        d31: -171.0e-12,
        g33: 24.8e-3,
        epsilon_r: 1700.0,
        youngs_modulus: 52.0e9,
        k33: 0.72,
        curie_temperature: 638.0,
        quality_factor: 75.0,
    }
}

/// PZT-5H (Navy Type VI) soft PZT preset.
pub fn pzt5h_preset() -> PiezoelectricCeramic {
    PiezoelectricCeramic {
        d33: 593.0e-12,
        d31: -274.0e-12,
        g33: 19.7e-3,
        epsilon_r: 3400.0,
        youngs_modulus: 50.0e9,
        k33: 0.75,
        curie_temperature: 466.0,
        quality_factor: 65.0,
    }
}

/// Compute piezoelectric strain from applied electric field.
///
/// `S_33 = d33 * E_3`
///
/// # Arguments
/// * `pzt` -- piezoelectric ceramic properties.
/// * `e_field` -- applied electric field E_3 (V/m).
///
/// Returns strain (dimensionless).
pub fn piezo_strain_33(pzt: &PiezoelectricCeramic, e_field: f64) -> f64 {
    pzt.d33 * e_field
}

/// Compute piezoelectric transverse strain from applied electric field.
///
/// `S_11 = d31 * E_3`
///
/// # Arguments
/// * `pzt` -- piezoelectric ceramic properties.
/// * `e_field` -- applied electric field E_3 (V/m).
///
/// Returns transverse strain (dimensionless, typically negative for d31 < 0).
pub fn piezo_strain_31(pzt: &PiezoelectricCeramic, e_field: f64) -> f64 {
    pzt.d31 * e_field
}

/// Compute open-circuit voltage from applied stress on a piezoelectric element.
///
/// `V = g33 * sigma * t`
///
/// # Arguments
/// * `pzt` -- piezoelectric ceramic properties.
/// * `stress` -- applied mechanical stress (Pa).
/// * `thickness` -- element thickness (m).
///
/// Returns voltage (V).
pub fn piezo_open_circuit_voltage(pzt: &PiezoelectricCeramic, stress: f64, thickness: f64) -> f64 {
    pzt.g33 * stress * thickness
}

/// Compute charge generated by a piezoelectric element under stress.
///
/// `Q = d33 * F = d33 * sigma * A`
///
/// # Arguments
/// * `pzt` -- piezoelectric ceramic properties.
/// * `stress` -- applied stress (Pa).
/// * `area` -- electrode area (m^2).
///
/// Returns charge (C).
pub fn piezo_charge_generated(pzt: &PiezoelectricCeramic, stress: f64, area: f64) -> f64 {
    pzt.d33 * stress * area
}

/// Compute stored electrical energy density in a piezoelectric element.
///
/// `u_e = 0.5 * epsilon_0 * epsilon_r * E^2`
///
/// # Arguments
/// * `pzt` -- piezoelectric ceramic properties.
/// * `e_field` -- electric field (V/m).
///
/// Returns energy density (J/m^3).
pub fn piezo_energy_density(pzt: &PiezoelectricCeramic, e_field: f64) -> f64 {
    0.5 * EPS_0 * pzt.epsilon_r * e_field * e_field
}

/// Compute resonant frequency of a piezoelectric disc in thickness mode.
///
/// `f_r = v_s / (2 * t)`
///
/// where `v_s = sqrt(E / rho)` and we approximate density from the quality factor
/// relationship. For a simpler model we accept sound velocity directly.
///
/// # Arguments
/// * `sound_velocity` -- longitudinal sound velocity in the ceramic (m/s).
/// * `thickness` -- element thickness (m).
///
/// Returns resonant frequency (Hz).
pub fn piezo_resonant_frequency(sound_velocity: f64, thickness: f64) -> f64 {
    sound_velocity / (2.0 * thickness)
}

// ---------------------------------------------------------------------------
// Ferroelectric hysteresis
// ---------------------------------------------------------------------------

/// Ferroelectric hysteresis state for P-E loop modelling.
///
/// Uses a simplified hyperbolic tangent saturation model:
/// `P = P_s * tanh((E - E_c * sign) / E_a)`
///
/// where P_s is saturation polarisation, E_c is the coercive field,
/// and E_a is a shape parameter controlling the loop width.
#[derive(Debug, Clone)]
pub struct FerroelectricHysteresis {
    /// Saturation polarisation P_s (C/m^2).
    pub saturation_polarisation: f64,
    /// Coercive field E_c (V/m).
    pub coercive_field: f64,
    /// Shape parameter E_a (V/m), controls steepness of the transition.
    pub shape_param: f64,
    /// Remanent polarisation P_r (C/m^2).
    pub remanent_polarisation: f64,
}

impl FerroelectricHysteresis {
    /// Create a new ferroelectric hysteresis model.
    pub fn new(
        saturation_polarisation: f64,
        coercive_field: f64,
        shape_param: f64,
        remanent_polarisation: f64,
    ) -> Self {
        Self {
            saturation_polarisation,
            coercive_field,
            shape_param,
            remanent_polarisation,
        }
    }

    /// BaTiO3 typical values.
    pub fn barium_titanate() -> Self {
        Self::new(0.26, 1.0e5, 5.0e4, 0.22)
    }

    /// PZT (soft) typical values.
    pub fn pzt_soft() -> Self {
        Self::new(0.35, 1.4e5, 6.0e4, 0.30)
    }

    /// Compute polarisation on the upper branch (increasing E).
    ///
    /// `P = P_s * tanh((E + E_c) / E_a)`
    ///
    /// # Arguments
    /// * `e_field` -- applied electric field (V/m).
    ///
    /// Returns polarisation (C/m^2).
    pub fn polarisation_upper(&self, e_field: f64) -> f64 {
        self.saturation_polarisation * ((e_field + self.coercive_field) / self.shape_param).tanh()
    }

    /// Compute polarisation on the lower branch (decreasing E).
    ///
    /// `P = P_s * tanh((E - E_c) / E_a)`
    ///
    /// # Arguments
    /// * `e_field` -- applied electric field (V/m).
    ///
    /// Returns polarisation (C/m^2).
    pub fn polarisation_lower(&self, e_field: f64) -> f64 {
        self.saturation_polarisation * ((e_field - self.coercive_field) / self.shape_param).tanh()
    }

    /// Estimate hysteresis loss per cycle from the coercive field and remanent
    /// polarisation (simplified rectangular approximation).
    ///
    /// `W = 4 * E_c * P_r`
    ///
    /// Returns energy density per cycle (J/m^3).
    pub fn hysteresis_loss_per_cycle(&self) -> f64 {
        4.0 * self.coercive_field * self.remanent_polarisation
    }

    /// Curie-Weiss permittivity above the Curie temperature.
    ///
    /// `epsilon_r = C / (T - T_c)`
    ///
    /// # Arguments
    /// * `curie_constant` -- Curie constant C (K).
    /// * `curie_temp` -- Curie temperature T_c (K).
    /// * `temp` -- current temperature (K).
    ///
    /// Returns relative permittivity (dimensionless), or `f64::INFINITY` if T <= T_c.
    pub fn curie_weiss_permittivity(curie_constant: f64, curie_temp: f64, temp: f64) -> f64 {
        if temp <= curie_temp {
            f64::INFINITY
        } else {
            curie_constant / (temp - curie_temp)
        }
    }
}

// ---------------------------------------------------------------------------
// Zirconia phase transformation toughening
// ---------------------------------------------------------------------------

/// Zirconia tetragonal-to-monoclinic transformation toughening model.
///
/// Based on the Budiansky-Hutchinson-Lambropoulos framework.
#[derive(Debug, Clone)]
pub struct ZirconiaToughening {
    /// Dilatational transformation strain epsilon_T (dimensionless, ~0.04).
    pub transformation_strain: f64,
    /// Volume fraction of transformable tetragonal phase f_v (0-1).
    pub volume_fraction: f64,
    /// Shear modulus G (Pa).
    pub shear_modulus: f64,
    /// Base intrinsic fracture toughness K_IC0 (Pa sqrt(m)).
    pub base_toughness: f64,
    /// Transformation zone half-height h (m).
    pub zone_height: f64,
}

impl ZirconiaToughening {
    /// Create a new zirconia toughening model.
    pub fn new(
        transformation_strain: f64,
        volume_fraction: f64,
        shear_modulus: f64,
        base_toughness: f64,
        zone_height: f64,
    ) -> Self {
        Self {
            transformation_strain,
            volume_fraction,
            shear_modulus,
            base_toughness,
            zone_height,
        }
    }

    /// 3Y-TZP (3 mol% yttria-stabilised tetragonal zirconia polycrystal) preset.
    pub fn y_tzp_3mol() -> Self {
        Self::new(0.04, 0.30, 80.0e9, 5.0e6, 50.0e-6)
    }

    /// Toughening increment delta_K_IC (Pa sqrt(m)).
    ///
    /// `delta_K = 0.22 * E_eff * epsilon_T * f_v * sqrt(h)`
    ///
    /// where `E_eff = 2 * G`.
    pub fn toughening_increment(&self) -> f64 {
        let e_eff = 2.0 * self.shear_modulus;
        0.22 * e_eff * self.transformation_strain * self.volume_fraction * self.zone_height.sqrt()
    }

    /// Effective fracture toughness K_IC_eff (Pa sqrt(m)).
    pub fn effective_toughness(&self) -> f64 {
        self.base_toughness + self.toughening_increment()
    }

    /// Transformation stress sigma_t = G * epsilon_T / (1 - f_v) (Pa).
    pub fn transformation_stress(&self) -> f64 {
        self.shear_modulus * self.transformation_strain / (1.0 - self.volume_fraction)
    }

    /// Estimated transformation zone size from Irwin-type analysis (m).
    ///
    /// `h_est = (K_IC_eff / sigma_t)^2 / (3 pi)`
    pub fn zone_size_estimate(&self) -> f64 {
        let sigma_t = self.transformation_stress();
        let kic = self.effective_toughness();
        (kic / sigma_t).powi(2) / (3.0 * PI)
    }

    /// Critical grain size above which spontaneous t->m transformation occurs.
    ///
    /// Approximation: `d_c = 6 * gamma / (delta_G_chem + E * epsilon_T^2 / (1-nu))`
    ///
    /// # Arguments
    /// * `surface_energy` -- interface energy gamma (J/m^2).
    /// * `chem_driving_force` -- chemical free energy difference delta_G (J/m^3).
    /// * `youngs_modulus` -- E (Pa).
    /// * `poissons_ratio` -- nu.
    ///
    /// Returns critical grain size (m).
    pub fn critical_grain_size(
        &self,
        surface_energy: f64,
        chem_driving_force: f64,
        youngs_modulus: f64,
        poissons_ratio: f64,
    ) -> f64 {
        let strain_energy =
            youngs_modulus * self.transformation_strain.powi(2) / (1.0 - poissons_ratio);
        6.0 * surface_energy / (chem_driving_force + strain_energy)
    }
}

// ---------------------------------------------------------------------------
// Alumina property models
// ---------------------------------------------------------------------------

/// Temperature-dependent Young's modulus of alumina.
///
/// Empirical fit: `E(T) = E_0 * (1 - b * (T - T_ref))`
///
/// where b ~ 5e-5 1/K for alumina.
///
/// # Arguments
/// * `e0` -- room temperature modulus (Pa).
/// * `temp` -- current temperature (K).
/// * `t_ref` -- reference temperature (K), typically 300.
/// * `b_coeff` -- temperature coefficient (1/K), ~5e-5 for alumina.
///
/// Returns Young's modulus at temperature T (Pa).
pub fn alumina_modulus_vs_temp(e0: f64, temp: f64, t_ref: f64, b_coeff: f64) -> f64 {
    (e0 * (1.0 - b_coeff * (temp - t_ref))).max(0.0)
}

/// Alumina flexural strength size effect (Weibull volume scaling).
///
/// `sigma(V) = sigma_ref * (V_ref / V)^(1/m)`
///
/// # Arguments
/// * `sigma_ref` -- reference strength (Pa).
/// * `v_ref` -- reference test volume (m^3).
/// * `v_target` -- target component volume (m^3).
/// * `weibull_m` -- Weibull modulus.
///
/// Returns scaled strength (Pa).
pub fn alumina_size_effect(sigma_ref: f64, v_ref: f64, v_target: f64, weibull_m: f64) -> f64 {
    weibull_volume_scaling(sigma_ref, v_ref, v_target, weibull_m)
}

/// Alumina thermal conductivity vs temperature (empirical inverse-T fit).
///
/// `k(T) = A / T + B`
///
/// # Arguments
/// * `a_coeff` -- fitting constant A (W/m).
/// * `b_coeff` -- fitting constant B (W/(m K)).
/// * `temp` -- temperature (K).
///
/// Returns thermal conductivity (W/(m K)).
pub fn alumina_thermal_conductivity_vs_temp(a_coeff: f64, b_coeff: f64, temp: f64) -> f64 {
    a_coeff / temp + b_coeff
}

// ---------------------------------------------------------------------------
// Indentation fracture toughness (Anstis method)
// ---------------------------------------------------------------------------

/// Indentation fracture toughness from Vickers indentation (Anstis method).
///
/// `K_IC = 0.016 * (E/H)^0.5 * P / c^1.5`
///
/// # Arguments
/// * `youngs_modulus` -- E (Pa).
/// * `hardness` -- Vickers hardness H (Pa).
/// * `load` -- indentation load P (N).
/// * `crack_length` -- radial crack half-length c (m).
///
/// Returns fracture toughness K_IC (Pa sqrt(m)).
pub fn indentation_fracture_toughness(
    youngs_modulus: f64,
    hardness: f64,
    load: f64,
    crack_length: f64,
) -> f64 {
    0.016 * (youngs_modulus / hardness).sqrt() * load / crack_length.powf(1.5)
}

// ---------------------------------------------------------------------------
// Ceramic fatigue (subcritical crack growth)
// ---------------------------------------------------------------------------

/// Subcritical crack growth rate (Paris-type law for ceramics).
///
/// `da/dN = A * (K_I / K_IC)^n`
///
/// where n is typically 10-50 for ceramics (much higher than metals).
///
/// # Arguments
/// * `k_i` -- stress intensity factor (Pa sqrt(m)).
/// * `k_ic` -- fracture toughness (Pa sqrt(m)).
/// * `paris_a` -- Paris law coefficient A (m/cycle).
/// * `paris_n` -- Paris law exponent n.
///
/// Returns crack growth rate (m/cycle).
pub fn subcritical_crack_growth_rate(k_i: f64, k_ic: f64, paris_a: f64, paris_n: f64) -> f64 {
    paris_a * (k_i / k_ic).powf(paris_n)
}

/// Estimate cycles to failure from initial crack size to critical crack size.
///
/// Uses numerical integration of the Paris law.
///
/// # Arguments
/// * `a0` -- initial crack half-length (m).
/// * `stress` -- applied stress (Pa).
/// * `geometry_factor` -- Y factor.
/// * `k_ic` -- fracture toughness (Pa sqrt(m)).
/// * `paris_a` -- Paris coefficient (m/cycle).
/// * `paris_n` -- Paris exponent.
/// * `steps` -- number of integration steps.
///
/// Returns estimated number of cycles to failure.
pub fn cycles_to_failure(
    a0: f64,
    stress: f64,
    geometry_factor: f64,
    k_ic: f64,
    paris_a: f64,
    paris_n: f64,
    steps: usize,
) -> f64 {
    let a_c = critical_flaw_size(k_ic, stress, geometry_factor);
    if a0 >= a_c {
        return 0.0;
    }
    let da = (a_c - a0) / steps as f64;
    let mut cycles = 0.0;
    let mut a = a0;
    for _ in 0..steps {
        let k_i = geometry_factor * stress * (PI * a).sqrt();
        let rate = subcritical_crack_growth_rate(k_i, k_ic, paris_a, paris_n);
        if rate < 1e-30 {
            cycles += da / 1e-30;
        } else {
            cycles += da / rate;
        }
        a += da;
    }
    cycles
}

// ---------------------------------------------------------------------------
// Ceramic R-curve behaviour
// ---------------------------------------------------------------------------

/// Compute rising R-curve toughness as function of crack extension.
///
/// `K_R(da) = K_0 + (K_ss - K_0) * (1 - exp(-da / lambda))`
///
/// where K_0 is the initiation toughness, K_ss the steady-state toughness,
/// and lambda is a characteristic bridging length.
///
/// # Arguments
/// * `da` -- crack extension (m).
/// * `k_init` -- initiation toughness K_0 (Pa sqrt(m)).
/// * `k_steady` -- steady-state toughness K_ss (Pa sqrt(m)).
/// * `lambda` -- characteristic bridging length (m).
///
/// Returns toughness at crack extension da (Pa sqrt(m)).
pub fn r_curve_toughness(da: f64, k_init: f64, k_steady: f64, lambda: f64) -> f64 {
    k_init + (k_steady - k_init) * (1.0 - (-da / lambda).exp())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // --- ceramic_preset ---

    #[test]
    fn test_preset_alumina_modulus() {
        let m = ceramic_preset(CeramicType::Alumina);
        assert!((m.youngs_modulus - 380.0e9).abs() < 1.0, "Alumina E");
    }

    #[test]
    fn test_preset_zirconia_toughness() {
        let m = ceramic_preset(CeramicType::ZirconiaYSZ);
        assert!(m.fracture_toughness_kic > 5.0e6, "YSZ KIC > 5 MPa sqrt(m)");
    }

    #[test]
    fn test_preset_sic_hardness() {
        let m = ceramic_preset(CeramicType::SiliconCarbide);
        assert!(m.hardness_vickers > 2000.0, "SiC HV > 2000");
    }

    #[test]
    fn test_preset_si3n4_conductivity() {
        let m = ceramic_preset(CeramicType::SiliconNitride);
        assert!(m.thermal_conductivity > 20.0, "Si3N4 k > 20 W/mK");
    }

    #[test]
    fn test_preset_b4c_melting_point() {
        let m = ceramic_preset(CeramicType::BoronCarbide);
        assert!(m.melting_point > 2500.0, "B4C melting > 2500 K");
    }

    #[test]
    fn test_preset_ha_modulus_low() {
        let m = ceramic_preset(CeramicType::HydroxyApatite);
        assert!(m.youngs_modulus < 120.0e9, "HA E < 120 GPa");
    }

    #[test]
    fn test_preset_all_positive_fields() {
        for ct in [
            CeramicType::Alumina,
            CeramicType::ZirconiaYSZ,
            CeramicType::SiliconCarbide,
            CeramicType::SiliconNitride,
            CeramicType::BoronCarbide,
            CeramicType::HydroxyApatite,
        ] {
            let m = ceramic_preset(ct);
            assert!(m.youngs_modulus > 0.0);
            assert!(m.fracture_toughness_kic > 0.0);
            assert!(m.hardness_vickers > 0.0);
            assert!(m.thermal_conductivity > 0.0);
            assert!(m.grain_size_um > 0.0);
        }
    }

    // --- hall_petch_strength ---

    #[test]
    fn test_hall_petch_basic() {
        let s = hall_petch_strength(0.5e6, 100.0e6, 1e-6);
        assert!((s - (100.0e6 + 0.5e6 / 1e-3)).abs() < 1.0, "HP value");
    }

    #[test]
    fn test_hall_petch_larger_grain_lower_strength() {
        let s_fine = hall_petch_strength(0.5e6, 100.0e6, 1e-6);
        let s_coarse = hall_petch_strength(0.5e6, 100.0e6, 100e-6);
        assert!(s_fine > s_coarse, "finer grain -> higher strength");
    }

    #[test]
    fn test_hall_petch_zero_k() {
        let s = hall_petch_strength(0.0, 200.0e6, 5e-6);
        assert!((s - 200.0e6).abs() < EPS, "k=0 -> sigma_0");
    }

    // --- critical_flaw_size ---

    #[test]
    fn test_critical_flaw_size_positive() {
        let a = critical_flaw_size(3.5e6, 200.0e6, 1.0);
        assert!(a > 0.0, "flaw size must be positive");
    }

    #[test]
    fn test_critical_flaw_size_decreases_with_stress() {
        let a1 = critical_flaw_size(3.5e6, 100.0e6, 1.0);
        let a2 = critical_flaw_size(3.5e6, 200.0e6, 1.0);
        assert!(a1 > a2, "higher stress -> smaller critical flaw");
    }

    #[test]
    fn test_critical_flaw_size_formula() {
        let kic = 3.5e6_f64;
        let sig = 200.0e6_f64;
        let y = 1.12_f64;
        let expected = (kic / (y * sig)).powi(2) / PI;
        let got = critical_flaw_size(kic, sig, y);
        assert!((got - expected).abs() < 1e-20, "formula");
    }

    // --- thermal_shock_resistance ---

    #[test]
    fn test_tsr_positive() {
        let m = ceramic_preset(CeramicType::Alumina);
        let r = thermal_shock_resistance(&m, 200.0);
        assert!(r > 0.0, "TSR must be positive");
    }

    #[test]
    fn test_tsr_scales_with_delta_t() {
        let m = ceramic_preset(CeramicType::SiliconCarbide);
        let r1 = thermal_shock_resistance(&m, 100.0);
        let r2 = thermal_shock_resistance(&m, 200.0);
        assert!((r2 - 2.0 * r1).abs() < EPS * 1e6, "TSR linear in delta_T");
    }

    // --- hasselman_parameters ---

    #[test]
    fn test_hasselman_r_first_positive() {
        let h = hasselman_parameters(400.0e6, 380.0e9, 0.22, 8.1e-6, 25.0, 3.5e6, 1.0);
        assert!(h.r_first > 0.0, "R must be positive");
    }

    #[test]
    fn test_hasselman_r_second_equals_k_times_r() {
        let k = 25.0;
        let h = hasselman_parameters(400.0e6, 380.0e9, 0.22, 8.1e-6, k, 3.5e6, 1.0);
        assert!((h.r_second - k * h.r_first).abs() < 1e-6, "R' = k * R");
    }

    #[test]
    fn test_hasselman_r_fourth_positive() {
        let h = hasselman_parameters(400.0e6, 380.0e9, 0.22, 8.1e-6, 25.0, 3.5e6, 1.0);
        assert!(h.r_fourth > 0.0, "R'''' must be positive");
    }

    // --- creep_rate_ceramic ---

    #[test]
    fn test_creep_rate_positive() {
        let rate = creep_rate_ceramic(1500.0, 100.0e6, 5e-6, 400_000.0);
        assert!(rate > 0.0, "creep rate must be positive");
    }

    #[test]
    fn test_creep_rate_higher_temp_faster() {
        let r1 = creep_rate_ceramic(1200.0, 100.0e6, 5e-6, 300_000.0);
        let r2 = creep_rate_ceramic(1600.0, 100.0e6, 5e-6, 300_000.0);
        assert!(r2 > r1, "higher T -> faster creep");
    }

    // --- sintering_densification ---

    #[test]
    fn test_sintering_density_increases_with_time() {
        let rho0 = 0.60;
        let d1 = sintering_densification(1600.0, 100.0, rho0);
        let d2 = sintering_densification(1600.0, 3600.0, rho0);
        assert!(d2 >= d1, "density should increase with time");
    }

    #[test]
    fn test_sintering_density_bounded() {
        let d = sintering_densification(1800.0, 1e9, 0.50);
        assert!(d <= 1.0, "density cannot exceed 1.0");
        assert!(d >= 0.50, "density cannot decrease from initial");
    }

    #[test]
    fn test_sintering_zero_time() {
        let rho0 = 0.65;
        let d = sintering_densification(1600.0, 0.0, rho0);
        assert!((d - rho0).abs() < 1e-12, "no sintering at t=0");
    }

    // --- grain_growth ---

    #[test]
    fn test_grain_growth_increases_with_time() {
        let d1 = grain_growth_size(1e-6, 2.0, 1e-14, 300e3, 1500.0, 100.0);
        let d2 = grain_growth_size(1e-6, 2.0, 1e-14, 300e3, 1500.0, 10000.0);
        assert!(d2 > d1, "grain size should increase with time");
    }

    #[test]
    fn test_grain_growth_zero_time_equals_initial() {
        let d0 = 1.5e-6;
        let d = grain_growth_size(d0, 2.0, 1e-14, 300e3, 1500.0, 0.0);
        assert!((d - d0).abs() < 1e-15, "d(t=0) = d0");
    }

    #[test]
    fn test_grain_growth_rate_positive() {
        let rate = grain_growth_rate(2e-6, 2.0, 1e-14, 300e3, 1500.0);
        assert!(rate > 0.0, "growth rate must be positive");
    }

    // --- vickers_hardness_to_mpa ---

    #[test]
    fn test_hv_conversion_known_value() {
        let mpa = vickers_hardness_to_mpa(1000.0);
        assert!((mpa - 9807.0).abs() < 0.01, "1000 HV = 9807 MPa");
    }

    // --- ceramic_fracture_probability ---

    #[test]
    fn test_weibull_at_scale_param() {
        let p = ceramic_fracture_probability(100.0e6, 100.0e6, 10.0);
        let expected = 1.0 - (-1.0_f64).exp();
        assert!((p - expected).abs() < 1e-12, "P_f at scale param");
    }

    #[test]
    fn test_weibull_zero_stress() {
        let p = ceramic_fracture_probability(0.0, 100.0e6, 10.0);
        assert!(p.abs() < EPS, "zero stress -> zero probability");
    }

    #[test]
    fn test_weibull_monotone_increasing() {
        let p1 = ceramic_fracture_probability(80.0e6, 100.0e6, 10.0);
        let p2 = ceramic_fracture_probability(120.0e6, 100.0e6, 10.0);
        assert!(p2 > p1, "P_f increases with stress");
    }

    // --- weibull_strength_at_probability ---

    #[test]
    fn test_weibull_strength_at_50_percent() {
        let sigma_0 = 400.0e6;
        let m = 10.0;
        let sigma = weibull_strength_at_probability(0.5, sigma_0, m);
        // Check round-trip
        let pf = ceramic_fracture_probability(sigma, sigma_0, m);
        assert!((pf - 0.5).abs() < 1e-10, "round-trip P_f = 0.5");
    }

    // --- weibull_volume_scaling ---

    #[test]
    fn test_volume_scaling_larger_weaker() {
        let sigma = weibull_volume_scaling(400.0e6, 1e-6, 1e-3, 10.0);
        assert!(sigma < 400.0e6, "larger volume -> lower strength");
    }

    #[test]
    fn test_volume_scaling_same_volume() {
        let sigma = weibull_volume_scaling(400.0e6, 1e-6, 1e-6, 10.0);
        assert!(
            (sigma - 400.0e6).abs() < EPS,
            "same volume -> same strength"
        );
    }

    // --- thermal_conductivity models ---

    #[test]
    fn test_debye_conductivity_positive() {
        let k = debye_thermal_conductivity(2.5e6, 5000.0, 10e-9);
        assert!(k > 0.0, "conductivity must be positive");
    }

    #[test]
    fn test_debye_conductivity_proportional_to_mfp() {
        let k1 = debye_thermal_conductivity(2.5e6, 5000.0, 10e-9);
        let k2 = debye_thermal_conductivity(2.5e6, 5000.0, 20e-9);
        assert!((k2 / k1 - 2.0).abs() < 1e-10, "linear in mfp");
    }

    #[test]
    fn test_porous_conductivity_decreases_with_porosity() {
        let k0 = porous_thermal_conductivity(25.0, 0.0);
        let k1 = porous_thermal_conductivity(25.0, 0.3);
        assert!(k1 < k0, "porosity reduces conductivity");
    }

    #[test]
    fn test_porous_conductivity_zero_porosity_equals_solid() {
        let k = porous_thermal_conductivity(25.0, 0.0);
        assert!((k - 25.0).abs() < 1e-10, "zero porosity -> solid value");
    }

    #[test]
    fn test_umklapp_conductivity_inverse_t() {
        let k1 = umklapp_conductivity(25.0, 300.0, 300.0);
        let k2 = umklapp_conductivity(25.0, 300.0, 600.0);
        assert!((k1 - 25.0).abs() < 1e-10, "at T_ref -> k_ref");
        assert!((k2 - 12.5).abs() < 1e-10, "at 2*T_ref -> k_ref/2");
    }

    // --- piezoelectric ceramics ---

    #[test]
    fn test_pzt4_d33_positive() {
        let pzt = pzt4_preset();
        assert!(pzt.d33 > 0.0, "d33 must be positive for PZT");
    }

    #[test]
    fn test_pzt5h_higher_d33_than_pzt4() {
        let p4 = pzt4_preset();
        let p5h = pzt5h_preset();
        assert!(p5h.d33 > p4.d33, "PZT-5H is softer -> higher d33");
    }

    #[test]
    fn test_piezo_strain_sign() {
        let pzt = pzt5a_preset();
        let s33 = piezo_strain_33(&pzt, 1e6);
        let s31 = piezo_strain_31(&pzt, 1e6);
        assert!(s33 > 0.0, "d33 > 0 -> positive axial strain");
        assert!(s31 < 0.0, "d31 < 0 -> negative transverse strain");
    }

    #[test]
    fn test_piezo_voltage_proportional_to_stress() {
        let pzt = pzt4_preset();
        let v1 = piezo_open_circuit_voltage(&pzt, 1.0e6, 1e-3);
        let v2 = piezo_open_circuit_voltage(&pzt, 2.0e6, 1e-3);
        assert!((v2 / v1 - 2.0).abs() < 1e-10, "V linear in stress");
    }

    #[test]
    fn test_piezo_charge_proportional_to_area() {
        let pzt = pzt4_preset();
        let q1 = piezo_charge_generated(&pzt, 1e6, 1e-4);
        let q2 = piezo_charge_generated(&pzt, 1e6, 2e-4);
        assert!((q2 / q1 - 2.0).abs() < 1e-10, "charge linear in area");
    }

    #[test]
    fn test_piezo_energy_density_positive() {
        let pzt = pzt5a_preset();
        let u = piezo_energy_density(&pzt, 1e6);
        assert!(u > 0.0, "energy density must be positive");
    }

    #[test]
    fn test_piezo_resonant_frequency() {
        let f = piezo_resonant_frequency(4000.0, 1e-3);
        assert!((f - 2.0e6).abs() < 1.0, "f_r = v/(2t) = 4000/0.002 = 2 MHz");
    }

    // --- ferroelectric hysteresis ---

    #[test]
    fn test_hysteresis_upper_at_zero_field() {
        let h = FerroelectricHysteresis::barium_titanate();
        let p = h.polarisation_upper(0.0);
        // At E=0, upper branch: tanh(E_c / E_a) > 0
        assert!(p > 0.0, "upper branch at E=0 should be positive (remanent)");
    }

    #[test]
    fn test_hysteresis_lower_at_zero_field() {
        let h = FerroelectricHysteresis::barium_titanate();
        let p = h.polarisation_lower(0.0);
        // At E=0, lower branch: tanh(-E_c / E_a) < 0
        assert!(p < 0.0, "lower branch at E=0 should be negative");
    }

    #[test]
    fn test_hysteresis_saturation() {
        let h = FerroelectricHysteresis::pzt_soft();
        let p_up = h.polarisation_upper(1e7);
        // At very high field, should approach P_s
        assert!(
            (p_up - h.saturation_polarisation).abs() < 0.01,
            "should saturate at high field"
        );
    }

    #[test]
    fn test_hysteresis_loss_positive() {
        let h = FerroelectricHysteresis::barium_titanate();
        assert!(h.hysteresis_loss_per_cycle() > 0.0, "loss must be positive");
    }

    #[test]
    fn test_curie_weiss_diverges_at_tc() {
        let eps = FerroelectricHysteresis::curie_weiss_permittivity(1.7e5, 393.0, 393.0);
        assert!(eps.is_infinite(), "should diverge at T_c");
    }

    #[test]
    fn test_curie_weiss_decreases_above_tc() {
        let eps1 = FerroelectricHysteresis::curie_weiss_permittivity(1.7e5, 393.0, 400.0);
        let eps2 = FerroelectricHysteresis::curie_weiss_permittivity(1.7e5, 393.0, 500.0);
        assert!(
            eps1 > eps2,
            "permittivity decreases as T increases above T_c"
        );
    }

    // --- zirconia toughening ---

    #[test]
    fn test_zirconia_toughening_increment_positive() {
        let z = ZirconiaToughening::y_tzp_3mol();
        assert!(z.toughening_increment() > 0.0, "delta_K must be positive");
    }

    #[test]
    fn test_zirconia_effective_gt_base() {
        let z = ZirconiaToughening::y_tzp_3mol();
        assert!(z.effective_toughness() > z.base_toughness, "K_eff > K_0");
    }

    #[test]
    fn test_zirconia_toughening_increases_with_fv() {
        let z1 = ZirconiaToughening::new(0.04, 0.1, 80e9, 5e6, 50e-6);
        let z2 = ZirconiaToughening::new(0.04, 0.4, 80e9, 5e6, 50e-6);
        assert!(
            z2.toughening_increment() > z1.toughening_increment(),
            "more transformable phase -> more toughening"
        );
    }

    #[test]
    fn test_zirconia_zone_size_positive() {
        let z = ZirconiaToughening::y_tzp_3mol();
        assert!(z.zone_size_estimate() > 0.0, "zone size must be positive");
    }

    #[test]
    fn test_zirconia_transformation_stress_positive() {
        let z = ZirconiaToughening::y_tzp_3mol();
        assert!(z.transformation_stress() > 0.0, "sigma_t must be positive");
    }

    #[test]
    fn test_zirconia_critical_grain_size_positive() {
        let z = ZirconiaToughening::y_tzp_3mol();
        let dc = z.critical_grain_size(1.0, 100e6, 210e9, 0.31);
        assert!(dc > 0.0, "critical grain size must be positive");
    }

    // --- alumina models ---

    #[test]
    fn test_alumina_modulus_decreases_with_temp() {
        let e_rt = alumina_modulus_vs_temp(380e9, 300.0, 300.0, 5e-5);
        let e_ht = alumina_modulus_vs_temp(380e9, 1000.0, 300.0, 5e-5);
        assert!(e_ht < e_rt, "modulus should decrease at high temperature");
    }

    #[test]
    fn test_alumina_modulus_at_ref_temp() {
        let e = alumina_modulus_vs_temp(380e9, 300.0, 300.0, 5e-5);
        assert!((e - 380e9).abs() < 1.0, "E(T_ref) = E_0");
    }

    // --- indentation fracture toughness ---

    #[test]
    fn test_indentation_toughness_positive() {
        let k = indentation_fracture_toughness(380e9, 18e9, 10.0, 50e-6);
        assert!(k > 0.0, "K_IC must be positive");
    }

    // --- subcritical crack growth ---

    #[test]
    fn test_scg_rate_increases_with_ki() {
        let r1 = subcritical_crack_growth_rate(2e6, 4e6, 1e-15, 20.0);
        let r2 = subcritical_crack_growth_rate(3e6, 4e6, 1e-15, 20.0);
        assert!(r2 > r1, "higher K_I -> faster growth");
    }

    #[test]
    fn test_cycles_to_failure_positive() {
        let n = cycles_to_failure(10e-6, 200e6, 1.12, 3.5e6, 1e-15, 20.0, 1000);
        assert!(n > 0.0, "cycles must be positive for subcritical crack");
    }

    // --- R-curve ---

    #[test]
    fn test_r_curve_at_zero_extension_equals_k0() {
        let k = r_curve_toughness(0.0, 3.0e6, 8.0e6, 100e-6);
        assert!((k - 3.0e6).abs() < 1e-6, "K_R(0) = K_0");
    }

    #[test]
    fn test_r_curve_approaches_steady_state() {
        let k = r_curve_toughness(1.0, 3.0e6, 8.0e6, 100e-6);
        assert!((k - 8.0e6).abs() < 1.0, "K_R at large extension -> K_ss");
    }

    #[test]
    fn test_r_curve_monotonically_increasing() {
        let k1 = r_curve_toughness(10e-6, 3.0e6, 8.0e6, 100e-6);
        let k2 = r_curve_toughness(50e-6, 3.0e6, 8.0e6, 100e-6);
        assert!(k2 > k1, "R-curve should be monotonically increasing");
    }
}
