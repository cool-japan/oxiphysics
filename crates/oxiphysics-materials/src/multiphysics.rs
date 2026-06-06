// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coupled multiphysics material models.
//!
//! Provides constitutive models coupling mechanical, thermal, electrical,
//! magnetic, electrochemical, and radiation physics including operator
//! splitting and staggered Newton solvers.

// ─────────────────────────────────────────────────────────────────────────────
// ThermoMechanical — thermal stress coupling
// ─────────────────────────────────────────────────────────────────────────────

/// Thermo-mechanical coupling: thermal stress and heat generation from deformation.
///
/// σ = C : (ε - α * ΔT * I), coupled with ρ*Cp*dT/dt = k*∇²T + σ:dε/dt.
#[derive(Debug, Clone)]
pub struct ThermoMechanical {
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Coefficient of thermal expansion α (1/K).
    pub cte: f64,
    /// Thermal conductivity k (W/m/K).
    pub thermal_conductivity: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Specific heat capacity Cp (J/kg/K).
    pub specific_heat: f64,
    /// Reference temperature T0 (K).
    pub t_ref: f64,
}

impl ThermoMechanical {
    /// Creates a thermo-mechanical model for steel.
    pub fn steel() -> Self {
        Self {
            youngs_modulus: 200e9,
            poisson_ratio: 0.3,
            cte: 12e-6,
            thermal_conductivity: 50.0,
            density: 7850.0,
            specific_heat: 500.0,
            t_ref: 293.15,
        }
    }

    /// Creates a custom thermo-mechanical material.
    pub fn new(
        youngs_modulus: f64,
        poisson_ratio: f64,
        cte: f64,
        thermal_conductivity: f64,
        density: f64,
        specific_heat: f64,
        t_ref: f64,
    ) -> Self {
        Self {
            youngs_modulus,
            poisson_ratio,
            cte,
            thermal_conductivity,
            density,
            specific_heat,
            t_ref,
        }
    }

    /// Thermal strain: ε_th = α * ΔT.
    pub fn thermal_strain(&self, temperature: f64) -> f64 {
        self.cte * (temperature - self.t_ref)
    }

    /// Uniaxial thermal stress: σ = E * α * ΔT (constrained expansion).
    pub fn thermal_stress_uniaxial(&self, temperature: f64) -> f64 {
        -self.youngs_modulus * self.thermal_strain(temperature)
    }

    /// Thermal diffusivity α_d = k / (ρ * Cp) (m²/s).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.thermal_conductivity / (self.density * self.specific_heat)
    }

    /// Lamé parameter λ.
    pub fn lame_lambda(&self) -> f64 {
        self.youngs_modulus * self.poisson_ratio
            / ((1.0 + self.poisson_ratio) * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Lamé parameter μ (shear modulus).
    pub fn lame_mu(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Hydrostatic thermal stress (3D isotropic): σ_v = -3K*α*ΔT.
    pub fn hydrostatic_thermal_stress(&self, temperature: f64) -> f64 {
        let k_bulk = self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio));
        -3.0 * k_bulk * self.cte * (temperature - self.t_ref)
    }

    /// Thermoelastic dissipation rate per unit volume (Gough-Joule effect).
    pub fn thermoelastic_dissipation(&self, strain_rate: f64, temperature: f64) -> f64 {
        let k_bulk = self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio));
        -3.0 * k_bulk * self.cte * temperature * strain_rate
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PiezoElectric — piezoelectric coupling
// ─────────────────────────────────────────────────────────────────────────────

/// Piezoelectric constitutive model.
///
/// Couples mechanical strain and electric displacement via the d-matrix:
/// S = sE : T + d^t : E_field;  D = d : T + ε^T : E_field.
#[derive(Debug, Clone)]
pub struct PiezoElectric {
    /// Mechanical compliance at constant E-field \[s11, s12, s33, s44\] (Pa^-1).
    pub compliance: [f64; 4],
    /// Piezoelectric d-matrix coefficients \[d31, d33, d15\] (C/N = m/V).
    pub d_matrix: [f64; 3],
    /// Permittivity at constant stress \[ε11, ε33\] (F/m).
    pub permittivity: [f64; 2],
    /// Material name.
    pub name: String,
}

impl PiezoElectric {
    /// Creates a model for PZT-5H (lead zirconate titanate).
    pub fn pzt5h() -> Self {
        Self {
            compliance: [16.5e-12, -4.78e-12, 20.7e-12, 43.5e-12],
            d_matrix: [-274e-12, 593e-12, 741e-12],
            permittivity: [3130.0 * 8.854e-12, 3400.0 * 8.854e-12],
            name: "PZT-5H".to_string(),
        }
    }

    /// Creates a model for PVDF (polyvinylidene fluoride).
    pub fn pvdf() -> Self {
        Self {
            compliance: [365e-12, -145e-12, 472e-12, 1600e-12],
            d_matrix: [23e-12, -33e-12, -27e-12],
            permittivity: [12.0 * 8.854e-12, 12.0 * 8.854e-12],
            name: "PVDF".to_string(),
        }
    }

    /// Direct piezoelectric effect: charge density D3 = d33 * T3.
    pub fn direct_effect_d33(&self, stress_33: f64) -> f64 {
        self.d_matrix[1] * stress_33
    }

    /// Converse piezoelectric effect: strain S3 = d33 * E3.
    pub fn converse_effect_d33(&self, e_field_3: f64) -> f64 {
        self.d_matrix[1] * e_field_3
    }

    /// Electromechanical coupling coefficient k33.
    pub fn coupling_k33(&self) -> f64 {
        let d33 = self.d_matrix[1];
        let s33 = self.compliance[2];
        let eps33 = self.permittivity[1];
        d33 / (s33 * eps33).sqrt()
    }

    /// Open-circuit resonance shift due to piezoelectric stiffening.
    pub fn stiffened_modulus_33(&self) -> f64 {
        let d33 = self.d_matrix[1];
        let s33 = self.compliance[2];
        let eps33 = self.permittivity[1];
        (s33 - d33 * d33 / eps33).recip()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Magnetostrictive — Joule magnetostriction and Villari effect
// ─────────────────────────────────────────────────────────────────────────────

/// Magnetostrictive material model (Joule and Villari effects).
///
/// Models Terfenol-D-like materials with magnetization curve and Preisach hysteresis.
#[derive(Debug, Clone)]
pub struct Magnetostrictive {
    /// Saturation magnetostriction λ_s (dimensionless, e.g. 1750e-6 for Terfenol-D).
    pub lambda_s: f64,
    /// Saturation magnetization M_s (A/m).
    pub m_saturation: f64,
    /// Young's modulus at H=0 (Pa).
    pub youngs_modulus: f64,
    /// Magnetomechanical coupling coefficient d (m/A).
    pub d_coeff: f64,
    /// Coercive field Hc (A/m) for hysteresis.
    pub coercive_field: f64,
}

impl Magnetostrictive {
    /// Creates a Terfenol-D model.
    pub fn terfenol_d() -> Self {
        Self {
            lambda_s: 1750e-6,
            m_saturation: 7.65e5,
            youngs_modulus: 30e9,
            d_coeff: 1.67e-8,
            coercive_field: 10e3,
        }
    }

    /// Joule magnetostriction: ε = (3/2) * λ_s * (M/M_s)^2.
    pub fn magnetostriction(&self, magnetization: f64) -> f64 {
        let m_norm = (magnetization / self.m_saturation).clamp(-1.0, 1.0);
        1.5 * self.lambda_s * m_norm * m_norm
    }

    /// Villari effect: change in magnetization due to stress.
    ///
    /// dM/dσ ≈ d * H_applied at constant field.
    pub fn villari_effect(&self, stress: f64, h_field: f64) -> f64 {
        // Simplified: dM ≈ d * σ * (dB/dH) at operating point
        self.d_coeff * stress * h_field / self.coercive_field
    }

    /// Preisach-like hysteresis: simple scalar Preisach model.
    ///
    /// Returns the magnetization for a triangular hysteron with thresholds (α, β).
    pub fn preisach_hysteron(&self, h: f64, alpha: f64, beta: f64, state: bool) -> bool {
        if h > alpha {
            true
        } else if h < beta {
            false
        } else {
            state
        }
    }

    /// Magnetization curve (Langevin function approximation).
    pub fn magnetization_curve(&self, h_field: f64) -> f64 {
        let a = self.coercive_field;
        let x = h_field / a;
        // Brillouin/Langevin approximation
        self.m_saturation * x.tanh()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ElectrochemicalMaterial — battery electrode material
// ─────────────────────────────────────────────────────────────────────────────

/// Electrochemical material model for battery electrodes.
///
/// Butler-Volmer kinetics, diffusion, and intercalation stress.
#[derive(Debug, Clone)]
pub struct ElectrochemicalMaterial {
    /// Exchange current density i0 (A/m²).
    pub exchange_current_density: f64,
    /// Charge transfer coefficient α (0 to 1).
    pub alpha: f64,
    /// Diffusion coefficient D (m²/s).
    pub diffusion_coefficient: f64,
    /// Partial molar volume Ω (m³/mol).
    pub partial_molar_volume: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Faraday constant.
    pub faraday: f64,
}

impl ElectrochemicalMaterial {
    /// Faraday constant (C/mol).
    const F: f64 = 96485.0;
    /// Gas constant (J/mol/K).
    const R: f64 = 8.314;

    /// Creates a lithium graphite anode model.
    pub fn lithium_graphite() -> Self {
        Self {
            exchange_current_density: 2.0,
            alpha: 0.5,
            diffusion_coefficient: 3.9e-14,
            partial_molar_volume: 3.497e-6,
            youngs_modulus: 10e9,
            poisson_ratio: 0.3,
            faraday: Self::F,
        }
    }

    /// Butler-Volmer current density: j = i0 * \[exp(α*F*η/RT) - exp(-(1-α)*F*η/RT)\].
    pub fn butler_volmer(&self, overpotential: f64, temperature: f64) -> f64 {
        let f_rt = Self::F / (Self::R * temperature);
        self.exchange_current_density
            * ((self.alpha * f_rt * overpotential).exp()
                - (-(1.0 - self.alpha) * f_rt * overpotential).exp())
    }

    /// Linear diffusion flux (1D, Fick's first law): J = -D * dc/dx.
    pub fn diffusion_flux(&self, concentration_gradient: f64) -> f64 {
        -self.diffusion_coefficient * concentration_gradient
    }

    /// Intercalation stress (hydrostatic): σ_h = -Ω * E / (3 * (1-ν)) * (c - c_avg).
    pub fn intercalation_stress(&self, concentration: f64, c_avg: f64) -> f64 {
        let prefactor =
            self.partial_molar_volume * self.youngs_modulus / (3.0 * (1.0 - self.poisson_ratio));
        -prefactor * (concentration - c_avg)
    }

    /// Open-circuit potential shift due to concentration (Nernst-like).
    pub fn nernst_shift(&self, c_norm: f64, temperature: f64) -> f64 {
        let c = c_norm.clamp(1e-10, 1.0 - 1e-10);
        Self::R * temperature / Self::F * (c / (1.0 - c)).ln()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PoroElastic — Biot consolidation
// ─────────────────────────────────────────────────────────────────────────────

/// Poroelastic material: Biot consolidation theory.
///
/// Couples effective stress with pore pressure: σ_eff = σ_total - b * p * I.
#[derive(Debug, Clone)]
pub struct PoroElastic {
    /// Drained Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν (drained).
    pub poisson_ratio: f64,
    /// Biot coefficient b (0 to 1).
    pub biot_coefficient: f64,
    /// Permeability k (m²).
    pub permeability: f64,
    /// Fluid viscosity μ (Pa·s).
    pub fluid_viscosity: f64,
    /// Biot modulus M (Pa).
    pub biot_modulus: f64,
}

impl PoroElastic {
    /// Creates a saturated clay model.
    pub fn saturated_clay() -> Self {
        Self {
            youngs_modulus: 20e6,
            poisson_ratio: 0.35,
            biot_coefficient: 0.9,
            permeability: 1e-13,
            fluid_viscosity: 1e-3,
            biot_modulus: 1e9,
        }
    }

    /// Creates a sandstone model.
    pub fn sandstone() -> Self {
        Self {
            youngs_modulus: 15e9,
            poisson_ratio: 0.25,
            biot_coefficient: 0.7,
            permeability: 1e-13,
            fluid_viscosity: 1e-3,
            biot_modulus: 20e9,
        }
    }

    /// Consolidation coefficient cv = k * M / μ (m²/s).
    pub fn consolidation_coefficient(&self) -> f64 {
        self.permeability * self.biot_modulus / self.fluid_viscosity
    }

    /// Effective stress: σ_eff = σ_total - b * p.
    pub fn effective_stress(&self, total_stress: f64, pore_pressure: f64) -> f64 {
        total_stress - self.biot_coefficient * pore_pressure
    }

    /// Darcy velocity: q = -(k/μ) * ∇p.
    pub fn darcy_velocity(&self, pressure_gradient: f64) -> f64 {
        -(self.permeability / self.fluid_viscosity) * pressure_gradient
    }

    /// Skempton coefficient B = b / (b + φ * M * (1/Kf - 1/Ks)).
    ///
    /// Approximate: B ≈ b*M / (K_drained + b²*M).
    pub fn skempton_coefficient(&self) -> f64 {
        let k_bulk = self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio));
        let b = self.biot_coefficient;
        let m = self.biot_modulus;
        b * m / (k_bulk + b * b * m)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SwellingMaterial — hygroscopic and gel swelling
// ─────────────────────────────────────────────────────────────────────────────

/// Swelling material: hygroscopic expansion, osmotic pressure, gel (Flory-Rehner).
#[derive(Debug, Clone)]
pub struct SwellingMaterial {
    /// Hygroscopic expansion coefficient β_h (strain per unit relative humidity).
    pub hygroscopic_coeff: f64,
    /// Osmotic pressure coefficient (Pa per mol/m³ concentration difference).
    pub osmotic_coeff: f64,
    /// Flory-Rehner cross-link density (mol/m³).
    pub crosslink_density: f64,
    /// Flory-Huggins parameter χ.
    pub flory_chi: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl SwellingMaterial {
    /// Gas constant (J/mol/K).
    const R: f64 = 8.314;

    /// Creates a polyacrylamide hydrogel model.
    pub fn hydrogel() -> Self {
        Self {
            hygroscopic_coeff: 0.002,
            osmotic_coeff: 2479.0, // RT at 25°C = 2479 J/mol
            crosslink_density: 100.0,
            flory_chi: 0.5,
            temperature: 298.15,
        }
    }

    /// Hygroscopic swelling strain: ε = β_h * ΔRH.
    pub fn hygroscopic_strain(&self, delta_rh: f64) -> f64 {
        self.hygroscopic_coeff * delta_rh
    }

    /// Osmotic pressure (van't Hoff): Π = c * R * T.
    pub fn osmotic_pressure(&self, concentration: f64) -> f64 {
        concentration * Self::R * self.temperature
    }

    /// Flory-Rehner elastic free energy: ΔF_el = (3/2) * ν * kT * (λ² - 1 - ln λ).
    ///
    /// λ = swelling ratio.
    pub fn flory_rehner_elastic_energy(&self, swelling_ratio: f64) -> f64 {
        let nu = self.crosslink_density;
        let kt = Self::R * self.temperature / 6.022e23;
        1.5 * nu * kt * (swelling_ratio * swelling_ratio - 1.0 - swelling_ratio.ln())
    }

    /// Total swelling driving force: mixing + elastic.
    pub fn swelling_equilibrium_condition(&self, phi: f64) -> f64 {
        // dF_mix/dφ + dF_el/dφ = 0 at equilibrium
        let phi = phi.clamp(1e-5, 1.0 - 1e-5);
        let mixing = (1.0 - phi).ln() + phi + self.flory_chi * phi * phi;
        let elastic = self.crosslink_density * (phi.powf(1.0 / 3.0) - 0.5 * phi);
        mixing + elastic
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RadiationDamage — displacement per atom and void swelling
// ─────────────────────────────────────────────────────────────────────────────

/// Radiation damage material model.
///
/// Tracks displacement per atom (dpa), void swelling, embrittlement, and He bubble growth.
#[derive(Debug, Clone)]
pub struct RadiationDamage {
    /// Displacement threshold energy Ed (eV).
    pub threshold_energy_ev: f64,
    /// Atomic density n (atoms/m³).
    pub atomic_density: f64,
    /// Initial yield strength σy0 (Pa).
    pub yield_strength_initial: f64,
    /// Accumulated dpa (displacement per atom).
    pub dpa: f64,
    /// Void swelling fraction (dimensionless).
    pub void_fraction: f64,
    /// Helium concentration (appm).
    pub helium_conc_appm: f64,
}

impl RadiationDamage {
    /// Creates a stainless steel 316L radiation damage model.
    pub fn ss316l() -> Self {
        Self {
            threshold_energy_ev: 40.0,
            atomic_density: 8.5e28,
            yield_strength_initial: 220e6,
            dpa: 0.0,
            void_fraction: 0.0,
            helium_conc_appm: 0.0,
        }
    }

    /// Irradiation dose in dpa from neutron flux and cross-section.
    ///
    /// dpa = φ * σ_d * t / (n * Ed * 2).
    pub fn compute_dpa(&self, flux: f64, cross_section: f64, time: f64) -> f64 {
        flux * cross_section * time
            / (self.atomic_density * 2.0 * self.threshold_energy_ev * 1.6e-19)
    }

    /// Yield strength increase due to irradiation hardening.
    ///
    /// Δσy ≈ α * μ * b * sqrt(N_d), where N_d ~ dpa density (simplified to Δσy = k * sqrt(dpa)).
    pub fn irradiation_hardening(&self, k_hardening: f64) -> f64 {
        k_hardening * self.dpa.sqrt()
    }

    /// Total yield strength (initial + hardening).
    pub fn yield_strength(&self, k_hardening: f64) -> f64 {
        self.yield_strength_initial + self.irradiation_hardening(k_hardening)
    }

    /// Void swelling rate: dS/ddpa ≈ A * exp(-B/T) (temperature-dependent).
    pub fn void_swelling_rate(&self, temperature_k: f64, a_coeff: f64, b_coeff: f64) -> f64 {
        a_coeff * (-b_coeff / temperature_k).exp()
    }

    /// Helium bubble growth: r ~ (3 * C_He / (4π * n_b))^(1/3).
    pub fn bubble_radius(&self, bubble_density: f64) -> f64 {
        let c_he = self.helium_conc_appm * 1e-6 * self.atomic_density;
        (3.0 * c_he / (4.0 * std::f64::consts::PI * bubble_density)).cbrt()
    }

    /// Ductile-to-brittle transition temperature shift (DBTT shift ~ 20°C/dpa).
    pub fn dbtt_shift_k(&self) -> f64 {
        20.0 * self.dpa
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PhaseTransformMaterial — shape memory alloy
// ─────────────────────────────────────────────────────────────────────────────

/// Shape memory alloy (SMA) phase transformation material.
///
/// Implements the Tanaka model for martensite fraction and pseudoelasticity.
#[derive(Debug, Clone)]
pub struct PhaseTransformMaterial {
    /// Young's modulus of austenite EA (Pa).
    pub ea: f64,
    /// Young's modulus of martensite EM (Pa).
    pub em: f64,
    /// Maximum transformation strain εL.
    pub max_strain: f64,
    /// Austenite start temperature As (K).
    pub as_temp: f64,
    /// Austenite finish temperature Af (K).
    pub af_temp: f64,
    /// Martensite start temperature Ms (K).
    pub ms_temp: f64,
    /// Martensite finish temperature Mf (K).
    pub mf_temp: f64,
    /// Current martensite volume fraction ξ ∈ \[0, 1\].
    pub xi: f64,
}

impl PhaseTransformMaterial {
    /// Creates a Nitinol (NiTi) SMA model.
    pub fn nitinol() -> Self {
        Self {
            ea: 75e9,
            em: 30e9,
            max_strain: 0.08,
            as_temp: 291.0,
            af_temp: 307.0,
            ms_temp: 291.0,
            mf_temp: 271.0,
            xi: 0.0,
        }
    }

    /// Effective Young's modulus (mixture rule): E = EA + ξ*(EM - EA).
    pub fn effective_modulus(&self) -> f64 {
        self.ea + self.xi * (self.em - self.ea)
    }

    /// Martensite fraction from temperature (cooling, Tanaka model).
    pub fn martensite_fraction_cooling(&self, temperature: f64) -> f64 {
        if temperature >= self.ms_temp {
            0.0
        } else if temperature <= self.mf_temp {
            1.0
        } else {
            let b_m = std::f64::consts::PI / (self.ms_temp - self.mf_temp);
            0.5 * (1.0 - (b_m * (temperature - (self.ms_temp + self.mf_temp) / 2.0)).cos())
        }
    }

    /// Martensite fraction from temperature (heating, Tanaka model).
    pub fn martensite_fraction_heating(&self, temperature: f64) -> f64 {
        if temperature >= self.af_temp {
            0.0
        } else if temperature <= self.as_temp {
            self.xi // no change if below As
        } else {
            let b_a = std::f64::consts::PI / (self.af_temp - self.as_temp);
            self.xi * 0.5 * (1.0 + (b_a * (temperature - self.as_temp)).cos())
        }
    }

    /// Pseudoelastic transformation stress (at constant T).
    pub fn transformation_stress(&self, strain: f64) -> f64 {
        let eps = strain.clamp(0.0, self.max_strain);
        let xi_new = eps / self.max_strain;
        let e_eff = self.ea + xi_new * (self.em - self.ea);
        e_eff * (eps - xi_new * self.max_strain)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ElectroMagnetic — electromagnetic material properties
// ─────────────────────────────────────────────────────────────────────────────

/// Electromagnetic material properties: permeability, losses, skin depth.
#[derive(Debug, Clone)]
pub struct ElectroMagnetic {
    /// Relative permeability μr.
    pub rel_permeability: f64,
    /// Electrical conductivity σ (S/m).
    pub conductivity: f64,
    /// Steinmetz coefficient k_h for hysteresis loss.
    pub steinmetz_k_h: f64,
    /// Steinmetz exponent n.
    pub steinmetz_n: f64,
    /// Eddy current coefficient k_e.
    pub eddy_current_k_e: f64,
    /// Operating frequency (Hz).
    pub frequency: f64,
}

impl ElectroMagnetic {
    /// Permeability of free space μ0 (H/m).
    const MU0: f64 = 1.256_637_062e-6;

    /// Creates a silicon steel model (electrical steel).
    pub fn silicon_steel() -> Self {
        Self {
            rel_permeability: 5000.0,
            conductivity: 2e6,
            steinmetz_k_h: 40.0,
            steinmetz_n: 1.8,
            eddy_current_k_e: 0.5,
            frequency: 50.0,
        }
    }

    /// Absolute permeability μ = μ0 * μr (H/m).
    pub fn permeability(&self) -> f64 {
        Self::MU0 * self.rel_permeability
    }

    /// Skin depth δ = sqrt(2 / (ω * σ * μ)) (m).
    pub fn skin_depth(&self) -> f64 {
        let omega = 2.0 * std::f64::consts::PI * self.frequency;
        (2.0 / (omega * self.conductivity * self.permeability())).sqrt()
    }

    /// Hysteresis loss (Steinmetz): P_h = k_h * f * B^n (W/m³ per cycle normalized).
    pub fn hysteresis_loss(&self, b_max: f64) -> f64 {
        self.steinmetz_k_h * self.frequency * b_max.powf(self.steinmetz_n)
    }

    /// Eddy current loss: P_e = k_e * f² * B² (W/m³).
    pub fn eddy_current_loss(&self, b_max: f64) -> f64 {
        self.eddy_current_k_e * self.frequency.powi(2) * b_max.powi(2)
    }

    /// Total core loss: P_total = P_h + P_e.
    pub fn total_core_loss(&self, b_max: f64) -> f64 {
        self.hysteresis_loss(b_max) + self.eddy_current_loss(b_max)
    }

    /// Inductance per unit length for a solenoid (H/m).
    pub fn inductance_per_length(&self, n_turns_per_m: f64, area_m2: f64) -> f64 {
        self.permeability() * n_turns_per_m * n_turns_per_m * area_m2
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CoupledSolver — multiphysics coupling strategies
// ─────────────────────────────────────────────────────────────────────────────

/// Multiphysics coupled solver strategies.
///
/// Provides operator splitting, staggered scheme, and monolithic Newton iteration.
pub struct CoupledSolver {
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Relaxation parameter ω for successive substitution (0 < ω ≤ 1).
    pub relaxation: f64,
}

impl CoupledSolver {
    /// Creates a new coupled solver.
    pub fn new(tolerance: f64, max_iter: usize, relaxation: f64) -> Self {
        Self {
            tolerance,
            max_iter,
            relaxation,
        }
    }

    /// Operator splitting: alternately solve field A then field B.
    ///
    /// `solve_a(b) -> a` and `solve_b(a) -> b` are the individual field solvers.
    /// Returns the converged (a, b) pair and number of iterations.
    pub fn operator_split<FA, FB>(
        &self,
        a0: f64,
        b0: f64,
        solve_a: FA,
        solve_b: FB,
    ) -> (f64, f64, usize)
    where
        FA: Fn(f64) -> f64,
        FB: Fn(f64) -> f64,
    {
        let mut a = a0;
        let mut b = b0;
        for iter in 0..self.max_iter {
            let a_new = solve_a(b);
            let b_new = solve_b(a_new);
            let err = ((a_new - a).powi(2) + (b_new - b).powi(2)).sqrt();
            a = a * (1.0 - self.relaxation) + a_new * self.relaxation;
            b = b * (1.0 - self.relaxation) + b_new * self.relaxation;
            if err < self.tolerance {
                return (a, b, iter + 1);
            }
        }
        (a, b, self.max_iter)
    }

    /// Staggered scheme: sequential solve with fixed-point iteration.
    ///
    /// Returns converged state vector and iteration count.
    pub fn staggered<F>(&self, state0: &[f64], step_fn: F) -> (Vec<f64>, usize)
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let mut state = state0.to_vec();
        for iter in 0..self.max_iter {
            let state_new = step_fn(&state);
            let err = state
                .iter()
                .zip(state_new.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            // Relaxed update
            let state_relaxed: Vec<f64> = state
                .iter()
                .zip(state_new.iter())
                .map(|(a, b)| a * (1.0 - self.relaxation) + b * self.relaxation)
                .collect();
            state = state_relaxed;
            if err < self.tolerance {
                return (state, iter + 1);
            }
        }
        (state, self.max_iter)
    }

    /// Monolithic Newton iteration for a coupled residual R(x) = 0.
    ///
    /// Uses finite-difference Jacobian and direct solve via Gaussian elimination.
    pub fn monolithic_newton<F>(&self, x0: &[f64], residual: F) -> (Vec<f64>, usize)
    where
        F: Fn(&[f64]) -> Vec<f64>,
    {
        let n = x0.len();
        let mut x = x0.to_vec();
        let eps = 1e-7;

        for iter in 0..self.max_iter {
            let r = residual(&x);
            let r_norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
            if r_norm < self.tolerance {
                return (x, iter);
            }
            // Finite-difference Jacobian J_ij = (R_i(x + eps*e_j) - R_i(x)) / eps
            let mut j = vec![0.0_f64; n * n];
            for j_col in 0..n {
                let mut x_pert = x.clone();
                x_pert[j_col] += eps;
                let r_pert = residual(&x_pert);
                for i in 0..n {
                    j[i * n + j_col] = (r_pert[i] - r[i]) / eps;
                }
            }
            // Solve J * dx = -r using Gaussian elimination
            let dx = self.gaussian_solve(&j, &r.iter().map(|&v| -v).collect::<Vec<_>>(), n);
            for i in 0..n {
                x[i] += self.relaxation * dx[i];
            }
        }
        (x, self.max_iter)
    }

    /// Gaussian elimination for Ax = b.
    fn gaussian_solve(&self, a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut aug: Vec<f64> = (0..n)
            .flat_map(|i| {
                let mut row: Vec<f64> = a[i * n..(i + 1) * n].to_vec();
                row.push(b[i]);
                row
            })
            .collect();

        for col in 0..n {
            // Find pivot
            let mut max_row = col;
            let mut max_val = aug[col * (n + 1) + col].abs();
            for row in col + 1..n {
                let v = aug[row * (n + 1) + col].abs();
                if v > max_val {
                    max_val = v;
                    max_row = row;
                }
            }
            // Swap rows
            for j in 0..=n {
                aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
            let pivot = aug[col * (n + 1) + col];
            if pivot.abs() < 1e-15 {
                continue;
            }
            for row in col + 1..n {
                let factor = aug[row * (n + 1) + col] / pivot;
                for j in col..=n {
                    let v = aug[col * (n + 1) + j];
                    aug[row * (n + 1) + j] -= factor * v;
                }
            }
        }

        // Back substitution
        let mut x = vec![0.0_f64; n];
        for i in (0..n).rev() {
            let mut sum = aug[i * (n + 1) + n];
            for j in i + 1..n {
                sum -= aug[i * (n + 1) + j] * x[j];
            }
            let diag = aug[i * (n + 1) + i];
            x[i] = if diag.abs() > 1e-15 { sum / diag } else { 0.0 };
        }
        x
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thermo_mechanical_thermal_strain() {
        let mat = ThermoMechanical::steel();
        let eps = mat.thermal_strain(393.15); // 120°C
        assert!((eps - 12e-6 * 100.0).abs() < 1e-15);
    }

    #[test]
    fn test_thermo_mechanical_diffusivity() {
        let mat = ThermoMechanical::steel();
        let diff = mat.thermal_diffusivity();
        assert!((diff - 50.0 / (7850.0 * 500.0)).abs() < 1e-15);
    }

    #[test]
    fn test_thermo_mechanical_lame() {
        let mat = ThermoMechanical::steel();
        assert!(mat.lame_lambda() > 0.0);
        assert!(mat.lame_mu() > 0.0);
    }

    #[test]
    fn test_thermo_mechanical_hydrostatic_stress() {
        let mat = ThermoMechanical::steel();
        let sigma = mat.hydrostatic_thermal_stress(393.15);
        assert!(sigma < 0.0); // compression for constrained heating
    }

    #[test]
    fn test_piezo_pzt5h_coupling() {
        let pz = PiezoElectric::pzt5h();
        let k33 = pz.coupling_k33();
        assert!(k33 > 0.0 && k33 < 1.0);
    }

    #[test]
    fn test_piezo_direct_effect() {
        let pz = PiezoElectric::pzt5h();
        let d = pz.direct_effect_d33(1e6); // 1 MPa
        assert!(d.abs() > 0.0);
    }

    #[test]
    fn test_piezo_converse_effect() {
        let pz = PiezoElectric::pzt5h();
        let s = pz.converse_effect_d33(1e6); // 1 MV/m
        assert!(s.abs() > 0.0);
    }

    #[test]
    fn test_piezo_pvdf_name() {
        let pz = PiezoElectric::pvdf();
        assert_eq!(pz.name, "PVDF");
    }

    #[test]
    fn test_magnetostrictive_magnetostriction() {
        let ms = Magnetostrictive::terfenol_d();
        let eps = ms.magnetostriction(ms.m_saturation);
        assert!((eps - 1.5 * ms.lambda_s).abs() < 1e-15);
    }

    #[test]
    fn test_magnetostrictive_curve() {
        let ms = Magnetostrictive::terfenol_d();
        let m = ms.magnetization_curve(0.0);
        assert!(m.abs() < 1e-10); // at H=0, M≈0
    }

    #[test]
    fn test_electrochemical_butler_volmer() {
        let ec = ElectrochemicalMaterial::lithium_graphite();
        let j = ec.butler_volmer(0.0, 298.0);
        assert!(j.abs() < 1e-10); // at zero overpotential, j=0
    }

    #[test]
    fn test_electrochemical_butler_volmer_positive_eta() {
        let ec = ElectrochemicalMaterial::lithium_graphite();
        let j = ec.butler_volmer(0.05, 298.0);
        assert!(j > 0.0); // anodic
    }

    #[test]
    fn test_electrochemical_intercalation_stress() {
        let ec = ElectrochemicalMaterial::lithium_graphite();
        let sigma = ec.intercalation_stress(1000.0, 800.0);
        assert!(sigma < 0.0); // higher than avg → compressive
    }

    #[test]
    fn test_poroelastic_consolidation_coeff() {
        let pe = PoroElastic::saturated_clay();
        let cv = pe.consolidation_coefficient();
        assert!(cv > 0.0);
    }

    #[test]
    fn test_poroelastic_effective_stress() {
        let pe = PoroElastic::saturated_clay();
        let sigma_eff = pe.effective_stress(-100e3, 50e3);
        // -100kPa - 0.9*50kPa = -145kPa
        assert!((sigma_eff - (-100e3 - 0.9 * 50e3)).abs() < 1.0);
    }

    #[test]
    fn test_poroelastic_darcy_velocity() {
        let pe = PoroElastic::saturated_clay();
        let q = pe.darcy_velocity(1000.0); // pressure gradient 1000 Pa/m
        assert!(q < 0.0); // flows in direction of decreasing pressure
    }

    #[test]
    fn test_swelling_hygroscopic_strain() {
        let sw = SwellingMaterial::hydrogel();
        let eps = sw.hygroscopic_strain(0.5); // 50% RH increase
        assert!((eps - 0.002 * 0.5).abs() < 1e-15);
    }

    #[test]
    fn test_swelling_osmotic_pressure() {
        let sw = SwellingMaterial::hydrogel();
        let pi = sw.osmotic_pressure(1.0); // 1 mol/m³
        assert!((pi - 1.0 * 8.314 * 298.15).abs() < 1.0);
    }

    #[test]
    fn test_radiation_damage_dpa() {
        let rd = RadiationDamage::ss316l();
        let dpa = rd.compute_dpa(1e15, 1e-28, 3.15e7); // 1 year
        assert!(dpa > 0.0);
    }

    #[test]
    fn test_radiation_damage_hardening() {
        let mut rd = RadiationDamage::ss316l();
        rd.dpa = 10.0;
        let delta_sy = rd.irradiation_hardening(100e6);
        assert!((delta_sy - 100e6 * 10.0_f64.sqrt()).abs() < 1.0);
    }

    #[test]
    fn test_radiation_damage_dbtt() {
        let mut rd = RadiationDamage::ss316l();
        rd.dpa = 5.0;
        let shift = rd.dbtt_shift_k();
        assert!((shift - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_phase_transform_martensite_fraction() {
        let sma = PhaseTransformMaterial::nitinol();
        let xi = sma.martensite_fraction_cooling(261.0); // below Mf
        assert!((xi - 1.0).abs() < 1e-10);
        let xi2 = sma.martensite_fraction_cooling(310.0); // above Ms
        assert!(xi2.abs() < 1e-10);
    }

    #[test]
    fn test_phase_transform_effective_modulus() {
        let mut sma = PhaseTransformMaterial::nitinol();
        sma.xi = 0.0;
        assert!((sma.effective_modulus() - sma.ea).abs() < 1.0);
        sma.xi = 1.0;
        assert!((sma.effective_modulus() - sma.em).abs() < 1.0);
    }

    #[test]
    fn test_electromagnetic_skin_depth() {
        let em = ElectroMagnetic::silicon_steel();
        let delta = em.skin_depth();
        assert!(delta > 0.0 && delta < 1.0); // mm range for steel at 50 Hz
    }

    #[test]
    fn test_electromagnetic_hysteresis_loss() {
        let em = ElectroMagnetic::silicon_steel();
        let ph = em.hysteresis_loss(1.5); // 1.5 T
        assert!(ph > 0.0);
    }

    #[test]
    fn test_electromagnetic_total_loss_increases_with_b() {
        let em = ElectroMagnetic::silicon_steel();
        let p1 = em.total_core_loss(1.0);
        let p2 = em.total_core_loss(1.5);
        assert!(p2 > p1);
    }

    #[test]
    fn test_coupled_solver_operator_split() {
        let solver = CoupledSolver::new(1e-6, 100, 1.0);
        // Fixed point: a = b + 1, b = a - 1  => converged at any a
        let (a, b, _iter) = solver.operator_split(0.0, 0.0, |b| b + 1.0, |a| a - 1.0);
        assert!((a - b - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_coupled_solver_staggered() {
        let solver = CoupledSolver::new(1e-8, 100, 0.5);
        let (state, _iter) = solver.staggered(&[1.0, 1.0], |s| vec![s[0] * 0.5, s[1] * 0.5]);
        // Should converge toward [0, 0]
        assert!(state[0].abs() < 0.01);
        assert!(state[1].abs() < 0.01);
    }

    #[test]
    fn test_coupled_solver_newton_linear() {
        let solver = CoupledSolver::new(1e-10, 20, 1.0);
        // R(x) = x - 3 => x* = 3
        let (x, _) = solver.monolithic_newton(&[0.0], |s| vec![s[0] - 3.0]);
        assert!((x[0] - 3.0).abs() < 1e-6);
    }
}
