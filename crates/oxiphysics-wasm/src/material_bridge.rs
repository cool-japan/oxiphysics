// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly materials bridge.
//!
//! Provides pure-Rust material models for hyperelastic, plastic, damage,
//! fatigue, viscoelastic, thermal, and equation-of-state computations,
//! ready for serialisation across the WASM boundary.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::{err_to_jsvalue, to_js_value};

// ---------------------------------------------------------------------------
// WasmMaterial
// ---------------------------------------------------------------------------

/// Isotropic linear elastic material with failure properties.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmMaterial {
    /// Human-readable material name.
    #[wasm_bindgen(skip)]
    pub name: String,
    /// Mass density (kg/m³).
    pub density: f64,
    /// Young's modulus (Pa).
    pub elastic_modulus: f64,
    /// Poisson ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Yield strength (Pa).
    pub yield_strength: f64,
    /// Ultimate tensile strength (Pa).
    pub ultimate_strength: f64,
    /// Fracture toughness K_Ic (Pa·√m).
    pub fracture_toughness: f64,
    /// Thermal conductivity (W/m·K).
    pub thermal_conductivity: f64,
    /// Specific heat capacity (J/kg·K).
    pub specific_heat: f64,
    /// Coefficient of linear thermal expansion (1/K).
    pub thermal_expansion: f64,
}

impl Default for WasmMaterial {
    fn default() -> Self {
        WasmMaterial {
            name: "generic".to_string(),
            density: 1000.0,
            elastic_modulus: 70e9,
            poisson_ratio: 0.33,
            yield_strength: 200e6,
            ultimate_strength: 300e6,
            fracture_toughness: 20.0,
            thermal_conductivity: 200.0,
            specific_heat: 900.0,
            thermal_expansion: 23e-6,
        }
    }
}

impl WasmMaterial {
    /// Create a material with the given name and density.
    pub fn new(name: impl Into<String>, density: f64) -> Self {
        WasmMaterial {
            name: name.into(),
            density,
            ..Default::default()
        }
    }

    /// Shear modulus G = E / (2(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.elastic_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Bulk modulus K = E / (3(1-2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.elastic_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Lamé first parameter λ = Eν / ((1+ν)(1-2ν)).
    pub fn lame_lambda(&self) -> f64 {
        let e = self.elastic_modulus;
        let nu = self.poisson_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Wave speed c = sqrt(E/ρ) (longitudinal, uniaxial approximation).
    pub fn wave_speed(&self) -> f64 {
        (self.elastic_modulus / self.density).sqrt()
    }

    /// Thermal diffusivity α = k / (ρ·c_p).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.thermal_conductivity / (self.density * self.specific_heat)
    }

    /// Validate basic ranges.
    pub fn validate(&self) -> Result<(), String> {
        if self.density <= 0.0 {
            return Err("density must be positive".to_string());
        }
        if self.elastic_modulus <= 0.0 {
            return Err("elastic_modulus must be positive".to_string());
        }
        if !(-1.0..0.5).contains(&self.poisson_ratio) {
            return Err("poisson_ratio must be in (-1, 0.5)".to_string());
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// WasmMaterial — wasm-bindgen impl
// ---------------------------------------------------------------------------

#[wasm_bindgen]
impl WasmMaterial {
    /// Construct a material with the given name and density.
    ///
    /// All other properties default to generic aluminium-like values.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(name: String, density: f64) -> WasmMaterial {
        WasmMaterial::new(name, density)
    }

    /// Name getter (String field requires explicit accessor).
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.name.clone()
    }

    /// Name setter.
    #[wasm_bindgen(setter)]
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    /// Shear modulus (Pa).
    #[wasm_bindgen(js_name = "shear_modulus")]
    pub fn shear_modulus_js(&self) -> f64 {
        self.shear_modulus()
    }

    /// Bulk modulus (Pa).
    #[wasm_bindgen(js_name = "bulk_modulus")]
    pub fn bulk_modulus_js(&self) -> f64 {
        self.bulk_modulus()
    }

    /// Lamé first parameter λ (Pa).
    #[wasm_bindgen(js_name = "lame_lambda")]
    pub fn lame_lambda_js(&self) -> f64 {
        self.lame_lambda()
    }

    /// Longitudinal wave speed (m/s).
    #[wasm_bindgen(js_name = "wave_speed")]
    pub fn wave_speed_js(&self) -> f64 {
        self.wave_speed()
    }

    /// Thermal diffusivity (m²/s).
    #[wasm_bindgen(js_name = "thermal_diffusivity")]
    pub fn thermal_diffusivity_js(&self) -> f64 {
        self.thermal_diffusivity()
    }

    /// Validate properties; returns error string on failure, empty string on success.
    #[wasm_bindgen(js_name = "validate")]
    pub fn validate_js(&self) -> String {
        match self.validate() {
            Ok(()) => String::new(),
            Err(e) => e,
        }
    }

    /// Serialise to a JSON string for consumption in JavaScript.
    #[wasm_bindgen(js_name = "to_json")]
    pub fn to_json_js(&self) -> Result<String, JsValue> {
        serde_json::to_string(self).map_err(err_to_jsvalue)
    }

    /// Serialise to a `JsValue` object.
    #[wasm_bindgen(js_name = "to_js_value")]
    pub fn to_js_value_js(&self) -> Result<JsValue, JsValue> {
        to_js_value(self)
    }
}

// ---------------------------------------------------------------------------
// WasmMaterialPreset
// ---------------------------------------------------------------------------

/// Factory for common engineering material presets.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum WasmMaterialPreset {
    /// Structural steel (AISI 1020).
    Steel,
    /// 6061-T6 aluminium alloy.
    Aluminum,
    /// Normal-strength concrete (C30/37).
    Concrete,
    /// Pine wood (along grain).
    Wood,
    /// Natural rubber.
    Rubber,
    /// Soda-lime silica glass.
    Glass,
    /// Grade 5 titanium (Ti-6Al-4V).
    Titanium,
    /// Unidirectional carbon-fibre composite (in-plane average).
    Carbon,
}

impl WasmMaterialPreset {
    /// Return the material data corresponding to the preset.
    pub fn from_preset(preset: WasmMaterialPreset) -> WasmMaterial {
        match preset {
            WasmMaterialPreset::Steel => WasmMaterial {
                name: "steel".to_string(),
                density: 7850.0,
                elastic_modulus: 200e9,
                poisson_ratio: 0.29,
                yield_strength: 250e6,
                ultimate_strength: 400e6,
                fracture_toughness: 50.0,
                thermal_conductivity: 50.0,
                specific_heat: 490.0,
                thermal_expansion: 12e-6,
            },
            WasmMaterialPreset::Aluminum => WasmMaterial {
                name: "aluminum".to_string(),
                density: 2700.0,
                elastic_modulus: 69e9,
                poisson_ratio: 0.33,
                yield_strength: 276e6,
                ultimate_strength: 310e6,
                fracture_toughness: 29.0,
                thermal_conductivity: 167.0,
                specific_heat: 900.0,
                thermal_expansion: 23e-6,
            },
            WasmMaterialPreset::Concrete => WasmMaterial {
                name: "concrete".to_string(),
                density: 2400.0,
                elastic_modulus: 30e9,
                poisson_ratio: 0.2,
                yield_strength: 30e6,
                ultimate_strength: 30e6,
                fracture_toughness: 1.0,
                thermal_conductivity: 1.7,
                specific_heat: 840.0,
                thermal_expansion: 12e-6,
            },
            WasmMaterialPreset::Wood => WasmMaterial {
                name: "wood".to_string(),
                density: 600.0,
                elastic_modulus: 12e9,
                poisson_ratio: 0.35,
                yield_strength: 40e6,
                ultimate_strength: 80e6,
                fracture_toughness: 0.5,
                thermal_conductivity: 0.15,
                specific_heat: 1700.0,
                thermal_expansion: 5e-6,
            },
            WasmMaterialPreset::Rubber => WasmMaterial {
                name: "rubber".to_string(),
                density: 1100.0,
                elastic_modulus: 0.01e9,
                poisson_ratio: 0.499,
                yield_strength: 10e6,
                ultimate_strength: 30e6,
                fracture_toughness: 0.1,
                thermal_conductivity: 0.16,
                specific_heat: 2000.0,
                thermal_expansion: 200e-6,
            },
            WasmMaterialPreset::Glass => WasmMaterial {
                name: "glass".to_string(),
                density: 2500.0,
                elastic_modulus: 70e9,
                poisson_ratio: 0.22,
                yield_strength: 0.0,
                ultimate_strength: 50e6,
                fracture_toughness: 0.7,
                thermal_conductivity: 1.0,
                specific_heat: 840.0,
                thermal_expansion: 9e-6,
            },
            WasmMaterialPreset::Titanium => WasmMaterial {
                name: "titanium".to_string(),
                density: 4430.0,
                elastic_modulus: 114e9,
                poisson_ratio: 0.34,
                yield_strength: 880e6,
                ultimate_strength: 950e6,
                fracture_toughness: 75.0,
                thermal_conductivity: 7.2,
                specific_heat: 560.0,
                thermal_expansion: 8.6e-6,
            },
            WasmMaterialPreset::Carbon => WasmMaterial {
                name: "carbon_cfrp".to_string(),
                density: 1600.0,
                elastic_modulus: 70e9,
                poisson_ratio: 0.1,
                yield_strength: 600e6,
                ultimate_strength: 700e6,
                fracture_toughness: 30.0,
                thermal_conductivity: 5.0,
                specific_heat: 700.0,
                thermal_expansion: 2e-6,
            },
        }
    }

    /// All available presets.
    pub fn all() -> &'static [WasmMaterialPreset] {
        &[
            WasmMaterialPreset::Steel,
            WasmMaterialPreset::Aluminum,
            WasmMaterialPreset::Concrete,
            WasmMaterialPreset::Wood,
            WasmMaterialPreset::Rubber,
            WasmMaterialPreset::Glass,
            WasmMaterialPreset::Titanium,
            WasmMaterialPreset::Carbon,
        ]
    }
}

// ---------------------------------------------------------------------------
// WasmMaterialPreset — wasm-bindgen free functions
// ---------------------------------------------------------------------------

/// Create a `WasmMaterial` from the named preset. Returns a `JsValue` object.
///
/// The `preset` argument is the JS enum value (integer discriminant).
#[wasm_bindgen(js_name = "material_from_preset")]
pub fn material_from_preset_js(preset: WasmMaterialPreset) -> Result<JsValue, JsValue> {
    let mat = WasmMaterialPreset::from_preset(preset);
    to_js_value(&mat)
}

/// Return all preset names as a JSON array of strings.
#[wasm_bindgen(js_name = "material_preset_names")]
pub fn material_preset_names_js() -> String {
    let names: Vec<&str> = vec![
        "Steel", "Aluminum", "Concrete", "Wood", "Rubber", "Glass", "Titanium", "Carbon",
    ];
    serde_json::to_string(&names).unwrap_or_default()
}

// ---------------------------------------------------------------------------
// WasmHyperelastic
// ---------------------------------------------------------------------------

/// Hyperelastic strain energy density model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmHyperelastic {
    /// Neo-Hookean model: W = (μ/2)(I₁-3) - μ ln(J) + (λ/2)(ln J)².
    NeoHookean {
        /// Shear modulus μ (Pa).
        mu: f64,
        /// First Lamé parameter λ (Pa).
        lambda: f64,
    },
    /// Mooney-Rivlin two-parameter model.
    MooneyRivlin {
        /// C10 coefficient (Pa).
        c10: f64,
        /// C01 coefficient (Pa).
        c01: f64,
    },
    /// Ogden model (up to N terms).
    Ogden {
        /// Stretch exponents αᵢ.
        alphas: Vec<f64>,
        /// Moduli μᵢ (Pa).
        mus: Vec<f64>,
    },
}

impl WasmHyperelastic {
    /// Strain energy for a uniaxial stretch λ.
    ///
    /// - **NeoHookean**: W = (μ/2)(I₁−3) − μ ln J + (λ/2)(ln J)²  (incompressible: J=1)
    /// - **MooneyRivlin**: W = C₁₀(I₁−3) + C₀₁(I₂−3)
    /// - **Ogden**: W = Σᵣ (μᵣ/αᵣ)(λ^αᵣ + 2λ̃^αᵣ − 3) where λ̃ = 1/sqrt(λ)
    pub fn strain_energy_uniaxial(&self, stretch: f64) -> f64 {
        match self {
            WasmHyperelastic::NeoHookean { mu, lambda } => {
                let i1 = stretch * stretch + 2.0 / stretch;
                let j: f64 = 1.0; // incompressible
                (mu / 2.0) * (i1 - 3.0) - mu * j.ln() + (lambda / 2.0) * j.ln().powi(2)
            }
            WasmHyperelastic::MooneyRivlin { c10, c01 } => {
                let i1 = stretch * stretch + 2.0 / stretch;
                let i2 = 2.0 * stretch + 1.0 / (stretch * stretch);
                c10 * (i1 - 3.0) + c01 * (i2 - 3.0)
            }
            WasmHyperelastic::Ogden { alphas, mus } => {
                let mut w = 0.0;
                let lam2 = 1.0 / stretch.sqrt();
                for (&alpha, &mu) in alphas.iter().zip(mus.iter()) {
                    w += (mu / alpha) * (stretch.powf(alpha) + 2.0 * lam2.powf(alpha) - 3.0);
                }
                w
            }
        }
    }

    /// Initial shear modulus.
    pub fn initial_shear_modulus(&self) -> f64 {
        match self {
            WasmHyperelastic::NeoHookean { mu, .. } => *mu,
            WasmHyperelastic::MooneyRivlin { c10, c01 } => 2.0 * (c10 + c01),
            WasmHyperelastic::Ogden { mus, .. } => mus.iter().sum::<f64>() / 2.0,
        }
    }
}

// ---------------------------------------------------------------------------
// WasmHyperelastic — wasm-bindgen free functions
// (WasmHyperelastic has payload variants; cannot be #[wasm_bindgen] directly.
//  Expose via serde-JSON free functions instead.)
// ---------------------------------------------------------------------------

/// Create a Neo-Hookean hyperelastic model as a `JsValue`.
#[wasm_bindgen(js_name = "hyperelastic_neo_hookean")]
pub fn hyperelastic_neo_hookean_js(mu: f64, lambda: f64) -> Result<JsValue, JsValue> {
    let m = WasmHyperelastic::NeoHookean { mu, lambda };
    to_js_value(&m)
}

/// Create a Mooney-Rivlin hyperelastic model as a `JsValue`.
#[wasm_bindgen(js_name = "hyperelastic_mooney_rivlin")]
pub fn hyperelastic_mooney_rivlin_js(c10: f64, c01: f64) -> Result<JsValue, JsValue> {
    let m = WasmHyperelastic::MooneyRivlin { c10, c01 };
    to_js_value(&m)
}

/// Strain energy at uniaxial stretch for a JSON-encoded `WasmHyperelastic`.
#[wasm_bindgen(js_name = "hyperelastic_strain_energy")]
pub fn hyperelastic_strain_energy_js(model_json: &str, stretch: f64) -> Result<f64, JsValue> {
    let m: WasmHyperelastic = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.strain_energy_uniaxial(stretch))
}

/// Initial shear modulus for a JSON-encoded `WasmHyperelastic`.
#[wasm_bindgen(js_name = "hyperelastic_shear_modulus")]
pub fn hyperelastic_shear_modulus_js(model_json: &str) -> Result<f64, JsValue> {
    let m: WasmHyperelastic = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.initial_shear_modulus())
}

// ---------------------------------------------------------------------------
// WasmPlasticModel
// ---------------------------------------------------------------------------

/// Plasticity model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmPlasticModel {
    /// J2 (von Mises) plasticity with linear isotropic hardening.
    J2 {
        /// Initial yield stress σ_y (Pa).
        yield_stress: f64,
        /// Linear hardening modulus H (Pa).
        hardening: f64,
    },
    /// Drucker-Prager cone (pressure-dependent yield).
    DruckerPrager {
        /// Cohesion c (Pa).
        cohesion: f64,
        /// Internal friction angle φ (radians).
        friction_angle: f64,
    },
}

impl WasmPlasticModel {
    /// Von Mises yield criterion: returns `true` if `sigma_eq` exceeds the yield surface.
    pub fn is_yielding(&self, sigma_eq: f64, eps_p: f64, pressure: f64) -> bool {
        match self {
            WasmPlasticModel::J2 {
                yield_stress,
                hardening,
            } => sigma_eq > yield_stress + hardening * eps_p,
            WasmPlasticModel::DruckerPrager {
                cohesion,
                friction_angle,
            } => {
                let phi = *friction_angle;
                let k = cohesion * phi.cos() - pressure * phi.sin();
                sigma_eq > k
            }
        }
    }

    /// Current yield stress given accumulated plastic strain `eps_p`.
    pub fn current_yield_stress(&self, eps_p: f64) -> f64 {
        match self {
            WasmPlasticModel::J2 {
                yield_stress,
                hardening,
            } => yield_stress + hardening * eps_p,
            WasmPlasticModel::DruckerPrager { cohesion, .. } => *cohesion,
        }
    }
}

// ---------------------------------------------------------------------------
// WasmPlasticModel — wasm-bindgen free functions
// (WasmPlasticModel has payload variants; cannot be #[wasm_bindgen] directly.)
// ---------------------------------------------------------------------------

/// Create a J2 plasticity model as a `JsValue`.
#[wasm_bindgen(js_name = "plastic_j2")]
pub fn plastic_j2_js(yield_stress: f64, hardening: f64) -> Result<JsValue, JsValue> {
    let m = WasmPlasticModel::J2 {
        yield_stress,
        hardening,
    };
    to_js_value(&m)
}

/// Create a Drucker-Prager plasticity model as a `JsValue`.
#[wasm_bindgen(js_name = "plastic_drucker_prager")]
pub fn plastic_drucker_prager_js(cohesion: f64, friction_angle: f64) -> Result<JsValue, JsValue> {
    let m = WasmPlasticModel::DruckerPrager {
        cohesion,
        friction_angle,
    };
    to_js_value(&m)
}

/// Check yielding for a JSON-encoded `WasmPlasticModel`.
#[wasm_bindgen(js_name = "plastic_is_yielding")]
pub fn plastic_is_yielding_js(
    model_json: &str,
    sigma_eq: f64,
    eps_p: f64,
    pressure: f64,
) -> Result<bool, JsValue> {
    let m: WasmPlasticModel = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.is_yielding(sigma_eq, eps_p, pressure))
}

// ---------------------------------------------------------------------------
// WasmDamageModel
// ---------------------------------------------------------------------------

/// Continuum damage model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmDamageModel {
    /// Brittle damage driven by fracture energy G_f.
    Brittle {
        /// Mode-I fracture energy (J/m²).
        fracture_energy: f64,
        /// Characteristic element length (m).
        element_length: f64,
    },
    /// Ductile damage model (Lemaitre-type).
    Ductile {
        /// Failure strain at fracture.
        eps_f: f64,
        /// Stress triaxiality at fracture.
        triaxiality: f64,
    },
    /// Gurson-Tvergaard-Needleman (GTN) void growth model.
    Gurson {
        /// Initial void fraction f₀.
        f0: f64,
        /// Void nucleation volume fraction f_n.
        fn_: f64,
        /// Nucleation strain standard deviation s_n.
        sn: f64,
    },
}

impl WasmDamageModel {
    /// Compute damage variable D ∈ \[0,1\] for the given equivalent strain.
    pub fn compute_damage(&self, eps_eq: f64) -> f64 {
        match self {
            WasmDamageModel::Brittle {
                fracture_energy,
                element_length,
            } => {
                let eps_0 = fracture_energy / element_length.max(1e-12);
                (1.0 - (-eps_eq / eps_0.max(1e-15)).exp()).clamp(0.0, 1.0)
            }
            WasmDamageModel::Ductile { eps_f, triaxiality } => {
                let threshold = eps_f * (-1.5 * triaxiality).exp();
                (eps_eq / threshold.max(1e-15)).clamp(0.0, 1.0)
            }
            WasmDamageModel::Gurson { f0, fn_, sn } => {
                // Nucleation increment (Gaussian) – simplified scalar version.
                let nucleation = fn_ / (sn * (2.0 * std::f64::consts::PI).sqrt())
                    * (-(eps_eq * eps_eq) / (2.0 * sn * sn)).exp();
                (f0 + nucleation * eps_eq).clamp(0.0, 1.0)
            }
        }
    }

    /// Returns `true` if the element is fully damaged (D ≥ 0.999).
    pub fn is_failed(&self, eps_eq: f64) -> bool {
        self.compute_damage(eps_eq) >= 0.999
    }
}

// ---------------------------------------------------------------------------
// WasmDamageModel — wasm-bindgen free functions
// (WasmDamageModel has payload variants; cannot be #[wasm_bindgen] directly.)
// ---------------------------------------------------------------------------

/// Create a brittle damage model as a `JsValue`.
#[wasm_bindgen(js_name = "damage_brittle")]
pub fn damage_brittle_js(fracture_energy: f64, element_length: f64) -> Result<JsValue, JsValue> {
    let m = WasmDamageModel::Brittle {
        fracture_energy,
        element_length,
    };
    to_js_value(&m)
}

/// Create a ductile damage model as a `JsValue`.
#[wasm_bindgen(js_name = "damage_ductile")]
pub fn damage_ductile_js(eps_f: f64, triaxiality: f64) -> Result<JsValue, JsValue> {
    let m = WasmDamageModel::Ductile { eps_f, triaxiality };
    to_js_value(&m)
}

/// Compute damage for a JSON-encoded `WasmDamageModel`.
#[wasm_bindgen(js_name = "damage_compute")]
pub fn damage_compute_js(model_json: &str, eps_eq: f64) -> Result<f64, JsValue> {
    let m: WasmDamageModel = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.compute_damage(eps_eq))
}

// ---------------------------------------------------------------------------
// WasmFatigue
// ---------------------------------------------------------------------------

/// Stress-life (S-N) Basquin model.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmBasquin {
    /// Fatigue strength coefficient σ'_f (Pa).
    pub sigma_f: f64,
    /// Fatigue strength exponent b (negative).
    pub b: f64,
}

impl WasmBasquin {
    /// Cycles to failure N_f for stress amplitude σ_a.
    pub fn cycles_to_failure(&self, sigma_a: f64) -> f64 {
        if sigma_a <= 0.0 {
            return f64::INFINITY;
        }
        (sigma_a / self.sigma_f).powf(1.0 / self.b)
    }
}

#[wasm_bindgen]
impl WasmBasquin {
    /// Construct a Basquin S-N model.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(sigma_f: f64, b: f64) -> WasmBasquin {
        WasmBasquin { sigma_f, b }
    }

    /// Cycles to failure for the given stress amplitude (Pa).
    #[wasm_bindgen(js_name = "cycles_to_failure")]
    pub fn cycles_to_failure_js(&self, sigma_a: f64) -> f64 {
        self.cycles_to_failure(sigma_a)
    }
}

/// Strain-life (ε-N) Coffin-Manson model.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCoffinManson {
    /// Fatigue ductility coefficient ε'_f.
    pub eps_f: f64,
    /// Fatigue ductility exponent c (negative).
    pub c: f64,
    /// Combined with Basquin for total strain amplitude.
    /// Skipped because nested `#[wasm_bindgen]` structs cannot be pub fields.
    #[wasm_bindgen(skip)]
    pub basquin: WasmBasquin,
    /// Young's modulus for elastic strain (Pa).
    pub elastic_modulus: f64,
}

impl WasmCoffinManson {
    /// Total strain amplitude for given cycles N.
    pub fn strain_amplitude(&self, n_cycles: f64) -> f64 {
        let elastic =
            self.basquin.sigma_f / self.elastic_modulus * (2.0 * n_cycles).powf(self.basquin.b);
        let plastic = self.eps_f * (2.0 * n_cycles).powf(self.c);
        elastic + plastic
    }
}

#[wasm_bindgen]
impl WasmCoffinManson {
    /// Construct a Coffin-Manson model.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new(
        eps_f: f64,
        c: f64,
        sigma_f: f64,
        b: f64,
        elastic_modulus: f64,
    ) -> WasmCoffinManson {
        WasmCoffinManson {
            eps_f,
            c,
            basquin: WasmBasquin { sigma_f, b },
            elastic_modulus,
        }
    }

    /// Total strain amplitude for N cycles.
    #[wasm_bindgen(js_name = "strain_amplitude")]
    pub fn strain_amplitude_js(&self, n_cycles: f64) -> f64 {
        self.strain_amplitude(n_cycles)
    }
}

/// Simple rainflow cycle counter output.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmFatigue {
    /// Accumulated damage D ∈ \[0,1\] via Miner's rule.
    pub miner_damage: f64,
    /// Number of counted cycles.
    pub cycle_count: u64,
    /// History of stress amplitudes.
    /// Skipped: `Vec<f64>` is not directly exposed as a pub field via wasm-bindgen.
    #[wasm_bindgen(skip)]
    pub amplitudes: Vec<f64>,
    /// History of mean stresses.
    #[wasm_bindgen(skip)]
    pub mean_stresses: Vec<f64>,
    /// Basquin model for life prediction.
    /// Skipped: `Option<WasmBasquin>` is not wasm-bindgen-compatible as a pub field.
    #[wasm_bindgen(skip)]
    pub basquin: Option<WasmBasquin>,
}

impl WasmFatigue {
    /// Create a new counter with optional Basquin model.
    pub fn new(basquin: Option<WasmBasquin>) -> Self {
        WasmFatigue {
            basquin,
            ..Default::default()
        }
    }

    /// Record a cycle with given amplitude and mean stress.
    pub fn record_cycle(&mut self, amplitude: f64, mean: f64) {
        self.amplitudes.push(amplitude);
        self.mean_stresses.push(mean);
        self.cycle_count += 1;
        if let Some(ref b) = self.basquin {
            let nf = b.cycles_to_failure(amplitude);
            if nf > 0.0 && nf.is_finite() {
                self.miner_damage += 1.0 / nf;
            }
        }
    }

    /// Returns `true` if Miner damage sum ≥ 1 (failure predicted).
    pub fn is_failed(&self) -> bool {
        self.miner_damage >= 1.0
    }

    /// Reset counter.
    pub fn reset(&mut self) {
        self.miner_damage = 0.0;
        self.cycle_count = 0;
        self.amplitudes.clear();
        self.mean_stresses.clear();
    }
}

#[wasm_bindgen]
impl WasmFatigue {
    /// Construct a fatigue counter without a Basquin model.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmFatigue {
        WasmFatigue::new(None)
    }

    /// Construct a fatigue counter with a Basquin S-N model.
    #[wasm_bindgen(js_name = "with_basquin")]
    pub fn wasm_with_basquin(sigma_f: f64, b: f64) -> WasmFatigue {
        WasmFatigue::new(Some(WasmBasquin { sigma_f, b }))
    }

    /// Record one cycle with the given amplitude and mean stress.
    #[wasm_bindgen(js_name = "record_cycle")]
    pub fn record_cycle_js(&mut self, amplitude: f64, mean: f64) {
        self.record_cycle(amplitude, mean);
    }

    /// Returns `true` if Miner's rule predicts failure (D ≥ 1).
    #[wasm_bindgen(js_name = "is_failed")]
    pub fn is_failed_js(&self) -> bool {
        self.is_failed()
    }

    /// Reset all counters and history.
    #[wasm_bindgen(js_name = "reset")]
    pub fn reset_js(&mut self) {
        self.reset();
    }

    /// Return amplitude history as a flat `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_amplitudes")]
    pub fn get_amplitudes_js(&self) -> Vec<f64> {
        self.amplitudes.clone()
    }

    /// Return mean stress history as a flat `Vec<f64>`.
    #[wasm_bindgen(js_name = "get_mean_stresses")]
    pub fn get_mean_stresses_js(&self) -> Vec<f64> {
        self.mean_stresses.clone()
    }
}

// ---------------------------------------------------------------------------
// WasmViscoelastic
// ---------------------------------------------------------------------------

/// Viscoelastic constitutive model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmViscoelastic {
    /// Maxwell model (spring + dashpot in series).
    Maxwell {
        /// Spring stiffness E (Pa).
        modulus: f64,
        /// Dashpot viscosity η (Pa·s).
        eta: f64,
    },
    /// Kelvin-Voigt model (spring ‖ dashpot).
    Kelvin {
        /// Spring stiffness E (Pa).
        modulus: f64,
        /// Dashpot viscosity η (Pa·s).
        eta: f64,
    },
    /// Generalised Maxwell (Prony series): E(t) = E_inf + Σ Eᵢ exp(-t/τᵢ).
    Prony {
        /// Equilibrium modulus E_∞ (Pa).
        modulus_inf: f64,
        /// Prony series moduli Eᵢ (Pa).
        moduli: Vec<f64>,
        /// Relaxation times τᵢ (s).
        relaxation_times: Vec<f64>,
    },
}

impl WasmViscoelastic {
    /// Relaxation modulus E(t) evaluated at time `t`.
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        match self {
            WasmViscoelastic::Maxwell { modulus, eta } => {
                let tau = eta / modulus.max(1e-30);
                modulus * (-t / tau).exp()
            }
            WasmViscoelastic::Kelvin { modulus, .. } => *modulus, // instantaneous
            WasmViscoelastic::Prony {
                modulus_inf,
                moduli,
                relaxation_times,
            } => {
                let transient: f64 = moduli
                    .iter()
                    .zip(relaxation_times.iter())
                    .map(|(&e, &tau)| e * (-t / tau.max(1e-30)).exp())
                    .sum();
                modulus_inf + transient
            }
        }
    }

    /// Creep compliance J(t) (simplified).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        match self {
            WasmViscoelastic::Maxwell { modulus, eta } => 1.0 / modulus + t / eta.max(1e-30),
            WasmViscoelastic::Kelvin { modulus, eta } => {
                let tau = eta / modulus.max(1e-30);
                (1.0 / modulus) * (1.0 - (-t / tau).exp())
            }
            WasmViscoelastic::Prony { modulus_inf, .. } => 1.0 / modulus_inf.max(1e-30),
        }
    }

    /// Relaxation time (longest).
    pub fn max_relaxation_time(&self) -> f64 {
        match self {
            WasmViscoelastic::Maxwell { modulus, eta } => eta / modulus.max(1e-30),
            WasmViscoelastic::Kelvin { modulus, eta } => eta / modulus.max(1e-30),
            WasmViscoelastic::Prony {
                relaxation_times, ..
            } => relaxation_times.iter().cloned().fold(0.0_f64, f64::max),
        }
    }
}

// ---------------------------------------------------------------------------
// WasmViscoelastic — wasm-bindgen free functions
// (WasmViscoelastic has payload variants with Vec; cannot be #[wasm_bindgen] directly.)
// ---------------------------------------------------------------------------

/// Create a Maxwell viscoelastic model as a `JsValue`.
#[wasm_bindgen(js_name = "visco_maxwell")]
pub fn visco_maxwell_js(modulus: f64, eta: f64) -> Result<JsValue, JsValue> {
    let m = WasmViscoelastic::Maxwell { modulus, eta };
    to_js_value(&m)
}

/// Create a Kelvin-Voigt viscoelastic model as a `JsValue`.
#[wasm_bindgen(js_name = "visco_kelvin")]
pub fn visco_kelvin_js(modulus: f64, eta: f64) -> Result<JsValue, JsValue> {
    let m = WasmViscoelastic::Kelvin { modulus, eta };
    to_js_value(&m)
}

/// Relaxation modulus E(t) for a JSON-encoded `WasmViscoelastic`.
#[wasm_bindgen(js_name = "visco_relaxation_modulus")]
pub fn visco_relaxation_modulus_js(model_json: &str, t: f64) -> Result<f64, JsValue> {
    let m: WasmViscoelastic = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.relaxation_modulus(t))
}

/// Creep compliance J(t) for a JSON-encoded `WasmViscoelastic`.
#[wasm_bindgen(js_name = "visco_creep_compliance")]
pub fn visco_creep_compliance_js(model_json: &str, t: f64) -> Result<f64, JsValue> {
    let m: WasmViscoelastic = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.creep_compliance(t))
}

// ---------------------------------------------------------------------------
// WasmThermalMaterial
// ---------------------------------------------------------------------------

/// Thermal material properties and simple heat-transfer utilities.
#[wasm_bindgen]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmThermalMaterial {
    /// Thermal conductivity k (W/m·K).
    pub conductivity: f64,
    /// Specific heat c_p (J/kg·K).
    pub specific_heat: f64,
    /// Density ρ (kg/m³).
    pub density: f64,
    /// Linear thermal expansion coefficient α (1/K).
    pub thermal_expansion: f64,
    /// Reference temperature T₀ (K) at which stress = 0.
    pub reference_temperature: f64,
    /// Emissivity ε ∈ \[0,1\] for radiation.
    pub emissivity: f64,
}

impl Default for WasmThermalMaterial {
    fn default() -> Self {
        WasmThermalMaterial {
            conductivity: 200.0,
            specific_heat: 900.0,
            density: 2700.0,
            thermal_expansion: 23e-6,
            reference_temperature: 293.15,
            emissivity: 0.1,
        }
    }
}

impl WasmThermalMaterial {
    /// Thermal diffusivity α = k/(ρ c_p).
    pub fn diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }

    /// Steady-state 1-D heat flux q = k·ΔT/L (W/m²).
    pub fn steady_state_flux(&self, delta_t: f64, length: f64) -> f64 {
        self.conductivity * delta_t / length.max(1e-30)
    }

    /// Transient temperature rise at time `t` for a semi-infinite body
    /// with surface heat flux q₀ (approximate analytical solution, °K).
    pub fn transient_temperature_rise(&self, q0: f64, t: f64, depth: f64) -> f64 {
        let alpha = self.diffusivity();
        if alpha <= 0.0 || t <= 0.0 {
            return 0.0;
        }
        // T(x,t) ≈ (2 q₀ / k) sqrt(αt/π) exp(-x²/(4αt))
        let sqrt_at = (alpha * t / std::f64::consts::PI).sqrt();
        let exp_term = (-depth * depth / (4.0 * alpha * t)).exp();
        (2.0 * q0 / self.conductivity) * sqrt_at * exp_term
    }

    /// Stefan-Boltzmann radiation heat flux (W/m²).
    pub fn radiation_flux(&self, surface_temp_k: f64, ambient_temp_k: f64) -> f64 {
        const SIGMA: f64 = 5.670374419e-8;
        self.emissivity * SIGMA * (surface_temp_k.powi(4) - ambient_temp_k.powi(4))
    }

    /// Thermal stress for 1-D constraint σ = -E α ΔT.
    pub fn thermal_stress(&self, elastic_modulus: f64, delta_t: f64) -> f64 {
        -elastic_modulus * self.thermal_expansion * delta_t
    }
}

#[wasm_bindgen]
impl WasmThermalMaterial {
    /// Construct a thermal material with default aluminium-like properties.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmThermalMaterial {
        WasmThermalMaterial::default()
    }

    /// Thermal diffusivity (m²/s).
    #[wasm_bindgen(js_name = "diffusivity")]
    pub fn diffusivity_js(&self) -> f64 {
        self.diffusivity()
    }

    /// Steady-state 1-D heat flux (W/m²).
    #[wasm_bindgen(js_name = "steady_state_flux")]
    pub fn steady_state_flux_js(&self, delta_t: f64, length: f64) -> f64 {
        self.steady_state_flux(delta_t, length)
    }

    /// Transient surface temperature rise (K).
    #[wasm_bindgen(js_name = "transient_temperature_rise")]
    pub fn transient_temperature_rise_js(&self, q0: f64, t: f64, depth: f64) -> f64 {
        self.transient_temperature_rise(q0, t, depth)
    }

    /// Radiation heat flux (W/m²).
    #[wasm_bindgen(js_name = "radiation_flux")]
    pub fn radiation_flux_js(&self, surface_temp_k: f64, ambient_temp_k: f64) -> f64 {
        self.radiation_flux(surface_temp_k, ambient_temp_k)
    }

    /// Thermal stress (Pa) for a constrained 1-D element.
    #[wasm_bindgen(js_name = "thermal_stress")]
    pub fn thermal_stress_js(&self, elastic_modulus: f64, delta_t: f64) -> f64 {
        self.thermal_stress(elastic_modulus, delta_t)
    }
}

// ---------------------------------------------------------------------------
// WasmMaterialCombiner
// ---------------------------------------------------------------------------

/// Composite material averaging schemes.
#[wasm_bindgen]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AveragingScheme {
    /// Voigt (iso-strain) upper bound.
    Voigt,
    /// Reuss (iso-stress) lower bound.
    Reuss,
    /// Hill (arithmetic average of Voigt and Reuss).
    Hill,
}

/// Utility for computing effective properties of composite materials.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmMaterialCombiner;

impl WasmMaterialCombiner {
    /// Create a new combiner.
    pub fn new() -> Self {
        Self
    }

    /// Compute the effective elastic modulus for a two-phase composite.
    ///
    /// `e1`, `e2` are moduli (Pa), `f1` is the volume fraction of phase 1.
    pub fn composite_modulus(&self, e1: f64, e2: f64, f1: f64, scheme: AveragingScheme) -> f64 {
        let f2 = 1.0 - f1;
        match scheme {
            AveragingScheme::Voigt => f1 * e1 + f2 * e2,
            AveragingScheme::Reuss => {
                let denom = f1 / e1.max(1e-30) + f2 / e2.max(1e-30);
                1.0 / denom.max(1e-30)
            }
            AveragingScheme::Hill => {
                let voigt = f1 * e1 + f2 * e2;
                let reuss = 1.0 / (f1 / e1.max(1e-30) + f2 / e2.max(1e-30)).max(1e-30);
                (voigt + reuss) / 2.0
            }
        }
    }

    /// Effective density (always Voigt = linear rule of mixtures).
    pub fn composite_density(&self, rho1: f64, rho2: f64, f1: f64) -> f64 {
        f1 * rho1 + (1.0 - f1) * rho2
    }

    /// Hashin-Shtrikman lower bound on bulk modulus.
    pub fn hs_lower_bulk(&self, k1: f64, k2: f64, g1: f64, f1: f64) -> f64 {
        let f2 = 1.0 - f1;
        k1 + f2 / (1.0 / (k2 - k1).max(1e-30) + 3.0 * f1 / (3.0 * k1 + 4.0 * g1).max(1e-30))
    }
}

#[wasm_bindgen]
impl WasmMaterialCombiner {
    /// Construct a material combiner.
    #[wasm_bindgen(constructor)]
    pub fn wasm_new() -> WasmMaterialCombiner {
        WasmMaterialCombiner::new()
    }

    /// Effective modulus for a two-phase composite using the given averaging scheme.
    #[wasm_bindgen(js_name = "composite_modulus")]
    pub fn composite_modulus_js(&self, e1: f64, e2: f64, f1: f64, scheme: AveragingScheme) -> f64 {
        self.composite_modulus(e1, e2, f1, scheme)
    }

    /// Effective density (rule of mixtures).
    #[wasm_bindgen(js_name = "composite_density")]
    pub fn composite_density_js(&self, rho1: f64, rho2: f64, f1: f64) -> f64 {
        self.composite_density(rho1, rho2, f1)
    }

    /// Hashin-Shtrikman lower bound on bulk modulus.
    #[wasm_bindgen(js_name = "hs_lower_bulk")]
    pub fn hs_lower_bulk_js(&self, k1: f64, k2: f64, g1: f64, f1: f64) -> f64 {
        self.hs_lower_bulk(k1, k2, g1, f1)
    }
}

// ---------------------------------------------------------------------------
// WasmEOS
// ---------------------------------------------------------------------------

/// Equation of State (EOS) model selection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum WasmEOS {
    /// Ideal gas: P = (γ-1)ρe.
    IdealGas {
        /// Adiabatic index γ.
        gamma: f64,
    },
    /// Mie-Grüneisen EOS: P = P_H + Γρ(e - e_H).
    MieGruneisen {
        /// Reference density ρ₀ (kg/m³).
        rho0: f64,
        /// Reference speed of sound c₀ (m/s).
        c0: f64,
        /// Slope of Us-Up Hugoniot s.
        s: f64,
        /// Grüneisen parameter Γ₀.
        gamma0: f64,
    },
    /// Tillotson EOS for high-velocity impact.
    Tillotson {
        /// Reference density ρ₀ (kg/m³).
        rho0: f64,
        /// Parameter a.
        a: f64,
        /// Parameter b.
        b: f64,
        /// Sublimation energy E_s (J/kg).
        e_sub: f64,
        /// Cold compression energy E_0 (J/kg).
        e0: f64,
        /// Parameter A (Pa).
        big_a: f64,
        /// Parameter B (Pa).
        big_b: f64,
        /// Adiabatic index α_t.
        alpha: f64,
        /// Adiabatic index β_t.
        beta: f64,
    },
}

impl WasmEOS {
    /// Compute pressure (Pa) given density `rho` (kg/m³) and specific energy `e` (J/kg).
    pub fn compute_pressure(&self, rho: f64, e: f64) -> f64 {
        match self {
            WasmEOS::IdealGas { gamma } => (gamma - 1.0) * rho * e,
            WasmEOS::MieGruneisen {
                rho0,
                c0,
                s,
                gamma0,
            } => {
                let mu = rho / rho0 - 1.0;
                let denom = (1.0 - (s - 1.0) * mu).powi(2).max(1e-10);
                let p_h = rho0 * c0 * c0 * mu * (1.0 + (1.0 - gamma0 / 2.0) * mu) / denom;
                let e_h = p_h * mu / (2.0 * rho0 * (1.0 + mu)).max(1e-10);
                p_h + gamma0 * rho * (e - e_h)
            }
            WasmEOS::Tillotson {
                rho0,
                a,
                b,
                e0,
                big_a,
                big_b,
                ..
            } => {
                let eta = rho / rho0;
                let mu = eta - 1.0;
                a * rho * e
                    + (b * rho * e / (e / e0 + 1.0).max(1e-30) + big_a * mu + big_b * mu * mu)
                        / (eta * eta).max(1e-30)
            }
        }
    }

    /// Speed of sound c = sqrt(∂P/∂ρ)|_s for each EOS model.
    ///
    /// - **IdealGas**: c = sqrt(γ P / ρ)
    /// - **MieGruneisen**: isentropic linearization c = sqrt(K_s / ρ), K_s = ρ₀·c₀²·(1+2s·μ)
    /// - **Tillotson**: cold-curve bulk modulus K = A + 2B·μ(1−μ), c = sqrt(K / ρ)
    pub fn sound_speed(&self, rho: f64, e: f64) -> f64 {
        match self {
            WasmEOS::IdealGas { gamma } => {
                let p = self.compute_pressure(rho, e);
                (gamma * p / rho.max(1e-30)).max(0.0).sqrt()
            }
            WasmEOS::MieGruneisen { rho0, c0, s, .. } => {
                let mu = rho / rho0 - 1.0;
                let k_s = rho0 * c0 * c0 * (1.0 + 2.0 * s * mu).max(0.0);
                (k_s / rho.max(1e-30)).max(0.0).sqrt()
            }
            WasmEOS::Tillotson {
                rho0, big_a, big_b, ..
            } => {
                let mu = rho / rho0 - 1.0;
                let k = (big_a + 2.0 * big_b * mu * (1.0 - mu)).max(0.0);
                (k / rho.max(1e-30)).max(0.0).sqrt()
            }
        }
    }
}

// ---------------------------------------------------------------------------
// WasmEOS — wasm-bindgen free functions
// (WasmEOS has payload variants; cannot be #[wasm_bindgen] directly.)
// ---------------------------------------------------------------------------

/// Create an ideal-gas EOS as a `JsValue`.
#[wasm_bindgen(js_name = "eos_ideal_gas")]
pub fn eos_ideal_gas_js(gamma: f64) -> Result<JsValue, JsValue> {
    let m = WasmEOS::IdealGas { gamma };
    to_js_value(&m)
}

/// Create a Mie-Grüneisen EOS as a `JsValue`.
#[wasm_bindgen(js_name = "eos_mie_gruneisen")]
pub fn eos_mie_gruneisen_js(rho0: f64, c0: f64, s: f64, gamma0: f64) -> Result<JsValue, JsValue> {
    let m = WasmEOS::MieGruneisen {
        rho0,
        c0,
        s,
        gamma0,
    };
    to_js_value(&m)
}

/// Compute pressure (Pa) for a JSON-encoded `WasmEOS`.
#[wasm_bindgen(js_name = "eos_compute_pressure")]
pub fn eos_compute_pressure_js(model_json: &str, rho: f64, e: f64) -> Result<f64, JsValue> {
    let m: WasmEOS = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.compute_pressure(rho, e))
}

/// Compute sound speed (m/s) for a JSON-encoded `WasmEOS`.
#[wasm_bindgen(js_name = "eos_sound_speed")]
pub fn eos_sound_speed_js(model_json: &str, rho: f64, e: f64) -> Result<f64, JsValue> {
    let m: WasmEOS = serde_json::from_str(model_json).map_err(err_to_jsvalue)?;
    Ok(m.sound_speed(rho, e))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(not(target_arch = "wasm32"))]
#[cfg(test)]
mod tests {
    use super::*;

    // --- WasmMaterial ---

    #[test]
    fn test_material_default_valid() {
        let m = WasmMaterial::default();
        assert!(m.validate().is_ok());
    }

    #[test]
    fn test_material_invalid_density() {
        let m = WasmMaterial {
            density: -1.0,
            ..Default::default()
        };
        assert!(m.validate().is_err());
    }

    #[test]
    fn test_material_shear_modulus() {
        let m = WasmMaterialPreset::from_preset(WasmMaterialPreset::Aluminum);
        let g = m.shear_modulus();
        assert!(g > 0.0);
        // G = E/(2(1+nu)) ≈ 25.9 GPa
        assert!((g - 69e9 / (2.0 * (1.0 + 0.33))).abs() < 1e6);
    }

    #[test]
    fn test_material_bulk_modulus() {
        let m = WasmMaterialPreset::from_preset(WasmMaterialPreset::Steel);
        assert!(m.bulk_modulus() > 0.0);
    }

    #[test]
    fn test_material_wave_speed() {
        let m = WasmMaterialPreset::from_preset(WasmMaterialPreset::Steel);
        let c = m.wave_speed();
        // ~5050 m/s for steel
        assert!(c > 4000.0 && c < 6000.0);
    }

    #[test]
    fn test_material_thermal_diffusivity() {
        let m = WasmMaterialPreset::from_preset(WasmMaterialPreset::Aluminum);
        let alpha = m.thermal_diffusivity();
        assert!(alpha > 0.0);
    }

    #[test]
    fn test_material_serialization() {
        let m = WasmMaterialPreset::from_preset(WasmMaterialPreset::Steel);
        let json = serde_json::to_string(&m).unwrap();
        let m2: WasmMaterial = serde_json::from_str(&json).unwrap();
        assert_eq!(m2.name, "steel");
    }

    // --- WasmMaterialPreset ---

    #[test]
    fn test_all_presets_valid() {
        for &preset in WasmMaterialPreset::all() {
            let m = WasmMaterialPreset::from_preset(preset);
            assert!(m.validate().is_ok(), "preset {:?} invalid", preset);
        }
    }

    #[test]
    fn test_rubber_preset_low_modulus() {
        let m = WasmMaterialPreset::from_preset(WasmMaterialPreset::Rubber);
        assert!(m.elastic_modulus < 1e9);
    }

    // --- WasmHyperelastic ---

    #[test]
    fn test_neohookean_strain_energy_at_one() {
        let nh = WasmHyperelastic::NeoHookean {
            mu: 1e6,
            lambda: 2e6,
        };
        // stretch = 1 should give zero energy
        let w = nh.strain_energy_uniaxial(1.0);
        assert!(
            w.abs() < 1.0,
            "energy at stretch=1 should be near zero, got {w}"
        );
    }

    #[test]
    fn test_neohookean_shear_modulus() {
        let nh = WasmHyperelastic::NeoHookean {
            mu: 1e6,
            lambda: 2e6,
        };
        assert!((nh.initial_shear_modulus() - 1e6).abs() < 1.0);
    }

    #[test]
    fn test_mooney_rivlin_energy_positive() {
        let mr = WasmHyperelastic::MooneyRivlin {
            c10: 0.5e6,
            c01: 0.1e6,
        };
        let w = mr.strain_energy_uniaxial(2.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_ogden_energy_positive() {
        let ogden = WasmHyperelastic::Ogden {
            alphas: vec![2.0, 2.0],
            mus: vec![0.5e6, 0.1e6],
        };
        let w = ogden.strain_energy_uniaxial(2.0);
        assert!(w > 0.0);
    }

    // --- WasmPlasticModel ---

    #[test]
    fn test_j2_not_yielding_below_threshold() {
        let j2 = WasmPlasticModel::J2 {
            yield_stress: 200e6,
            hardening: 10e9,
        };
        assert!(!j2.is_yielding(100e6, 0.0, 0.0));
    }

    #[test]
    fn test_j2_yielding_above_threshold() {
        let j2 = WasmPlasticModel::J2 {
            yield_stress: 200e6,
            hardening: 0.0,
        };
        assert!(j2.is_yielding(300e6, 0.0, 0.0));
    }

    #[test]
    fn test_drucker_prager_yield_stress() {
        let dp = WasmPlasticModel::DruckerPrager {
            cohesion: 5e6,
            friction_angle: 0.5,
        };
        // just check it returns a positive value
        assert!(dp.current_yield_stress(0.0) > 0.0);
    }

    // --- WasmDamageModel ---

    #[test]
    fn test_brittle_damage_zero_at_zero_strain() {
        let dm = WasmDamageModel::Brittle {
            fracture_energy: 100.0,
            element_length: 0.01,
        };
        let d = dm.compute_damage(0.0);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_brittle_damage_approaches_one() {
        let dm = WasmDamageModel::Brittle {
            fracture_energy: 1.0,
            element_length: 0.001,
        };
        // eps_0 = fracture_energy/element_length = 1000, so need eps_eq >> 1000
        let d = dm.compute_damage(1e6);
        assert!(d > 0.99);
    }

    #[test]
    fn test_ductile_damage_clamp() {
        let dm = WasmDamageModel::Ductile {
            eps_f: 0.1,
            triaxiality: 1.0 / 3.0,
        };
        let d = dm.compute_damage(1e6);
        assert!(d <= 1.0);
    }

    // --- WasmFatigue ---

    #[test]
    fn test_basquin_cycles_to_failure_finite() {
        let b = WasmBasquin {
            sigma_f: 900e6,
            b: -0.1,
        };
        let nf = b.cycles_to_failure(300e6);
        assert!(nf.is_finite() && nf > 0.0);
    }

    #[test]
    fn test_basquin_zero_amplitude_is_infinity() {
        let b = WasmBasquin {
            sigma_f: 900e6,
            b: -0.1,
        };
        assert_eq!(b.cycles_to_failure(0.0), f64::INFINITY);
    }

    #[test]
    fn test_fatigue_miner_accumulation() {
        let b = WasmBasquin {
            sigma_f: 1e9,
            b: -0.1,
        };
        let nf = b.cycles_to_failure(500e6);
        let mut fat = WasmFatigue::new(Some(b));
        for _ in 0..10 {
            fat.record_cycle(500e6, 0.0);
        }
        assert!((fat.miner_damage - 10.0 / nf).abs() < 1e-10);
    }

    #[test]
    fn test_fatigue_reset() {
        let mut fat = WasmFatigue::new(None);
        fat.record_cycle(100e6, 0.0);
        fat.reset();
        assert_eq!(fat.cycle_count, 0);
    }

    // --- WasmViscoelastic ---

    #[test]
    fn test_maxwell_relaxation_at_zero() {
        let m = WasmViscoelastic::Maxwell {
            modulus: 1e6,
            eta: 1e4,
        };
        let e0 = m.relaxation_modulus(0.0);
        assert!((e0 - 1e6).abs() < 1.0);
    }

    #[test]
    fn test_maxwell_relaxation_decays() {
        let m = WasmViscoelastic::Maxwell {
            modulus: 1e6,
            eta: 1e4,
        };
        let e0 = m.relaxation_modulus(0.0);
        let e1 = m.relaxation_modulus(1.0);
        assert!(e1 < e0);
    }

    #[test]
    fn test_prony_relaxation() {
        let p = WasmViscoelastic::Prony {
            modulus_inf: 1e6,
            moduli: vec![2e6, 3e6],
            relaxation_times: vec![0.1, 1.0],
        };
        let e_long = p.relaxation_modulus(1e10);
        assert!((e_long - 1e6).abs() < 1e3);
    }

    #[test]
    fn test_kelvin_creep_compliance_nonneg() {
        let k = WasmViscoelastic::Kelvin {
            modulus: 1e6,
            eta: 1e4,
        };
        assert!(k.creep_compliance(1.0) >= 0.0);
    }

    // --- WasmThermalMaterial ---

    #[test]
    fn test_thermal_diffusivity_positive() {
        let t = WasmThermalMaterial::default();
        assert!(t.diffusivity() > 0.0);
    }

    #[test]
    fn test_thermal_steady_flux() {
        let t = WasmThermalMaterial::default();
        let q = t.steady_state_flux(100.0, 0.01);
        assert!((q - 200.0 * 100.0 / 0.01).abs() < 1.0);
    }

    #[test]
    fn test_thermal_radiation_flux_positive() {
        let t = WasmThermalMaterial {
            emissivity: 1.0,
            ..Default::default()
        };
        let q = t.radiation_flux(500.0, 300.0);
        assert!(q > 0.0);
    }

    // --- WasmMaterialCombiner ---

    #[test]
    fn test_voigt_between_components() {
        let c = WasmMaterialCombiner::new();
        let e = c.composite_modulus(100.0, 200.0, 0.5, AveragingScheme::Voigt);
        assert!((e - 150.0).abs() < 1e-10);
    }

    #[test]
    fn test_reuss_between_voigt() {
        let c = WasmMaterialCombiner::new();
        let voigt = c.composite_modulus(100.0, 200.0, 0.5, AveragingScheme::Voigt);
        let reuss = c.composite_modulus(100.0, 200.0, 0.5, AveragingScheme::Reuss);
        assert!(reuss <= voigt);
    }

    #[test]
    fn test_hill_between_voigt_reuss() {
        let c = WasmMaterialCombiner::new();
        let voigt = c.composite_modulus(100.0, 200.0, 0.5, AveragingScheme::Voigt);
        let reuss = c.composite_modulus(100.0, 200.0, 0.5, AveragingScheme::Reuss);
        let hill = c.composite_modulus(100.0, 200.0, 0.5, AveragingScheme::Hill);
        assert!(hill >= reuss && hill <= voigt);
    }

    // --- WasmEOS ---

    #[test]
    fn test_ideal_gas_pressure() {
        let eos = WasmEOS::IdealGas { gamma: 1.4 };
        let p = eos.compute_pressure(1.2, 200_000.0);
        // P = (1.4-1)*1.2*200000 = 96000 Pa
        assert!((p - 96_000.0).abs() < 1.0);
    }

    #[test]
    fn test_ideal_gas_sound_speed() {
        let eos = WasmEOS::IdealGas { gamma: 1.4 };
        let c = eos.sound_speed(1.2, 200_000.0);
        assert!(c > 0.0);
    }

    #[test]
    fn test_mie_gruneisen_compressive() {
        let eos = WasmEOS::MieGruneisen {
            rho0: 2700.0,
            c0: 5300.0,
            s: 1.34,
            gamma0: 2.0,
        };
        let p = eos.compute_pressure(2700.0, 0.0);
        // At reference density, mu=0, P_H = 0
        assert!(p.abs() < 1e3);
    }

    // --- H1: Hugoniot internal energy ---

    #[test]
    fn hugoniot_energy_zero_at_reference() {
        // At rho == rho0, mu = 0, so p_H = 0 and e_H = 0.
        let eos = WasmEOS::MieGruneisen {
            rho0: 8930.0,
            c0: 3940.0,
            s: 1.489,
            gamma0: 1.99,
        };
        let p = eos.compute_pressure(8930.0, 0.0);
        // With e=0 and e_H=0, pressure should also be near zero.
        assert!(p.abs() < 1.0, "pressure at reference should be ~0, got {p}");
    }

    #[test]
    fn hugoniot_energy_positive_compression() {
        // mu = rho/rho0 - 1 = 0.1 > 0 → e_H > 0 and p_H > 0.
        let eos = WasmEOS::MieGruneisen {
            rho0: 8930.0,
            c0: 3940.0,
            s: 1.489,
            gamma0: 1.99,
        };
        let rho_compressed = 8930.0 * 1.1;
        let p = eos.compute_pressure(rho_compressed, 0.0);
        assert!(
            p > 0.0,
            "pressure under compression must be positive, got {p}"
        );
    }

    // --- H2: Sound speed ---

    #[test]
    fn sound_speed_positive_all_materials() {
        // IdealGas: air at roughly 300 K (e = c_v * T ≈ 716 * 300 = 214800 J/kg)
        let air = WasmEOS::IdealGas { gamma: 1.4 };
        let c_air = air.sound_speed(1.2, 214_800.0);
        assert!(c_air > 0.0, "air sound speed must be positive");

        // MieGruneisen: copper compressed 10%
        let copper = WasmEOS::MieGruneisen {
            rho0: 8930.0,
            c0: 3940.0,
            s: 1.489,
            gamma0: 1.99,
        };
        let c_cu = copper.sound_speed(9823.0, 0.0);
        assert!(c_cu > 0.0, "copper sound speed must be positive");

        // Tillotson: aluminium-like parameters
        let al = WasmEOS::Tillotson {
            rho0: 2700.0,
            a: 0.5,
            b: 1.63,
            e_sub: 3.4e8,
            e0: 3.0e8,
            big_a: 7.52e10,
            big_b: 6.5e10,
            alpha: 5.0,
            beta: 5.0,
        };
        let c_al = al.sound_speed(2700.0, 0.0);
        assert!(
            c_al > 0.0,
            "aluminium Tillotson sound speed must be positive"
        );
    }

    #[test]
    fn sound_speed_ideal_gas_air_at_300k() {
        // e = c_v * T; for diatomic ideal gas c_v = R/(M*(γ-1)) where M=0.029 kg/mol
        // c_v ≈ 8.314 / (0.029 * 0.4) ≈ 716 J/(kg·K)  → e ≈ 214800 J/kg
        let air = WasmEOS::IdealGas { gamma: 1.4 };
        let c = air.sound_speed(1.2, 214_800.0);
        // Expected ≈ 347 m/s ± 5 m/s
        assert!(
            (c - 347.0).abs() < 5.0,
            "air sound speed at 300K should be ~347 m/s, got {c}"
        );
    }

    // --- H3: Mooney-Rivlin uniaxial stretch ---

    #[test]
    fn mooney_rivlin_uniaxial_stretch() {
        let c10 = 0.5e6_f64;
        let c01 = 0.1e6_f64;
        let lambda = 2.0_f64;
        let mr = WasmHyperelastic::MooneyRivlin { c10, c01 };
        let w = mr.strain_energy_uniaxial(lambda);

        // I₁ = λ² + 2/λ,  I₂ = 2λ + 1/λ²
        let i1 = lambda * lambda + 2.0 / lambda;
        let i2 = 2.0 * lambda + 1.0 / (lambda * lambda);
        let expected = c10 * (i1 - 3.0) + c01 * (i2 - 3.0);

        assert!(
            (w - expected).abs() < 1.0,
            "Mooney-Rivlin energy mismatch: got {w}, expected {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // Integration tests: wasm-bindgen JS surface (non-wasm32 serde round-trips)
    // -----------------------------------------------------------------------

    /// Verify that WasmMaterial serialises to JSON and the name field survives the
    /// round-trip through the JS-surface getter pattern (name() accessor).
    #[test]
    fn integration_wasm_material_name_roundtrip() {
        let m = WasmMaterial::wasm_new("test_steel".to_string(), 7850.0);
        assert_eq!(m.name(), "test_steel");
        let json = m.to_json_js().expect("serialise to JSON string");
        let parsed: serde_json::Value = serde_json::from_str(&json).expect("parse JSON");
        assert_eq!(parsed["name"].as_str(), Some("test_steel"));
        assert!((parsed["density"].as_f64().expect("density") - 7850.0).abs() < 1e-6);
    }

    /// Verify that the `hyperelastic_strain_energy_js` JSON path produces the
    /// same result as calling the Rust method directly.
    #[test]
    fn integration_hyperelastic_json_path_matches_direct() {
        let mu = 2.0e6_f64;
        let lambda = 4.0e6_f64;
        let stretch = 1.5_f64;
        let nh = WasmHyperelastic::NeoHookean { mu, lambda };
        let expected = nh.strain_energy_uniaxial(stretch);

        let json = serde_json::to_string(&nh).expect("serialise");
        let from_json = hyperelastic_strain_energy_js(&json, stretch)
            .expect("hyperelastic_strain_energy_js should succeed");
        assert!(
            (from_json - expected).abs() < 1.0,
            "JSON path strain energy {from_json} != direct {expected}"
        );
    }
}
