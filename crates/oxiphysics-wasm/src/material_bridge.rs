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

// ---------------------------------------------------------------------------
// WasmMaterial
// ---------------------------------------------------------------------------

/// Isotropic linear elastic material with failure properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmMaterial {
    /// Human-readable material name.
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

    /// Wave speed c = √(E/ρ) (longitudinal, uniaxial approximation).
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
// WasmMaterialPreset
// ---------------------------------------------------------------------------

/// Factory for common engineering material presets.
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
    /// Strain energy for a uniaxial stretch λ (NeoHookean only, others return 0 as stub).
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
// WasmFatigue
// ---------------------------------------------------------------------------

/// Stress-life (S-N) Basquin model.
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

/// Strain-life (ε-N) Coffin-Manson model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WasmCoffinManson {
    /// Fatigue ductility coefficient ε'_f.
    pub eps_f: f64,
    /// Fatigue ductility exponent c (negative).
    pub c: f64,
    /// Combined with Basquin for total strain amplitude.
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

/// Simple rainflow cycle counter output.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct WasmFatigue {
    /// Accumulated damage D ∈ \[0,1\] via Miner's rule.
    pub miner_damage: f64,
    /// Number of counted cycles.
    pub cycle_count: u64,
    /// History of stress amplitudes.
    pub amplitudes: Vec<f64>,
    /// History of mean stresses.
    pub mean_stresses: Vec<f64>,
    /// Basquin model for life prediction.
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
// WasmThermalMaterial
// ---------------------------------------------------------------------------

/// Thermal material properties and simple heat-transfer utilities.
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
        // T(x,t) ≈ (2 q₀ / k) √(αt/π) exp(-x²/(4αt))
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

// ---------------------------------------------------------------------------
// WasmMaterialCombiner
// ---------------------------------------------------------------------------

/// Composite material averaging schemes.
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
                let p_h = rho0 * c0 * c0 * mu * (1.0 + (1.0 - gamma0 / 2.0) * mu)
                    / (1.0 - (s - 1.0) * mu).powi(2).max(1e-30);
                let e_h = 0.0; // reference Hugoniot energy (stub)
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

    /// Speed of sound c = √(∂P/∂ρ) for ideal gas (others return 0 as stub).
    pub fn sound_speed(&self, rho: f64, e: f64) -> f64 {
        match self {
            WasmEOS::IdealGas { gamma } => {
                let p = self.compute_pressure(rho, e);
                (gamma * p / rho.max(1e-30)).sqrt()
            }
            _ => 0.0,
        }
    }
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
}
