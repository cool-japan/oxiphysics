// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Materials API for Python interop.
//!
//! Provides a comprehensive set of material models for structural, thermal,
//! acoustic, and fatigue analysis. All types use plain `f64` and `Vec`f64`
//! — no nalgebra — for easy FFI transmission.

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// PyElasticMaterial
// ---------------------------------------------------------------------------

/// Linear elastic isotropic material.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyElasticMaterial {
    /// Young's modulus E (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio ν (dimensionless, 0..0.5).
    pub poisson_ratio: f64,
    /// Mass density ρ (kg/m³).
    pub density: f64,
}

impl PyElasticMaterial {
    /// Create a new linear elastic material.
    pub fn new(youngs_modulus: f64, poisson_ratio: f64, density: f64) -> Self {
        Self {
            youngs_modulus,
            poisson_ratio,
            density,
        }
    }

    /// Compute the Cauchy stress vector `\[σxx, σyy, τxy\]` from the
    /// engineering strain vector `\[εxx, εyy, γxy\]` using plane-stress
    /// constitutive relations.
    pub fn stress_from_strain(&self, strain: [f64; 3]) -> [f64; 3] {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        let c = e / (1.0 - nu * nu);
        [
            c * (strain[0] + nu * strain[1]),
            c * (nu * strain[0] + strain[1]),
            c * (1.0 - nu) * 0.5 * strain[2],
        ]
    }

    /// Compute the elastic strain energy density (J/m³) given a strain vector.
    pub fn strain_energy_density(&self, strain: [f64; 3]) -> f64 {
        let stress = self.stress_from_strain(strain);
        0.5 * (stress[0] * strain[0] + stress[1] * strain[1] + stress[2] * strain[2])
    }

    /// Shear modulus G derived from E and ν.
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Bulk modulus K derived from E and ν.
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
}

impl Default for PyElasticMaterial {
    fn default() -> Self {
        // Steel-like defaults
        Self::new(200e9, 0.3, 7850.0)
    }
}

// ---------------------------------------------------------------------------
// PyHyperelasticMaterial
// ---------------------------------------------------------------------------

/// Hyperelastic material supporting NeoHookean and Mooney–Rivlin models.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyHyperelasticMaterial {
    /// First Mooney–Rivlin constant C₁ (also the NeoHookean shear parameter).
    pub c1: f64,
    /// Second Mooney–Rivlin constant C₂.
    pub c2: f64,
    /// Bulk modulus κ (Pa).
    pub bulk_modulus: f64,
}

impl PyHyperelasticMaterial {
    /// Create a new hyperelastic material.
    pub fn new(c1: f64, c2: f64, bulk_modulus: f64) -> Self {
        Self {
            c1,
            c2,
            bulk_modulus,
        }
    }

    /// NeoHookean strain energy density W = C₁(Ī₁ − 3) + κ/2 (J − 1)².
    ///
    /// `i1_bar` is the first invariant of the isochoric right Cauchy–Green
    /// tensor; `j` is the volumetric Jacobian.
    pub fn strain_energy_neo_hookean(&self, i1_bar: f64, j: f64) -> f64 {
        self.c1 * (i1_bar - 3.0) + 0.5 * self.bulk_modulus * (j - 1.0).powi(2)
    }

    /// Mooney–Rivlin strain energy density W = C₁(Ī₁ − 3) + C₂(Ī₂ − 3) + κ/2 (J − 1)².
    ///
    /// `i1_bar` and `i2_bar` are the first and second invariants of the
    /// isochoric right Cauchy–Green tensor; `j` is the volumetric Jacobian.
    pub fn strain_energy_mooney_rivlin(&self, i1_bar: f64, i2_bar: f64, j: f64) -> f64 {
        self.c1 * (i1_bar - 3.0)
            + self.c2 * (i2_bar - 3.0)
            + 0.5 * self.bulk_modulus * (j - 1.0).powi(2)
    }

    /// Initial (small-strain) shear modulus μ = 2(C₁ + C₂).
    pub fn initial_shear_modulus(&self) -> f64 {
        2.0 * (self.c1 + self.c2)
    }
}

impl Default for PyHyperelasticMaterial {
    fn default() -> Self {
        // Soft rubber-like defaults
        Self::new(0.4e6, 0.1e6, 2.0e9)
    }
}

// ---------------------------------------------------------------------------
// PyPlasticMaterial
// ---------------------------------------------------------------------------

/// Elastoplastic material with linear isotropic hardening (von Mises).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPlasticMaterial {
    /// Initial yield stress σ_y (Pa).
    pub yield_stress: f64,
    /// Isotropic hardening modulus H (Pa).
    pub hardening_modulus: f64,
    /// Young's modulus E (Pa) — needed for elastic predictor.
    pub youngs_modulus: f64,
}

impl PyPlasticMaterial {
    /// Create a new plastic material.
    pub fn new(yield_stress: f64, hardening_modulus: f64, youngs_modulus: f64) -> Self {
        Self {
            yield_stress,
            hardening_modulus,
            youngs_modulus,
        }
    }

    /// Von Mises yield function f(σ_eq, α) = σ_eq − (σ_y + H α).
    ///
    /// Returns a positive value when yielding; zero or negative when elastic.
    pub fn yield_function(&self, equivalent_stress: f64, accumulated_plastic_strain: f64) -> f64 {
        let hardened_yield =
            self.yield_stress + self.hardening_modulus * accumulated_plastic_strain;
        equivalent_stress - hardened_yield
    }

    /// Current yield stress given the accumulated plastic strain.
    pub fn current_yield_stress(&self, accumulated_plastic_strain: f64) -> f64 {
        self.yield_stress + self.hardening_modulus * accumulated_plastic_strain
    }
}

impl Default for PyPlasticMaterial {
    fn default() -> Self {
        Self::new(250e6, 10e9, 200e9)
    }
}

// ---------------------------------------------------------------------------
// PyFatigueMaterial
// ---------------------------------------------------------------------------

/// Fatigue material described by the Basquin (S–N) power law.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFatigueMaterial {
    /// Basquin coefficient A (stress intercept at N = 1).
    pub basquin_coefficient: f64,
    /// Basquin exponent b (negative slope of the S–N curve).
    pub basquin_exponent: f64,
    /// Endurance limit σ_e (Pa).  Below this stress fatigue life is infinite.
    pub endurance_limit: f64,
}

impl PyFatigueMaterial {
    /// Create a new fatigue material.
    pub fn new(basquin_coefficient: f64, basquin_exponent: f64, endurance_limit: f64) -> Self {
        Self {
            basquin_coefficient,
            basquin_exponent,
            endurance_limit,
        }
    }

    /// Estimate cycles to failure N_f for a given stress amplitude σ_a (Pa).
    ///
    /// Returns `f64::INFINITY` when σ_a ≤ endurance_limit.
    /// Uses the Basquin relation: N_f = (A / σ_a)^(1/b).
    pub fn cycles_to_failure(&self, stress_amplitude: f64) -> f64 {
        if stress_amplitude <= self.endurance_limit {
            return f64::INFINITY;
        }
        // N_f = (A / σ_a)^(1/b)
        (self.basquin_coefficient / stress_amplitude).powf(1.0 / self.basquin_exponent)
    }

    /// Miner's rule damage increment for one block of `n_cycles` at `stress_amplitude`.
    pub fn miner_damage(&self, n_cycles: f64, stress_amplitude: f64) -> f64 {
        let nf = self.cycles_to_failure(stress_amplitude);
        if nf.is_infinite() { 0.0 } else { n_cycles / nf }
    }
}

impl Default for PyFatigueMaterial {
    fn default() -> Self {
        // Representative steel values
        Self::new(900e6, 0.1, 200e6)
    }
}

// ---------------------------------------------------------------------------
// PyThermalMaterial
// ---------------------------------------------------------------------------

/// Isotropic thermal material.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyThermalMaterial {
    /// Thermal conductivity λ (W/m·K).
    pub conductivity: f64,
    /// Specific heat capacity c_p (J/kg·K).
    pub specific_heat: f64,
    /// Mass density ρ (kg/m³).
    pub density: f64,
}

impl PyThermalMaterial {
    /// Create a new thermal material.
    pub fn new(conductivity: f64, specific_heat: f64, density: f64) -> Self {
        Self {
            conductivity,
            specific_heat,
            density,
        }
    }

    /// Thermal diffusivity α = λ / (ρ c_p)  (m²/s).
    pub fn thermal_diffusivity(&self) -> f64 {
        self.conductivity / (self.density * self.specific_heat)
    }

    /// Heat flux magnitude q = −λ ‖∇T‖  (W/m²) for a given temperature gradient magnitude.
    pub fn heat_flux(&self, grad_t_magnitude: f64) -> f64 {
        self.conductivity * grad_t_magnitude
    }
}

impl Default for PyThermalMaterial {
    fn default() -> Self {
        // Steel-like thermal properties
        Self::new(50.0, 480.0, 7850.0)
    }
}

// ---------------------------------------------------------------------------
// PyViscoelasticMaterial
// ---------------------------------------------------------------------------

/// Viscoelastic material described by a generalised Maxwell (Prony series) model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyViscoelasticMaterial {
    /// Long-term (equilibrium) modulus E∞ (Pa).
    pub equilibrium_modulus: f64,
    /// Prony series amplitudes (modulus contributions E_i).
    pub prony_moduli: Vec<f64>,
    /// Prony series relaxation times τ_i (s).
    pub prony_times: Vec<f64>,
}

impl PyViscoelasticMaterial {
    /// Create a new viscoelastic material.
    ///
    /// `prony_moduli` and `prony_times` must have equal length.
    pub fn new(equilibrium_modulus: f64, prony_moduli: Vec<f64>, prony_times: Vec<f64>) -> Self {
        Self {
            equilibrium_modulus,
            prony_moduli,
            prony_times,
        }
    }

    /// Relaxation modulus E(t) = E∞ + Σ E_i exp(−t/τ_i).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let sum: f64 = self
            .prony_moduli
            .iter()
            .zip(self.prony_times.iter())
            .map(|(&ei, &ti)| ei * (-t / ti).exp())
            .sum();
        self.equilibrium_modulus + sum
    }

    /// Instantaneous (glassy) modulus E(0) = E∞ + Σ E_i.
    pub fn glassy_modulus(&self) -> f64 {
        self.equilibrium_modulus + self.prony_moduli.iter().sum::<f64>()
    }
}

impl Default for PyViscoelasticMaterial {
    fn default() -> Self {
        Self::new(1e6, vec![2e6, 1e6], vec![0.1, 1.0])
    }
}

// ---------------------------------------------------------------------------
// PyDamageMaterial
// ---------------------------------------------------------------------------

/// Continuum damage mechanics material (isotropic scalar damage).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyDamageMaterial {
    /// Critical strain energy release rate G_c (J/m³).
    pub critical_energy_release: f64,
    /// Softening slope (negative, in Pa per unit damage).
    pub softening_slope: f64,
    /// Undamaged Young's modulus E₀ (Pa).
    pub undamaged_modulus: f64,
}

impl PyDamageMaterial {
    /// Create a new damage material.
    pub fn new(critical_energy_release: f64, softening_slope: f64, undamaged_modulus: f64) -> Self {
        Self {
            critical_energy_release,
            softening_slope,
            undamaged_modulus,
        }
    }

    /// Scalar damage variable d ∈ [0, 1] for a given strain energy density.
    ///
    /// d = 0 → undamaged; d = 1 → fully failed.
    pub fn damage_variable(&self, strain_energy: f64) -> f64 {
        if strain_energy <= 0.0 {
            return 0.0;
        }
        let d = (strain_energy - self.critical_energy_release)
            / (self.softening_slope.abs() + strain_energy);
        d.clamp(0.0, 1.0)
    }

    /// Effective (damaged) stiffness E_eff = (1 − d) E₀.
    pub fn effective_modulus(&self, strain_energy: f64) -> f64 {
        let d = self.damage_variable(strain_energy);
        (1.0 - d) * self.undamaged_modulus
    }
}

impl Default for PyDamageMaterial {
    fn default() -> Self {
        Self::new(1000.0, 1e8, 200e9)
    }
}

// ---------------------------------------------------------------------------
// PyCreepMaterial
// ---------------------------------------------------------------------------

/// Creep material described by Norton's power-law creep model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCreepMaterial {
    /// Norton creep coefficient A (1/s · Pa^{-n}).
    pub norton_coefficient: f64,
    /// Norton creep exponent n (dimensionless).
    pub norton_exponent: f64,
    /// Activation energy Q divided by universal gas constant R (K).
    pub activation_temperature: f64,
}

impl PyCreepMaterial {
    /// Create a new creep material.
    pub fn new(norton_coefficient: f64, norton_exponent: f64, activation_temperature: f64) -> Self {
        Self {
            norton_coefficient,
            norton_exponent,
            activation_temperature,
        }
    }

    /// Secondary creep rate dε/dt = A σ^n exp(−Q/RT).
    ///
    /// `stress` in Pa, `temp` in K.
    pub fn creep_rate(&self, stress: f64, temp: f64) -> f64 {
        let arrhenius = (-self.activation_temperature / temp).exp();
        self.norton_coefficient * stress.powf(self.norton_exponent) * arrhenius
    }

    /// Minimum creep rate (at reference temperature 293 K).
    pub fn reference_creep_rate(&self, stress: f64) -> f64 {
        self.creep_rate(stress, 293.0)
    }
}

impl Default for PyCreepMaterial {
    fn default() -> Self {
        // Representative high-temperature steel
        Self::new(1e-26, 5.0, 15_000.0)
    }
}

// ---------------------------------------------------------------------------
// PyAcousticMaterial
// ---------------------------------------------------------------------------

/// Acoustic (fluid) material.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyAcousticMaterial {
    /// Bulk modulus κ (Pa).
    pub bulk_modulus: f64,
    /// Mass density ρ (kg/m³).
    pub density: f64,
}

impl PyAcousticMaterial {
    /// Create a new acoustic material.
    pub fn new(bulk_modulus: f64, density: f64) -> Self {
        Self {
            bulk_modulus,
            density,
        }
    }

    /// Speed of sound c = sqrt(κ / ρ)  (m/s).
    pub fn speed_of_sound(&self) -> f64 {
        (self.bulk_modulus / self.density).sqrt()
    }

    /// Acoustic impedance Z = ρ c  (Pa·s/m = rayl).
    pub fn impedance(&self) -> f64 {
        self.density * self.speed_of_sound()
    }

    /// Reflection coefficient at an interface with another acoustic material.
    pub fn reflection_coefficient(&self, other: &PyAcousticMaterial) -> f64 {
        let z1 = self.impedance();
        let z2 = other.impedance();
        (z2 - z1) / (z2 + z1)
    }
}

impl Default for PyAcousticMaterial {
    fn default() -> Self {
        // Water at 20 °C
        Self::new(2.2e9, 998.0)
    }
}

// ---------------------------------------------------------------------------
// PyCompositeMaterial
// ---------------------------------------------------------------------------

/// Composite material with Voigt, Reuss, and Hill mixture rules.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCompositeMaterial {
    /// Modulus of phase 1 (matrix) E₁ (Pa).
    pub modulus1: f64,
    /// Modulus of phase 2 (fiber/particle) E₂ (Pa).
    pub modulus2: f64,
    /// Volume fraction of phase 1 (0..1).
    pub volume_fraction1: f64,
}

impl PyCompositeMaterial {
    /// Create a new composite material.
    pub fn new(modulus1: f64, modulus2: f64, volume_fraction1: f64) -> Self {
        Self {
            modulus1,
            modulus2,
            volume_fraction1,
        }
    }

    /// Voigt (isostrain) upper-bound modulus.
    pub fn voigt_modulus(&self, vf1: f64, e1: f64, e2: f64) -> f64 {
        vf1 * e1 + (1.0 - vf1) * e2
    }

    /// Reuss (isostress) lower-bound modulus.
    pub fn reuss_modulus(&self, vf1: f64, e1: f64, e2: f64) -> f64 {
        let vf2 = 1.0 - vf1;
        if vf1 / e1 + vf2 / e2 == 0.0 {
            0.0
        } else {
            1.0 / (vf1 / e1 + vf2 / e2)
        }
    }

    /// Hill (arithmetic mean of Voigt and Reuss) modulus.
    pub fn hill_modulus(&self) -> f64 {
        let v = self.voigt_modulus(self.volume_fraction1, self.modulus1, self.modulus2);
        let r = self.reuss_modulus(self.volume_fraction1, self.modulus1, self.modulus2);
        0.5 * (v + r)
    }

    /// Hashin–Shtrikman lower bound modulus (spherical inclusions).
    pub fn hashin_shtrikman_lower(&self) -> f64 {
        let vf2 = 1.0 - self.volume_fraction1;
        let e1 = self.modulus1.min(self.modulus2);
        let e2 = self.modulus1.max(self.modulus2);
        e1 + vf2 / (1.0 / (e2 - e1) + self.volume_fraction1 / (3.0 * e1))
    }
}

impl Default for PyCompositeMaterial {
    fn default() -> Self {
        Self::new(70e9, 200e9, 0.6)
    }
}

// ---------------------------------------------------------------------------
// Registration helper
// ---------------------------------------------------------------------------

/// Register all material classes into a Python sub-module named `"materials"`.
///
/// This is a no-op placeholder that documents the intended PyO3 registration
/// point. When PyO3 is enabled as a dependency the body should call
/// `m.add_class::`PyElasticMaterial`()` etc.
pub fn register_materials_module(_m: &str) {
    // Placeholder: actual PyO3 registration would happen here.
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- PyElasticMaterial ---

    #[test]
    fn test_elastic_new() {
        let mat = PyElasticMaterial::new(200e9, 0.3, 7850.0);
        assert_eq!(mat.youngs_modulus, 200e9);
        assert_eq!(mat.poisson_ratio, 0.3);
        assert_eq!(mat.density, 7850.0);
    }

    #[test]
    fn test_elastic_default() {
        let mat = PyElasticMaterial::default();
        assert!(mat.youngs_modulus > 0.0);
    }

    #[test]
    fn test_elastic_shear_modulus() {
        let mat = PyElasticMaterial::new(200e9, 0.3, 7850.0);
        let g = mat.shear_modulus();
        // G = E / (2(1+nu)) = 200e9 / 2.6 ≈ 76.9 GPa
        assert!((g - 76.923e9).abs() < 1e8);
    }

    #[test]
    fn test_elastic_bulk_modulus() {
        let mat = PyElasticMaterial::new(200e9, 0.3, 7850.0);
        let k = mat.bulk_modulus();
        // K = E / (3(1-2nu)) = 200e9 / 1.2 ≈ 166.7 GPa
        assert!((k - 166.667e9).abs() < 1e9);
    }

    #[test]
    fn test_elastic_stress_from_strain_uniaxial() {
        let mat = PyElasticMaterial::new(200e9, 0.3, 7850.0);
        let strain = [1e-3, 0.0, 0.0];
        let stress = mat.stress_from_strain(strain);
        // σxx = E/(1-nu²) * εxx ≈ 219.78 MPa
        assert!(stress[0] > 0.0);
        assert!(stress[2].abs() < 1.0);
    }

    #[test]
    fn test_elastic_strain_energy_positive() {
        let mat = PyElasticMaterial::new(200e9, 0.3, 7850.0);
        let e = mat.strain_energy_density([1e-3, 0.0, 0.0]);
        assert!(e > 0.0);
    }

    #[test]
    fn test_elastic_zero_strain_zero_stress() {
        let mat = PyElasticMaterial::default();
        let s = mat.stress_from_strain([0.0, 0.0, 0.0]);
        assert_eq!(s, [0.0, 0.0, 0.0]);
    }

    // --- PyHyperelasticMaterial ---

    #[test]
    fn test_hyperelastic_new() {
        let mat = PyHyperelasticMaterial::new(0.4e6, 0.1e6, 2.0e9);
        assert_eq!(mat.c1, 0.4e6);
        assert_eq!(mat.c2, 0.1e6);
        assert_eq!(mat.bulk_modulus, 2.0e9);
    }

    #[test]
    fn test_hyperelastic_neo_hookean_zero() {
        let mat = PyHyperelasticMaterial::default();
        // At identity deformation: Ī₁ = 3, J = 1 → W = 0
        let w = mat.strain_energy_neo_hookean(3.0, 1.0);
        assert!(w.abs() < 1e-10);
    }

    #[test]
    fn test_hyperelastic_mooney_rivlin_zero() {
        let mat = PyHyperelasticMaterial::default();
        let w = mat.strain_energy_mooney_rivlin(3.0, 3.0, 1.0);
        assert!(w.abs() < 1e-10);
    }

    #[test]
    fn test_hyperelastic_neo_hookean_positive_stretch() {
        let mat = PyHyperelasticMaterial::new(0.4e6, 0.0, 2.0e9);
        let w = mat.strain_energy_neo_hookean(4.0, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_hyperelastic_initial_shear_modulus() {
        let mat = PyHyperelasticMaterial::new(0.4e6, 0.1e6, 2.0e9);
        let mu = mat.initial_shear_modulus();
        assert!((mu - 1.0e6).abs() < 1.0);
    }

    #[test]
    fn test_hyperelastic_default() {
        let mat = PyHyperelasticMaterial::default();
        assert!(mat.c1 > 0.0);
        assert!(mat.bulk_modulus > 0.0);
    }

    // --- PyPlasticMaterial ---

    #[test]
    fn test_plastic_new() {
        let mat = PyPlasticMaterial::new(250e6, 10e9, 200e9);
        assert_eq!(mat.yield_stress, 250e6);
    }

    #[test]
    fn test_plastic_yield_function_elastic() {
        let mat = PyPlasticMaterial::new(250e6, 10e9, 200e9);
        let f = mat.yield_function(200e6, 0.0);
        assert!(f < 0.0);
    }

    #[test]
    fn test_plastic_yield_function_yielding() {
        let mat = PyPlasticMaterial::new(250e6, 10e9, 200e9);
        let f = mat.yield_function(300e6, 0.0);
        assert!(f > 0.0);
    }

    #[test]
    fn test_plastic_hardening_increases_yield() {
        let mat = PyPlasticMaterial::new(250e6, 10e9, 200e9);
        let y0 = mat.current_yield_stress(0.0);
        let y1 = mat.current_yield_stress(0.01);
        assert!(y1 > y0);
    }

    #[test]
    fn test_plastic_default() {
        let mat = PyPlasticMaterial::default();
        assert!(mat.yield_stress > 0.0);
    }

    // --- PyFatigueMaterial ---

    #[test]
    fn test_fatigue_new() {
        let mat = PyFatigueMaterial::new(900e6, 0.1, 200e6);
        assert_eq!(mat.basquin_coefficient, 900e6);
    }

    #[test]
    fn test_fatigue_below_endurance_infinite_life() {
        let mat = PyFatigueMaterial::default();
        let n = mat.cycles_to_failure(100e6); // below 200 MPa endurance limit
        assert!(n.is_infinite());
    }

    #[test]
    fn test_fatigue_above_endurance_finite_life() {
        let mat = PyFatigueMaterial::default();
        let n = mat.cycles_to_failure(500e6);
        assert!(n.is_finite() && n > 0.0);
    }

    #[test]
    fn test_fatigue_miner_damage_zero_below_endurance() {
        let mat = PyFatigueMaterial::default();
        let d = mat.miner_damage(1000.0, 100e6);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_fatigue_miner_damage_positive_above_endurance() {
        let mat = PyFatigueMaterial::default();
        let d = mat.miner_damage(1000.0, 500e6);
        assert!(d > 0.0);
    }

    #[test]
    fn test_fatigue_default() {
        let mat = PyFatigueMaterial::default();
        assert!(mat.endurance_limit > 0.0);
    }

    // --- PyThermalMaterial ---

    #[test]
    fn test_thermal_new() {
        let mat = PyThermalMaterial::new(50.0, 480.0, 7850.0);
        assert_eq!(mat.conductivity, 50.0);
    }

    #[test]
    fn test_thermal_diffusivity() {
        let mat = PyThermalMaterial::new(50.0, 480.0, 7850.0);
        let alpha = mat.thermal_diffusivity();
        assert!(alpha > 0.0);
        // α = 50 / (7850 * 480) ≈ 1.328e-5 m²/s
        assert!((alpha - 1.328e-5).abs() < 1e-7);
    }

    #[test]
    fn test_thermal_heat_flux() {
        let mat = PyThermalMaterial::default();
        let q = mat.heat_flux(10.0); // 10 K/m gradient
        assert!((q - 500.0).abs() < 1e-6);
    }

    #[test]
    fn test_thermal_default() {
        let mat = PyThermalMaterial::default();
        assert!(mat.specific_heat > 0.0);
    }

    // --- PyViscoelasticMaterial ---

    #[test]
    fn test_viscoelastic_new() {
        let mat = PyViscoelasticMaterial::new(1e6, vec![2e6], vec![0.1]);
        assert_eq!(mat.equilibrium_modulus, 1e6);
    }

    #[test]
    fn test_viscoelastic_relaxation_t0() {
        let mat = PyViscoelasticMaterial::default();
        let e0 = mat.relaxation_modulus(0.0);
        let eg = mat.glassy_modulus();
        assert!((e0 - eg).abs() < 1.0);
    }

    #[test]
    fn test_viscoelastic_relaxation_large_t() {
        let mat = PyViscoelasticMaterial::default();
        let einf = mat.relaxation_modulus(1e12);
        assert!((einf - mat.equilibrium_modulus).abs() < 1.0);
    }

    #[test]
    fn test_viscoelastic_glassy_modulus() {
        let mat = PyViscoelasticMaterial::new(1e6, vec![2e6, 1e6], vec![0.1, 1.0]);
        assert!((mat.glassy_modulus() - 4e6).abs() < 1.0);
    }

    #[test]
    fn test_viscoelastic_default() {
        let mat = PyViscoelasticMaterial::default();
        assert!(mat.prony_moduli.len() == mat.prony_times.len());
    }

    // --- PyDamageMaterial ---

    #[test]
    fn test_damage_new() {
        let mat = PyDamageMaterial::new(1000.0, 1e8, 200e9);
        assert_eq!(mat.critical_energy_release, 1000.0);
    }

    #[test]
    fn test_damage_zero_strain_energy() {
        let mat = PyDamageMaterial::default();
        assert_eq!(mat.damage_variable(0.0), 0.0);
    }

    #[test]
    fn test_damage_below_critical() {
        let mat = PyDamageMaterial::new(1000.0, 1e8, 200e9);
        assert_eq!(mat.damage_variable(500.0), 0.0);
    }

    #[test]
    fn test_damage_above_critical() {
        let mat = PyDamageMaterial::new(100.0, 1e8, 200e9);
        let d = mat.damage_variable(1000.0);
        assert!(d > 0.0 && d <= 1.0);
    }

    #[test]
    fn test_damage_effective_modulus_undamaged() {
        let mat = PyDamageMaterial::new(1000.0, 1e8, 200e9);
        let e = mat.effective_modulus(0.0);
        assert!((e - 200e9).abs() < 1.0);
    }

    #[test]
    fn test_damage_default() {
        let mat = PyDamageMaterial::default();
        assert!(mat.undamaged_modulus > 0.0);
    }

    // --- PyCreepMaterial ---

    #[test]
    fn test_creep_new() {
        let mat = PyCreepMaterial::new(1e-26, 5.0, 15_000.0);
        assert_eq!(mat.norton_exponent, 5.0);
    }

    #[test]
    fn test_creep_rate_positive() {
        let mat = PyCreepMaterial::default();
        let rate = mat.creep_rate(100e6, 800.0);
        assert!(rate >= 0.0);
    }

    #[test]
    fn test_creep_rate_increases_with_stress() {
        let mat = PyCreepMaterial::default();
        let r1 = mat.creep_rate(100e6, 800.0);
        let r2 = mat.creep_rate(200e6, 800.0);
        assert!(r2 > r1);
    }

    #[test]
    fn test_creep_rate_decreases_with_temperature_drop() {
        let mat = PyCreepMaterial::default();
        let r_hot = mat.creep_rate(100e6, 1000.0);
        let r_cold = mat.creep_rate(100e6, 500.0);
        assert!(r_hot > r_cold);
    }

    #[test]
    fn test_creep_reference_rate() {
        let mat = PyCreepMaterial::default();
        let r = mat.reference_creep_rate(100e6);
        assert!(r >= 0.0);
    }

    // --- PyAcousticMaterial ---

    #[test]
    fn test_acoustic_new() {
        let mat = PyAcousticMaterial::new(2.2e9, 998.0);
        assert_eq!(mat.density, 998.0);
    }

    #[test]
    fn test_acoustic_speed_of_sound_water() {
        let mat = PyAcousticMaterial::default();
        let c = mat.speed_of_sound();
        // Water: c ≈ 1484 m/s
        assert!((c - 1484.0).abs() < 10.0);
    }

    #[test]
    fn test_acoustic_impedance_positive() {
        let mat = PyAcousticMaterial::default();
        assert!(mat.impedance() > 0.0);
    }

    #[test]
    fn test_acoustic_reflection_same_material() {
        let mat = PyAcousticMaterial::default();
        let r = mat.reflection_coefficient(&mat.clone());
        assert!(r.abs() < 1e-10);
    }

    #[test]
    fn test_acoustic_default() {
        let mat = PyAcousticMaterial::default();
        assert!(mat.bulk_modulus > 0.0);
    }

    // --- PyCompositeMaterial ---

    #[test]
    fn test_composite_new() {
        let mat = PyCompositeMaterial::new(70e9, 200e9, 0.6);
        assert_eq!(mat.volume_fraction1, 0.6);
    }

    #[test]
    fn test_composite_voigt_bounds() {
        let mat = PyCompositeMaterial::default();
        let v = mat.voigt_modulus(0.6, 70e9, 200e9);
        assert!((70e9..=200e9).contains(&v));
    }

    #[test]
    fn test_composite_reuss_bounds() {
        let mat = PyCompositeMaterial::default();
        let r = mat.reuss_modulus(0.6, 70e9, 200e9);
        assert!((70e9..=200e9).contains(&r));
    }

    #[test]
    fn test_composite_voigt_ge_reuss() {
        let mat = PyCompositeMaterial::default();
        let v = mat.voigt_modulus(0.6, 70e9, 200e9);
        let r = mat.reuss_modulus(0.6, 70e9, 200e9);
        assert!(v >= r);
    }

    #[test]
    fn test_composite_hill_between_voigt_reuss() {
        let mat = PyCompositeMaterial::default();
        let h = mat.hill_modulus();
        let v = mat.voigt_modulus(mat.volume_fraction1, mat.modulus1, mat.modulus2);
        let r = mat.reuss_modulus(mat.volume_fraction1, mat.modulus1, mat.modulus2);
        assert!(h >= r && h <= v);
    }

    #[test]
    fn test_composite_hashin_shtrikman_lower_positive() {
        let mat = PyCompositeMaterial::default();
        let hs = mat.hashin_shtrikman_lower();
        assert!(hs > 0.0);
    }

    #[test]
    fn test_composite_pure_phase1() {
        let mat = PyCompositeMaterial::new(70e9, 200e9, 1.0);
        let v = mat.voigt_modulus(1.0, 70e9, 200e9);
        assert!((v - 70e9).abs() < 1.0);
    }

    #[test]
    fn test_composite_default() {
        let mat = PyCompositeMaterial::default();
        assert!(mat.modulus1 > 0.0 && mat.modulus2 > 0.0);
    }

    // --- register helper ---

    #[test]
    fn test_register_materials_module_no_panic() {
        register_materials_module("materials");
    }
}
