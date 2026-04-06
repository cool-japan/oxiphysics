// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Materials Library.

#![allow(missing_docs)]

// ===========================================================================
// Materials Library
// ===========================================================================

/// Classification of material behaviour.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum MaterialClass {
    /// Linear-elastic (Hookean) material.
    Elastic,
    /// Elasto-plastic material (von Mises yield criterion).
    Plastic,
    /// Viscous fluid material.
    Viscous,
}

/// A material definition holding both elastic and optional plastic properties.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyMaterial {
    /// Human-readable name.
    pub name: String,
    /// Classification.
    pub class: MaterialClass,
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless, -1 < ν < 0.5).
    pub poisson_ratio: f64,
    /// Mass density (kg/m^3).
    pub density: f64,
    /// Yield stress for plastic materials (Pa); ignored for elastic.
    pub yield_stress: Option<f64>,
    /// Isotropic hardening modulus (Pa); used after yield.
    pub hardening_modulus: Option<f64>,
    /// Dynamic viscosity (Pa·s) for viscous fluids.
    pub viscosity: Option<f64>,
}

impl PyMaterial {
    /// Create a generic elastic material.
    pub fn elastic(name: &str, young: f64, poisson: f64, density: f64) -> Self {
        Self {
            name: name.to_owned(),
            class: MaterialClass::Elastic,
            young_modulus: young,
            poisson_ratio: poisson,
            density,
            yield_stress: None,
            hardening_modulus: None,
            viscosity: None,
        }
    }

    /// Create an elasto-plastic material.
    pub fn plastic(
        name: &str,
        young: f64,
        poisson: f64,
        density: f64,
        yield_stress: f64,
        hardening: f64,
    ) -> Self {
        Self {
            name: name.to_owned(),
            class: MaterialClass::Plastic,
            young_modulus: young,
            poisson_ratio: poisson,
            density,
            yield_stress: Some(yield_stress),
            hardening_modulus: Some(hardening),
            viscosity: None,
        }
    }

    /// Create a viscous fluid material.
    pub fn viscous_fluid(name: &str, density: f64, viscosity: f64) -> Self {
        Self {
            name: name.to_owned(),
            class: MaterialClass::Viscous,
            young_modulus: 0.0,
            poisson_ratio: 0.0,
            density,
            yield_stress: None,
            hardening_modulus: None,
            viscosity: Some(viscosity),
        }
    }

    /// Return steel (structural carbon steel, approximate).
    pub fn steel() -> Self {
        Self::elastic("Steel", 200.0e9, 0.3, 7850.0)
    }

    /// Return aluminium alloy 6061.
    pub fn aluminium() -> Self {
        Self::elastic("Aluminium 6061", 69.0e9, 0.33, 2700.0)
    }

    /// Return natural rubber (vulcanised, approximate).
    ///
    /// Very low Young's modulus (~0.05 GPa) and near-incompressible (ν≈0.49).
    pub fn rubber() -> Self {
        Self::elastic("Natural Rubber", 0.05e9, 0.49, 920.0)
    }

    /// Return concrete (normal-weight, approximate).
    pub fn concrete() -> Self {
        Self::elastic("Concrete", 30.0e9, 0.2, 2400.0)
    }

    /// Return titanium alloy Ti-6Al-4V.
    pub fn titanium() -> Self {
        Self::elastic("Titanium Ti-6Al-4V", 114.0e9, 0.34, 4430.0)
    }

    /// Return water at 20 °C.
    pub fn water() -> Self {
        Self::viscous_fluid("Water (20°C)", 998.2, 1.002e-3)
    }

    /// Return air at standard conditions.
    pub fn air() -> Self {
        Self::viscous_fluid("Air (20°C, 1 atm)", 1.204, 1.81e-5)
    }

    /// Shear modulus G = E / (2(1 + ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }

    /// Bulk modulus K = E / (3(1 - 2ν)).
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Lamé's first parameter λ.
    pub fn lame_lambda(&self) -> f64 {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Speed of sound (longitudinal, P-wave) in the material.
    pub fn p_wave_speed(&self) -> f64 {
        let k = self.bulk_modulus();
        let g = self.shear_modulus();
        let m_modulus = k + 4.0 / 3.0 * g;
        if self.density > 0.0 {
            (m_modulus / self.density).sqrt()
        } else {
            0.0
        }
    }
}
