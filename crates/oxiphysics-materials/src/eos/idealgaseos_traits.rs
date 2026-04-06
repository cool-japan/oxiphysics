//! # IdealGasEos - Trait Implementations
//!
//! This module contains trait implementations for `IdealGasEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//! - `EosWithEnergy`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{EosWithEnergy, EquationOfState};
use super::types::IdealGasEos;

impl EquationOfState for IdealGasEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_at_temperature(density, 288.15)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let p = self.pressure(density);
        if density.abs() < f64::EPSILON {
            return 0.0;
        }
        (self.gamma * p / density).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        pressure / (self.specific_gas_constant * 288.15)
    }
}

impl EosWithEnergy for IdealGasEos {
    fn pressure_energy(&self, density: f64, energy: f64) -> f64 {
        (self.gamma - 1.0) * density * energy
    }
    fn internal_energy(&self, _density: f64, temperature: f64) -> f64 {
        self.specific_gas_constant * temperature / (self.gamma - 1.0)
    }
}
