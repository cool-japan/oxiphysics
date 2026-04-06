//! # StiffenedGasEos - Trait Implementations
//!
//! This module contains trait implementations for `StiffenedGasEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
#[allow(unused_imports)]
use super::functions::*;
use super::types::StiffenedGasEos;

impl EquationOfState for StiffenedGasEos {
    fn pressure(&self, density: f64) -> f64 {
        let cv = 4186.0;
        let e = cv * 288.15;
        self.pressure_from_energy(density, e)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let p = self.pressure(density);
        if density.abs() < f64::EPSILON {
            return 0.0;
        }
        (self.gamma * (p + self.p_inf) / density).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        let cv = 4186.0;
        let e = cv * 288.15;
        if (self.gamma - 1.0).abs() < f64::EPSILON {
            return self.reference_density;
        }
        (pressure + self.p_inf) / ((self.gamma - 1.0) * (e - self.reference_energy))
    }
}
