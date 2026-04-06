//! # NobleAbelEos - Trait Implementations
//!
//! This module contains trait implementations for `NobleAbelEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
#[allow(unused_imports)]
use super::functions::*;
use super::types::NobleAbelEos;

impl EquationOfState for NobleAbelEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_at_temperature(density, 300.0)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let h = density * 1e-6;
        let dp = self.pressure(density + h) - self.pressure(density - h);
        ((dp / (2.0 * h)).max(0.0)).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        self.density_at_temperature(pressure, 300.0)
    }
}
