//! # VinetEos - Trait Implementations
//!
//! This module contains trait implementations for `VinetEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
#[allow(unused_imports)]
use super::functions::*;
use super::types::VinetEos;

impl EquationOfState for VinetEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_from_volume(1.0 / density)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let h = density * 1e-6;
        let p_hi = self.pressure_from_volume(1.0 / (density + h));
        let p_lo = self.pressure_from_volume(1.0 / (density - h));
        ((p_lo - p_hi) / (2.0 * h)).max(0.0).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        1.0 / self.volume(pressure)
    }
}
