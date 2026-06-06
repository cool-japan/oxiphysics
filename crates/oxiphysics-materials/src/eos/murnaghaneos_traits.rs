//! # MurnaghanEos - Trait Implementations
//!
//! This module contains trait implementations for `MurnaghanEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
use super::types::MurnaghanEos;

impl EquationOfState for MurnaghanEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_from_volume(1.0 / density)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let v = 1.0 / density;
        let k = self.bulk_modulus(v);
        (k / density).max(0.0).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        1.0 / self.volume_from_pressure(pressure)
    }
}
