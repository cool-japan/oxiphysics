//! # BirchMurnaghan3Eos - Trait Implementations
//!
//! This module contains trait implementations for `BirchMurnaghan3Eos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
use super::types::BirchMurnaghan3Eos;

impl EquationOfState for BirchMurnaghan3Eos {
    fn pressure(&self, density: f64) -> f64 {
        let v = 1.0 / density;
        self.pressure_from_volume(v)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let v = 1.0 / density;
        let k = self.bulk_modulus(v);
        (k / density).max(0.0).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        1.0 / self.volume(pressure)
    }
}
