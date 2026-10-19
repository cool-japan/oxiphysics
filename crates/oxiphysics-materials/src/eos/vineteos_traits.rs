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
use super::types::VinetEos;

impl EquationOfState for VinetEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_from_volume(1.0 / density)
    }
    /// Sound speed `c = √(K/ρ)` from the analytic isothermal bulk modulus
    /// [`VinetEos::bulk_modulus`]; 0 for non-positive densities or where the
    /// modulus is non-positive (strong expansion).
    fn sound_speed(&self, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        let bulk_modulus = self.bulk_modulus(1.0 / density);
        if bulk_modulus <= 0.0 {
            0.0
        } else {
            (bulk_modulus / density).sqrt()
        }
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        1.0 / self.volume(pressure)
    }
}
