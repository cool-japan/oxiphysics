//! # TaitEos - Trait Implementations
//!
//! This module contains trait implementations for `TaitEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
use super::types::TaitEos;

impl EquationOfState for TaitEos {
    fn pressure(&self, density: f64) -> f64 {
        let b = self.b_coefficient();
        let ratio = density / self.reference_density;
        b * (ratio.powf(self.gamma) - 1.0) + self.reference_pressure
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let ratio = density / self.reference_density;
        self.speed_of_sound * ratio.powf((self.gamma - 1.0) * 0.5)
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        let b = self.b_coefficient();
        let x = (pressure - self.reference_pressure) / b + 1.0;
        self.reference_density * x.powf(1.0 / self.gamma)
    }
}
