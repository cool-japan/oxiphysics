//! # TillotsonEos - Trait Implementations
//!
//! This module contains trait implementations for `TillotsonEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
#[allow(unused_imports)]
use super::functions::*;
use super::types::TillotsonEos;

impl EquationOfState for TillotsonEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_from_energy(density, 0.0)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let h = density * 1e-6;
        let dp = self.pressure(density + h) - self.pressure(density - h);
        let drho = 2.0 * h;
        (dp / drho).max(0.0).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        let mut lo = self.rho0 * 0.1;
        let mut hi = self.rho0 * 5.0;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.pressure(mid) < pressure {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}
