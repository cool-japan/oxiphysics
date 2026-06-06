//! # VanDerWaalsEos - Trait Implementations
//!
//! This module contains trait implementations for `VanDerWaalsEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
use super::types::VanDerWaalsEos;

impl EquationOfState for VanDerWaalsEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_at_temperature(density, 300.0)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let h = density * 1e-6;
        let dp = self.pressure(density + h) - self.pressure(density - h);
        let drho = 2.0 * h;
        (dp / drho).max(0.0).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        let t = 300.0;
        let mut lo = 1e-5_f64;
        let mut hi = 10000.0_f64;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.pressure_at_temperature(mid, t) < pressure {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}
