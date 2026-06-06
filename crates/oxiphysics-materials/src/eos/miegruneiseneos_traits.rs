//! # MieGruneisenEos - Trait Implementations
//!
//! This module contains trait implementations for `MieGruneisenEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
use super::types::MieGruneisenEos;

impl EquationOfState for MieGruneisenEos {
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_from_energy(density, 0.0)
    }
    fn sound_speed(&self, density: f64) -> f64 {
        let mu = density / self.rho0 - 1.0;
        let denom = (1.0 - self.s * mu).powi(2);
        if denom.abs() < f64::EPSILON {
            return self.c0;
        }
        let ph = self.hugoniot_pressure(density);
        let dp_drho =
            self.c0 * self.c0 * (1.0 + 2.0 * self.s * mu) / (denom * denom) + self.gamma_0 * ph;
        dp_drho.max(0.0).sqrt()
    }
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        let mut lo = self.rho0 * 0.5;
        let mut hi = self.rho0 * 10.0;
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
