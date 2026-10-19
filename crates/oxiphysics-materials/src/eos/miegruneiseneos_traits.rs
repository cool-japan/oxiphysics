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
    /// Sound speed `c = √(dP/dρ)` of the zero-energy pressure curve `P(ρ)`.
    ///
    /// The slope is evaluated analytically (see
    /// [`MieGruneisenEos::pressure_density_derivative`]), so `c² = P_H′ − Γ₀ (e_H + ρ e_H′)`
    /// has the dimension of a squared velocity and `c(ρ₀) = c₀`.  For
    /// `ρ ≤ ρ₀` (where the curve is identically zero) the reference sound speed
    /// `c₀` is returned; where the curve softens (`dP/dρ ≤ 0`) the result is 0
    /// and beyond the Hugoniot pole it is `+∞`.
    fn sound_speed(&self, density: f64) -> f64 {
        if density <= self.rho0 {
            return self.c0;
        }
        let slope = self.pressure_density_derivative(density);
        if slope > 0.0 { slope.sqrt() } else { 0.0 }
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
