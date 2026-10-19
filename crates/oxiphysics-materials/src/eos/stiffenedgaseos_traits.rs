//! # StiffenedGasEos - Trait Implementations
//!
//! This module contains trait implementations for `StiffenedGasEos`.
//!
//! ## Implemented Traits
//!
//! - `EquationOfState`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::EquationOfState;
use super::types::StiffenedGasEos;

impl EquationOfState for StiffenedGasEos {
    /// Pressure on the isentrope through the reference state `(ρ_ref, P_ref)`:
    /// `P = (P_ref + p∞) (ρ/ρ_ref)^γ - p∞`.
    ///
    /// Non-positive densities give the zero-density limit `-p∞`.
    fn pressure(&self, density: f64) -> f64 {
        let stiffened_reference = self.reference_pressure + self.p_inf;
        let compression = if density <= 0.0 {
            0.0
        } else {
            density / self.reference_density
        };
        stiffened_reference * compression.powf(self.gamma) - self.p_inf
    }
    /// Isentropic sound speed `c = √(γ (P + p∞) / ρ)`; 0 for non-positive densities.
    fn sound_speed(&self, density: f64) -> f64 {
        if density <= 0.0 {
            return 0.0;
        }
        let stiffened_pressure = self.pressure(density) + self.p_inf;
        (self.gamma * stiffened_pressure / density).sqrt()
    }
    /// Exact inverse of [`pressure`](Self::pressure):
    /// `ρ = ρ_ref ((P + p∞)/(P_ref + p∞))^(1/γ)`.
    ///
    /// Pressures at or below `-p∞` lie beyond the zero-density limit of the
    /// isentrope and give density 0.
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        let stiffened_pressure = pressure + self.p_inf;
        if stiffened_pressure <= 0.0 {
            return 0.0;
        }
        let stiffened_reference = self.reference_pressure + self.p_inf;
        self.reference_density * (stiffened_pressure / stiffened_reference).powf(1.0 / self.gamma)
    }
}
