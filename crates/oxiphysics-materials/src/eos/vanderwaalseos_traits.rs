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
    /// Pressure on the isotherm at [`temperature`](VanDerWaalsEos::temperature).
    fn pressure(&self, density: f64) -> f64 {
        self.pressure_at_temperature(density, self.temperature)
    }
    /// Isothermal sound speed `√((∂P/∂ρ)_T)` on the isotherm at
    /// [`temperature`](VanDerWaalsEos::temperature); 0 inside the spinodal
    /// region where `(∂P/∂ρ)_T ≤ 0`.
    fn sound_speed(&self, density: f64) -> f64 {
        let slope = self.pressure_density_derivative_at_temperature(density, self.temperature);
        if slope <= 0.0 { 0.0 } else { slope.sqrt() }
    }
    /// Vapour-branch density at `pressure` on the isotherm at
    /// [`temperature`](VanDerWaalsEos::temperature), see
    /// [`VanDerWaalsEos::vapour_molar_volume`].
    fn density_from_pressure(&self, pressure: f64) -> f64 {
        self.density_at_temperature(pressure, self.temperature)
    }
}
