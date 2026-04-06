//! # LangevinThermostat - Trait Implementations
//!
//! This module contains trait implementations for `LangevinThermostat`.
//!
//! ## Implemented Traits
//!
//! - `Thermostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;

use super::functions::Thermostat;
#[allow(unused_imports)]
use super::functions::*;
use super::types::LangevinThermostat;

impl Thermostat for LangevinThermostat {
    /// Apply a single Ornstein-Uhlenbeck (O) step to velocities.
    ///
    /// In the BAOAB context, the B and A steps are handled by the integrator;
    /// here we expose only the O step via the `Thermostat` interface.
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        self.apply_ou_step(atoms, target_temp, dt, boltzmann_k);
    }
}
