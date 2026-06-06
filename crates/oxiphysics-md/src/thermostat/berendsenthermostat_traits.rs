//! # BerendsenThermostat - Trait Implementations
//!
//! This module contains trait implementations for `BerendsenThermostat`.
//!
//! ## Implemented Traits
//!
//! - `Thermostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;

use super::functions::Thermostat;
use super::types::BerendsenThermostat;

impl Thermostat for BerendsenThermostat {
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        let current_temp = atoms.temperature(boltzmann_k);
        if current_temp < 1e-30 {
            return;
        }
        let ratio = target_temp / current_temp;
        let scale_sq = 1.0 + dt / self.tau * (ratio - 1.0);
        let scale = scale_sq.max(0.0).sqrt();
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }
}
