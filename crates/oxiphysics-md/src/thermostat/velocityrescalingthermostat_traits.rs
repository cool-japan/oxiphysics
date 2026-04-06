//! # VelocityRescalingThermostat - Trait Implementations
//!
//! This module contains trait implementations for `VelocityRescalingThermostat`.
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
use super::types::VelocityRescalingThermostat;

impl Thermostat for VelocityRescalingThermostat {
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, _dt: f64, boltzmann_k: f64) {
        let current_temp = atoms.temperature(boltzmann_k);
        if current_temp < 1e-30 || target_temp < 1e-30 {
            return;
        }
        let scale = (target_temp / current_temp).sqrt();
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }
}
