//! # NoseHooverThermostat - Trait Implementations
//!
//! This module contains trait implementations for `NoseHooverThermostat`.
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
use super::types::NoseHooverThermostat;

impl Thermostat for NoseHooverThermostat {
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        let n = atoms.len();
        if n == 0 {
            return;
        }
        let ndof = (3 * n) as f64;
        let ke = atoms.kinetic_energy();
        let dxi_dt = (2.0 * ke - ndof * boltzmann_k * target_temp) / self.q;
        self.xi += dxi_dt * dt;
        let scale = (-self.xi * dt).exp();
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }
}
