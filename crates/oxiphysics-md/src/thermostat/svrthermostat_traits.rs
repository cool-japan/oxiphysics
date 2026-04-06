//! # SvrThermostat - Trait Implementations
//!
//! This module contains trait implementations for `SvrThermostat`.
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
use super::types::SvrThermostat;

impl Thermostat for SvrThermostat {
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        let n = atoms.len();
        if n == 0 {
            return;
        }
        let current_ke = atoms.kinetic_energy();
        if current_ke < 1e-30 {
            return;
        }
        let target_ke = 0.5 * (self.n_dof as f64) * boltzmann_k * target_temp;
        let scale = self.rescale_factor(current_ke, target_ke, dt);
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }
}
