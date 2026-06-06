//! # AndersenThermostat - Trait Implementations
//!
//! This module contains trait implementations for `AndersenThermostat`.
//!
//! ## Implemented Traits
//!
//! - `Thermostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use oxiphysics_core::math::Vec3;

use super::functions::Thermostat;
use super::types::AndersenThermostat;

impl Thermostat for AndersenThermostat {
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        let n = atoms.len();
        let collision_prob = (self.nu * dt).min(1.0);
        for i in 0..n {
            if self.next_uniform() < collision_prob {
                let sigma = (boltzmann_k * target_temp / atoms.masses[i]).sqrt();
                let vx = sigma * self.next_normal();
                let vy = sigma * self.next_normal();
                let vz = sigma * self.next_normal();
                atoms.velocities[i] = Vec3::new(vx, vy, vz);
            }
        }
    }
}
