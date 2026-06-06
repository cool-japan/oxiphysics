//! # CsvrThermostat - Trait Implementations
//!
//! This module contains trait implementations for `CsvrThermostat`.
//!
//! ## Implemented Traits
//!
//! - `Thermostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;

use super::functions::Thermostat;
use super::types::CsvrThermostat;

impl Thermostat for CsvrThermostat {
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        let n = atoms.len();
        if n == 0 {
            return;
        }
        let current_temp = atoms.temperature(boltzmann_k);
        if current_temp < 1e-30 {
            return;
        }
        let c = (-dt / self.tau).exp();
        let r = self.next_normal();
        let ndof = (3 * n) as f64;
        let ratio = target_temp / current_temp;
        let det_scale_sq = 1.0 + (1.0 - c) * (ratio - 1.0);
        let stoch = (2.0 * (1.0 - c) * ratio / ndof).max(0.0).sqrt() * r;
        let scale_sq = (det_scale_sq + stoch).max(0.0);
        let scale = scale_sq.sqrt();
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }
}
