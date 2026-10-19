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

use super::bussi::{bussi_decay, sample_alpha2};
use super::functions::Thermostat;
use super::types::CsvrThermostat;

impl Thermostat for CsvrThermostat {
    /// One CSVR step: draw `K'` from the exact Bussi update (eq. A7 of
    /// Bussi, Donadio & Parrinello 2007) and scale every velocity by
    /// `sqrt(K'/K)`.  `N_f = 3·N` and `K̄ = N_f·k_B·T/2`.
    ///
    /// A system with (numerically) zero kinetic energy cannot be rescaled and
    /// is left untouched, as is a call with a non-positive `boltzmann_k` or a
    /// negative / non-finite `target_temp`.
    fn apply(&mut self, atoms: &mut AtomSet, target_temp: f64, dt: f64, boltzmann_k: f64) {
        let n = atoms.len();
        if n == 0 {
            return;
        }
        let current_temp = atoms.temperature(boltzmann_k);
        if !current_temp.is_finite() || current_temp < 1e-30 {
            return;
        }
        if !target_temp.is_finite() || target_temp < 0.0 {
            return;
        }
        let kinetic_energy = atoms.kinetic_energy();
        let ndof = (3 * n) as f64;
        let target_ke = 0.5 * ndof * boltzmann_k * target_temp;
        let c = bussi_decay(dt, self.tau);
        let scale_sq = sample_alpha2(&mut self.rng_state, kinetic_energy, target_ke, ndof, c);
        if !scale_sq.is_finite() || scale_sq <= 0.0 {
            return;
        }
        let scale = scale_sq.sqrt();
        for v in &mut atoms.velocities {
            *v *= scale;
        }
    }
}
