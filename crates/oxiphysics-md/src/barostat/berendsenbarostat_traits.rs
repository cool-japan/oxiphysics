//! # BerendsenBarostat - Trait Implementations
//!
//! This module contains trait implementations for `BerendsenBarostat`.
//!
//! ## Implemented Traits
//!
//! - `Barostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;

use super::functions::Barostat;
use super::types::BerendsenBarostat;

impl Barostat for BerendsenBarostat {
    fn apply(
        &mut self,
        atoms: &mut AtomSet,
        pbox: &mut PeriodicBox,
        _target_pressure: f64,
        dt: f64,
    ) {
        let volume = pbox.volume();
        if volume < 1e-30 {
            return;
        }
        let ke = atoms.kinetic_energy();
        let p_kinetic = 2.0 * ke / (3.0 * volume);
        let mut virial = 0.0;
        for i in 0..atoms.len() {
            virial += atoms.positions[i].dot(&atoms.forces[i]);
        }
        let p_current = p_kinetic + virial / (3.0 * volume);
        let mu = self.scale_factor(p_current, dt);
        pbox.dims *= mu;
        for pos in &mut atoms.positions {
            *pos *= mu;
        }
    }
}
