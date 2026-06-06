//! # ParrinelloRahmanBarostat - Trait Implementations
//!
//! This module contains trait implementations for `ParrinelloRahmanBarostat`.
//!
//! ## Implemented Traits
//!
//! - `Barostat`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;

use super::functions::Barostat;
use super::types::ParrinelloRahmanBarostat;

impl Barostat for ParrinelloRahmanBarostat {
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
        let l = volume.cbrt();
        let force_box = 3.0 * volume * (p_current - self.target_pressure);
        let a_box = force_box / self.mass;
        let v_box = a_box * dt;
        let dl = v_box * dt;
        let new_l = l + dl;
        if new_l > 0.0 {
            let scale = new_l / l;
            pbox.dims *= scale;
            for pos in &mut atoms.positions {
                *pos *= scale;
            }
        }
    }
}
