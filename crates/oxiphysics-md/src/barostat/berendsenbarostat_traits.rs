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

use super::functions::{Barostat, virial_pressure};
use super::types::BerendsenBarostat;

impl Barostat for BerendsenBarostat {
    /// Berendsen weak-coupling step.
    ///
    /// The instantaneous pressure is `P = (2·KE + W) / (3·V)` with the pair
    /// virial `W` supplied by the caller (it is invariant under the wrapping
    /// of coordinates into the periodic cell).  The box and the positions are
    /// scaled by `μ = (1 − β·dt/τ·(P₀ − P))^(1/3)` towards the target pressure
    /// `P₀`: the `target_pressure` argument if it is finite, the barostat's own
    /// `target_pressure` otherwise.
    fn apply(
        &mut self,
        atoms: &mut AtomSet,
        pbox: &mut PeriodicBox,
        target_pressure: f64,
        virial: f64,
        dt: f64,
    ) {
        let volume = pbox.volume();
        if volume < 1e-30 {
            return;
        }
        let ke = atoms.kinetic_energy();
        let p_current = virial_pressure(volume, ke, virial);
        let p_target = if target_pressure.is_finite() {
            target_pressure
        } else {
            self.target_pressure
        };
        let mu = self.scale_factor_to(p_current, p_target, dt);
        pbox.dims *= mu;
        for pos in &mut atoms.positions {
            *pos *= mu;
        }
    }
}
