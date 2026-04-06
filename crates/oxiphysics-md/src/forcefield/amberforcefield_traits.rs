//! # AmberForceField - Trait Implementations
//!
//! This module contains trait implementations for `AmberForceField`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `ForceField`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::PeriodicBox;

use super::functions::ForceField;
#[allow(unused_imports)]
use super::functions::*;
use super::types::AmberForceField;

impl Default for AmberForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceField for AmberForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut energy = 0.0;
        energy += self.pair.compute_forces(atoms, pbox);
        energy += self.bonded.compute_forces(atoms, pbox);
        energy += self.angle.compute_forces(atoms, pbox);
        energy += self.dihedral.compute_forces(atoms, pbox);
        energy
    }
}
