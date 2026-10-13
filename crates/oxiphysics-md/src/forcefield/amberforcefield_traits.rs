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

    /// Sum of the component virials: the pair and bond terms use their
    /// analytic pair sums, the angle and dihedral terms the energy-based
    /// default of [`ForceField::compute_forces_virial`].
    fn compute_forces_virial(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> (f64, f64) {
        let (e_pair, w_pair) = self.pair.compute_forces_virial(atoms, pbox);
        let (e_bond, w_bond) = self.bonded.compute_forces_virial(atoms, pbox);
        let (e_angle, w_angle) = self.angle.compute_forces_virial(atoms, pbox);
        let (e_dihedral, w_dihedral) = self.dihedral.compute_forces_virial(atoms, pbox);
        (
            e_pair + e_bond + e_angle + e_dihedral,
            w_pair + w_bond + w_angle + w_dihedral,
        )
    }
}
