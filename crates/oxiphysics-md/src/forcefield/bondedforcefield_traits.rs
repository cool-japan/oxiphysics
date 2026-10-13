//! # BondedForceField - Trait Implementations
//!
//! This module contains trait implementations for `BondedForceField`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `ForceField`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::{PeriodicBox, distance_pbc};

use super::functions::ForceField;
use super::types::BondedForceField;

impl Default for BondedForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl BondedForceField {
    /// Accumulate harmonic bond forces and return `(energy, virial)`.
    ///
    /// With the minimum-image bond vector `dr = r_j - r_i` the force on atom `i`
    /// is `k·δ·dr/|dr|` (`δ = |dr| - r0`), so each bond contributes
    /// `(r_i - r_j)·F_i = -k·δ·|dr|` to the virial.
    fn accumulate_bond_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> (f64, f64) {
        let mut energy = 0.0;
        let mut virial = 0.0;
        for bond in &self.bonds {
            let (dr, dist) = distance_pbc(&atoms.positions[bond.i], &atoms.positions[bond.j], pbox);
            if dist < 1e-15 {
                continue;
            }
            let delta = dist - bond.r0;
            energy += 0.5 * bond.k * delta * delta;
            let f_on_i = dr * (bond.k * delta / dist);
            virial -= dr.dot(&f_on_i);
            atoms.forces[bond.i] += f_on_i;
            atoms.forces[bond.j] -= f_on_i;
        }
        (energy, virial)
    }
}

impl ForceField for BondedForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        self.accumulate_bond_forces(atoms, pbox).0
    }

    fn compute_forces_virial(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> (f64, f64) {
        self.accumulate_bond_forces(atoms, pbox)
    }
}
