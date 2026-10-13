//! # PairForceField - Trait Implementations
//!
//! This module contains trait implementations for `PairForceField`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `ForceField`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::atom::AtomSet;
use crate::neighbor::{CellList, PeriodicBox};

use super::functions::ForceField;
use super::types::PairForceField;

impl Default for PairForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl PairForceField {
    /// Accumulate pair forces and return `(energy, virial)`.
    ///
    /// For each minimum-image pair `dr = r_j - r_i` with repulsive-positive
    /// force magnitude `f = -dV/dr`, the force on `j` is `f·dr/|dr|` and the
    /// pair contributes `(r_i - r_j)·F_i = f·|dr|` to the virial.  Only
    /// separations enter, so the virial does not depend on how the atoms are
    /// wrapped into the periodic cell.
    fn accumulate_pair_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> (f64, f64) {
        let cutoff = self.max_cutoff();
        if cutoff <= 0.0 || atoms.is_empty() {
            return (0.0, 0.0);
        }
        let cell_list = CellList::build(&atoms.positions, pbox, cutoff);
        let pairs = cell_list.neighbor_pairs(&atoms.positions, pbox);
        let mut energy = 0.0;
        let mut virial = 0.0;
        for (i, j, dr, dist) in pairs {
            let ti = atoms.types[i];
            let tj = atoms.types[j];
            if let Some(pot) = self.find_potential(ti, tj)
                && dist < pot.cutoff()
                && dist > 1e-15
            {
                energy += pot.energy(dist);
                let f_mag = pot.force(dist);
                let f_vec = dr * (f_mag / dist);
                virial += dr.dot(&f_vec);
                atoms.forces[i] -= f_vec;
                atoms.forces[j] += f_vec;
            }
        }
        (energy, virial)
    }
}

impl ForceField for PairForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        self.accumulate_pair_forces(atoms, pbox).0
    }

    fn compute_forces_virial(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> (f64, f64) {
        self.accumulate_pair_forces(atoms, pbox)
    }
}
