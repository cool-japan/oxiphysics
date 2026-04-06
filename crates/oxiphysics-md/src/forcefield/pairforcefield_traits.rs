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
#[allow(unused_imports)]
use super::functions::*;
use super::types::PairForceField;

impl Default for PairForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceField for PairForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let cutoff = self.max_cutoff();
        if cutoff <= 0.0 || atoms.is_empty() {
            return 0.0;
        }
        let cell_list = CellList::build(&atoms.positions, pbox, cutoff);
        let pairs = cell_list.neighbor_pairs(&atoms.positions, pbox);
        let mut energy = 0.0;
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
                atoms.forces[i] -= f_vec;
                atoms.forces[j] += f_vec;
            }
        }
        energy
    }
}
