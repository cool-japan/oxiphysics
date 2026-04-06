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
#[allow(unused_imports)]
use super::functions::*;
use super::types::BondedForceField;

impl Default for BondedForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceField for BondedForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut energy = 0.0;
        for bond in &self.bonds {
            let (dr, dist) = distance_pbc(&atoms.positions[bond.i], &atoms.positions[bond.j], pbox);
            if dist < 1e-15 {
                continue;
            }
            let delta = dist - bond.r0;
            energy += 0.5 * bond.k * delta * delta;
            let f_on_i = dr * (bond.k * delta / dist);
            atoms.forces[bond.i] += f_on_i;
            atoms.forces[bond.j] -= f_on_i;
        }
        energy
    }
}
