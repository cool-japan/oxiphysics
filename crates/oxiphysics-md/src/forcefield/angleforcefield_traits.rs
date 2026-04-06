//! # AngleForceField - Trait Implementations
//!
//! This module contains trait implementations for `AngleForceField`.
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
use super::types::AngleForceField;

impl Default for AngleForceField {
    fn default() -> Self {
        Self::new()
    }
}

impl ForceField for AngleForceField {
    fn compute_forces(&self, atoms: &mut AtomSet, pbox: &PeriodicBox) -> f64 {
        let mut energy = 0.0;
        for angle in &self.angles {
            let (rji, dji) =
                distance_pbc(&atoms.positions[angle.j], &atoms.positions[angle.i], pbox);
            let (rjk, djk) =
                distance_pbc(&atoms.positions[angle.j], &atoms.positions[angle.k], pbox);
            if dji < 1e-15 || djk < 1e-15 {
                continue;
            }
            let rji_n = rji / dji;
            let rjk_n = rjk / djk;
            let cos_theta = rji_n.dot(&rjk_n).clamp(-1.0, 1.0);
            let theta = cos_theta.acos();
            let delta = theta - angle.theta0;
            energy += 0.5 * angle.k_theta * delta * delta;
            let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-15);
            let d_energy = angle.k_theta * delta;
            let fi = (rjk_n - rji_n * cos_theta) / (dji * sin_theta) * (-d_energy);
            let fk = (rji_n - rjk_n * cos_theta) / (djk * sin_theta) * (-d_energy);
            atoms.forces[angle.i] += fi;
            atoms.forces[angle.k] += fk;
            atoms.forces[angle.j] -= fi + fk;
        }
        energy
    }
}
