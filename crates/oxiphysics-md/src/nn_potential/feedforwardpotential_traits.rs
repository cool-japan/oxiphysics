//! # FeedForwardPotential - Trait Implementations
//!
//! This module contains trait implementations for `FeedForwardPotential`.
//!
//! ## Implemented Traits
//!
//! - `MlPotential`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::ml_potential::MlPotential;
use oxiphysics_core::math::Vec3;

use super::types::FeedForwardPotential;

impl MlPotential for FeedForwardPotential {
    fn energy_and_forces(&self, positions: &[Vec3], species: &[u32]) -> (f64, Vec<Vec3>) {
        let energy = self.total_energy(positions, species);
        let n = positions.len();
        let mut forces = vec![Vec3::zeros(); n];
        let h = self.fd_step;
        let mut pos_plus = positions.to_vec();
        let mut pos_minus = positions.to_vec();
        for i in 0..n {
            for a in 0..3_usize {
                pos_plus[i] = positions[i];
                pos_minus[i] = positions[i];
                pos_plus[i][a] += h;
                pos_minus[i][a] -= h;
                let e_plus = self.total_energy(&pos_plus, species);
                let e_minus = self.total_energy(&pos_minus, species);
                forces[i][a] = -(e_plus - e_minus) / (2.0 * h);
                pos_plus[i] = positions[i];
                pos_minus[i] = positions[i];
            }
        }
        (energy, forces)
    }
}
