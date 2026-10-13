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
    /// Total energy and analytic forces (back-propagation through the network
    /// times the closed-form symmetry-function derivatives, over neighbor
    /// lists); see [`FeedForwardPotential`].
    fn energy_and_forces(&self, positions: &[Vec3], species: &[u32]) -> (f64, Vec<Vec3>) {
        self.energy_and_forces_analytic(positions, species)
    }
}
