//! # ShapeMatchingConstraint - Trait Implementations
//!
//! This module contains trait implementations for `ShapeMatchingConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
use super::types::ShapeMatchingConstraint;

impl SoftConstraint for ShapeMatchingConstraint {
    /// Move every unpinned particle by `stiffness` times the vector to its goal
    /// position `R q_i + com`, where `com` is the current centre of mass and `R`
    /// the rotation of the polar decomposition of the cross-covariance matrix
    /// (see [`ShapeMatchingConstraint`]).
    ///
    /// All goals are derived from the positions at the start of the call.  The
    /// constraint does nothing for an empty group, when an index lies outside
    /// `particles`, or when the centre of mass is not finite (one broken
    /// particle must not spread to the rest of the group).
    fn project(&mut self, particles: &mut [SoftParticle], _dt_sub: Real) {
        if self.indices.is_empty() || self.indices.iter().any(|&idx| idx >= particles.len()) {
            return;
        }
        let com = self.current_com(particles);
        if com.iter().any(|coordinate| !coordinate.is_finite()) {
            return;
        }
        let rotation = self.compute_rotation(particles, &com);
        for (&idx, q) in self.indices.iter().zip(&self.rest_relative) {
            if particles[idx].inverse_mass == 0.0 {
                continue;
            }
            let goal = rotation * q + com;
            particles[idx].position += (goal - particles[idx].position) * self.stiffness;
        }
    }
}
