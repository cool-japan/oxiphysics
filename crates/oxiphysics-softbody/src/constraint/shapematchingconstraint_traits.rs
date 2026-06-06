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
use oxiphysics_core::math::{Real, Vec3};

use super::functions::SoftConstraint;
use super::types::ShapeMatchingConstraint;

impl SoftConstraint for ShapeMatchingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], _dt_sub: Real) {
        if self.indices.is_empty() {
            return;
        }
        let com = self.current_com(particles);
        let rot = self.compute_rotation(particles);
        for (k, &idx) in self.indices.iter().enumerate() {
            if particles[idx].inverse_mass == 0.0 {
                continue;
            }
            let q = &self.rest_relative[k];
            let goal = Vec3::new(
                rot[0] * q.x + rot[1] * q.y + rot[2] * q.z,
                rot[3] * q.x + rot[4] * q.y + rot[5] * q.z,
                rot[6] * q.x + rot[7] * q.y + rot[8] * q.z,
            ) + com;
            particles[idx].position += (goal - particles[idx].position) * self.stiffness;
        }
    }
}
