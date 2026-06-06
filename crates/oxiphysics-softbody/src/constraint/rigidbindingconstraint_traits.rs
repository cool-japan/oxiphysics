//! # RigidBindingConstraint - Trait Implementations
//!
//! This module contains trait implementations for `RigidBindingConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
use super::types::RigidBindingConstraint;

impl SoftConstraint for RigidBindingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let w = particles[self.particle].inverse_mass;
        if w < 1e-14 {
            return;
        }
        let diff = particles[self.particle].position - self.target;
        let dist = diff.norm();
        if dist < 1e-12 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-dist - alpha_tilde * self.lambda) / (w + alpha_tilde);
        self.lambda += delta_lambda;
        particles[self.particle].position += (diff / dist) * (delta_lambda * w);
    }
}
