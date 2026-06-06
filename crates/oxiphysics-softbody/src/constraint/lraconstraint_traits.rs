//! # LraConstraint - Trait Implementations
//!
//! This module contains trait implementations for `LraConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
use super::types::LraConstraint;

impl SoftConstraint for LraConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let w = particles[self.particle].inverse_mass;
        if w == 0.0 {
            return;
        }
        let diff = particles[self.particle].position - self.anchor;
        let dist = diff.norm();
        if dist <= self.max_distance || dist < 1e-12 {
            return;
        }
        let c = dist - self.max_distance;
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c - alpha_tilde * self.lambda) / (w + alpha_tilde);
        self.lambda += delta_lambda;
        particles[self.particle].position += (diff / dist) * (delta_lambda * w);
    }
}
