//! # MotorConstraint - Trait Implementations
//!
//! This module contains trait implementations for `MotorConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
use super::types::MotorConstraint;

impl SoftConstraint for MotorConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let w1 = particles[self.i].inverse_mass;
        let w2 = particles[self.j].inverse_mass;
        let w_sum = w1 + w2;
        if w_sum < 1e-14 {
            return;
        }
        let diff = particles[self.i].position - particles[self.j].position;
        let dist = diff.norm();
        if dist < 1e-12 {
            return;
        }
        let c = dist - self.target_length;
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        let correction = diff / dist * delta_lambda;
        particles[self.i].position += correction * w1;
        particles[self.j].position -= correction * w2;
    }
}
