//! # IsometricBendingConstraint - Trait Implementations
//!
//! This module contains trait implementations for `IsometricBendingConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

use super::functions::SoftConstraint;
#[allow(unused_imports)]
use super::functions::*;
use super::types::IsometricBendingConstraint;

impl SoftConstraint for IsometricBendingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2, i3] = self.indices;
        let positions = [
            particles[i0].position,
            particles[i1].position,
            particles[i2].position,
            particles[i3].position,
        ];
        let indices = [i0, i1, i2, i3];
        let mut c_val = 0.0_f64;
        for a in 0..4 {
            for b in 0..4 {
                c_val += self.q_matrix[a * 4 + b] * positions[a].dot(&positions[b]);
            }
        }
        c_val *= self.stiffness;
        if c_val.abs() < 1e-14 {
            return;
        }
        let mut grads = [Vec3::zeros(); 4];
        for a in 0..4 {
            for b in 0..4 {
                grads[a] += self.q_matrix[a * 4 + b] * positions[b];
            }
            grads[a] *= 2.0 * self.stiffness;
        }
        let mut w_sum = 0.0_f64;
        for a in 0..4 {
            w_sum += particles[indices[a]].inverse_mass * grads[a].norm_squared();
        }
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c_val - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        for a in 0..4 {
            particles[indices[a]].position +=
                grads[a] * (delta_lambda * particles[indices[a]].inverse_mass);
        }
    }
}
