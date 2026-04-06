//! # VolumeConstraint - Trait Implementations
//!
//! This module contains trait implementations for `VolumeConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
#[allow(unused_imports)]
use super::functions::*;
use super::types::VolumeConstraint;

impl SoftConstraint for VolumeConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2, i3] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let p3 = particles[i3].position;
        let vol = Self::compute_tet_volume(&p0, &p1, &p2, &p3);
        let c = vol - self.rest_volume;
        if c.abs() < 1e-14 {
            return;
        }
        let grad0 = (p1 - p2).cross(&(p3 - p2)) / 6.0;
        let grad1 = (p2 - p0).cross(&(p3 - p0)) / 6.0;
        let grad2 = (p0 - p1).cross(&(p3 - p1)) / 6.0;
        let grad3 = (p1 - p0).cross(&(p2 - p0)) / 6.0;
        let grads = [grad0, grad1, grad2, grad3];
        let idxs = [i0, i1, i2, i3];
        let mut w_sum = 0.0;
        for k in 0..4 {
            w_sum += particles[idxs[k]].inverse_mass * grads[k].norm_squared();
        }
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        for k in 0..4 {
            particles[idxs[k]].position +=
                grads[k] * (delta_lambda * particles[idxs[k]].inverse_mass);
        }
    }
}
