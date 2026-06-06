//! # AreaConstraint - Trait Implementations
//!
//! This module contains trait implementations for `AreaConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
use super::types::AreaConstraint;

impl SoftConstraint for AreaConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let area = Self::compute_triangle_area(&p0, &p1, &p2);
        let c = area - self.rest_area;
        if c.abs() < 1e-14 {
            return;
        }
        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let normal = e1.cross(&e2);
        let normal_len = normal.norm();
        if normal_len < 1e-14 {
            return;
        }
        let n_hat = normal / normal_len;
        let grad0 = 0.5 * (p1 - p2).cross(&n_hat);
        let grad1 = 0.5 * (p2 - p0).cross(&n_hat);
        let grad2 = 0.5 * (p0 - p1).cross(&n_hat);
        let grads = [grad0, grad1, grad2];
        let idxs = [i0, i1, i2];
        let mut w_sum = 0.0;
        for k in 0..3 {
            w_sum += particles[idxs[k]].inverse_mass * grads[k].norm_squared();
        }
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        for k in 0..3 {
            particles[idxs[k]].position +=
                grads[k] * (delta_lambda * particles[idxs[k]].inverse_mass);
        }
    }
}
