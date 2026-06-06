//! # BendingConstraint - Trait Implementations
//!
//! This module contains trait implementations for `BendingConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
use super::types::BendingConstraint;

impl SoftConstraint for BendingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2, i3] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let p3 = particles[i3].position;
        let current_angle = Self::compute_dihedral(&p0, &p1, &p2, &p3);
        let c = current_angle - self.rest_angle;
        if c.abs() < 1e-10 {
            return;
        }
        let e = p1 - p0;
        let e_len = e.norm();
        if e_len < 1e-12 {
            return;
        }
        let e_hat = e / e_len;
        let n1_raw = (p2 - p0).cross(&e);
        let n2_raw = (p3 - p0).cross(&e);
        let n1_len = n1_raw.norm();
        let n2_len = n2_raw.norm();
        if n1_len < 1e-12 || n2_len < 1e-12 {
            return;
        }
        let n1 = n1_raw / n1_len;
        let n2 = n2_raw / n2_len;
        let grad2 = n1 / (p2 - p0).cross(&e_hat).norm().max(1e-12);
        let grad3 = -n2 / (p3 - p0).cross(&e_hat).norm().max(1e-12);
        let w2 = particles[i2].inverse_mass;
        let w3 = particles[i3].inverse_mass;
        let w_sum = w2 * grad2.norm_squared() + w3 * grad3.norm_squared();
        if w_sum < 1e-12 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        particles[i2].position += grad2 * (delta_lambda * w2);
        particles[i3].position += grad3 * (delta_lambda * w3);
    }
}
