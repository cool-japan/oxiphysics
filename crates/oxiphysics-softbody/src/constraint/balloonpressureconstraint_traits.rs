//! # BalloonPressureConstraint - Trait Implementations
//!
//! This module contains trait implementations for `BalloonPressureConstraint`.
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
use super::types::BalloonPressureConstraint;

impl SoftConstraint for BalloonPressureConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let v_contrib = Self::triangle_volume_contribution(&p0, &p1, &p2);
        let c = v_contrib - self.rest_volume;
        let n = Self::triangle_normal_area(&p0, &p1, &p2);
        let n_len = n.norm().max(1e-12);
        let n_hat = n / n_len;
        let grad_scale = 1.0 / 6.0;
        let ws = [
            particles[i0].inverse_mass,
            particles[i1].inverse_mass,
            particles[i2].inverse_mass,
        ];
        let w_sum = ws.iter().sum::<Real>() * grad_scale * grad_scale * n_len * n_len;
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda = (-c * self.pressure - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        let correction = n_hat * (grad_scale * delta_lambda);
        particles[i0].position += correction * ws[0];
        particles[i1].position += correction * ws[1];
        particles[i2].position += correction * ws[2];
    }
}
