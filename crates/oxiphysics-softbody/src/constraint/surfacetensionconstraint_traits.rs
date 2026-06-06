//! # SurfaceTensionConstraint - Trait Implementations
//!
//! This module contains trait implementations for `SurfaceTensionConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
use super::types::SurfaceTensionConstraint;

impl SoftConstraint for SurfaceTensionConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let cross = (p1 - p0).cross(&(p2 - p0));
        let cross_len = cross.norm();
        if cross_len < 1e-12 {
            return;
        }
        let area = 0.5 * cross_len;
        let c = area - self.rest_area;
        if c.abs() < 1e-14 {
            return;
        }
        let n = cross / cross_len;
        let e01 = p1 - p0;
        let e02 = p2 - p0;
        let grad2 = 0.5 * e01.cross(&n);
        let grad1 = 0.5 * n.cross(&e02);
        let grad0 = -(grad1 + grad2);
        let ws = [
            particles[i0].inverse_mass,
            particles[i1].inverse_mass,
            particles[i2].inverse_mass,
        ];
        let grads = [grad0, grad1, grad2];
        let w_sum: Real = ws
            .iter()
            .zip(grads.iter())
            .map(|(&w, g)| w * g.norm_squared())
            .sum();
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let effective_c = c * self.gamma;
        let delta_lambda = (-effective_c - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        particles[i0].position += grad0 * (delta_lambda * ws[0]);
        particles[i1].position += grad1 * (delta_lambda * ws[1]);
        particles[i2].position += grad2 * (delta_lambda * ws[2]);
    }
}
