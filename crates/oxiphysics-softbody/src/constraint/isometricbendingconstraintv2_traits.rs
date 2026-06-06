//! # IsometricBendingConstraintV2 - Trait Implementations
//!
//! This module contains trait implementations for `IsometricBendingConstraintV2`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

use super::functions::SoftConstraint;
use super::types::IsometricBendingConstraintV2;

impl SoftConstraint for IsometricBendingConstraintV2 {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let ps: [Vec3; 4] = std::array::from_fn(|k| particles[self.indices[k]].position);
        let ws: [Real; 4] = std::array::from_fn(|k| particles[self.indices[k]].inverse_mass);
        let mut c = 0.0;
        for axis in 0..3 {
            let coords: [Real; 4] = std::array::from_fn(|k| match axis {
                0 => ps[k].x,
                1 => ps[k].y,
                _ => ps[k].z,
            });
            for i in 0..4 {
                for j in 0..4 {
                    c += self.q_matrix[i * 4 + j] * coords[i] * coords[j];
                }
            }
        }
        if c.abs() < 1e-14 {
            return;
        }
        let mut grads = [Vec3::zeros(); 4];
        for (k, grad_k) in grads.iter_mut().enumerate() {
            for axis in 0..3 {
                let mut g_k = 0.0;
                for (j, p_j) in ps.iter().enumerate() {
                    let pj_coord = match axis {
                        0 => p_j.x,
                        1 => p_j.y,
                        _ => p_j.z,
                    };
                    g_k += 2.0 * self.q_matrix[k * 4 + j] * pj_coord;
                }
                match axis {
                    0 => grad_k.x += g_k,
                    1 => grad_k.y += g_k,
                    _ => grad_k.z += g_k,
                }
            }
        }
        let w_sum: Real = (0..4).map(|k| ws[k] * grads[k].norm_squared()).sum();
        if w_sum < 1e-14 {
            return;
        }
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let delta_lambda =
            (-c * self.stiffness - alpha_tilde * self.lambda) / (w_sum + alpha_tilde);
        self.lambda += delta_lambda;
        for k in 0..4 {
            particles[self.indices[k]].position += grads[k] * (delta_lambda * ws[k]);
        }
    }
}
