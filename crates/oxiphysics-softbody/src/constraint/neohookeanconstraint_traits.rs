//! # NeoHookeanConstraint - Trait Implementations
//!
//! This module contains trait implementations for `NeoHookeanConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::{Real, Vec3};

use super::functions::SoftConstraint;
#[allow(unused_imports)]
use super::functions::*;
use super::types::NeoHookeanConstraint;

impl SoftConstraint for NeoHookeanConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let f = self.deformation_gradient(particles);
        let j = Self::compute_j(f);
        let i1 = Self::compute_i1(f);
        let i1_safe = (i1 / 3.0).max(0.0);
        let c_dev = i1_safe.sqrt() - 1.0;
        let c_vol = j - 1.0;
        let [i0, i1_idx, i2, i3] = self.indices;
        let ws = [
            particles[i0].inverse_mass,
            particles[i1_idx].inverse_mass,
            particles[i2].inverse_mass,
            particles[i3].inverse_mass,
        ];
        let w_sum: Real = ws.iter().sum();
        if w_sum < 1e-14 {
            return;
        }
        let alpha_dev = self.mu / (self.rest_volume * (dt_sub * dt_sub) + 1e-30);
        let denom_dev = w_sum + alpha_dev;
        if denom_dev > 1e-14 {
            let dl = (-c_dev - alpha_dev * self.lambda_dev) / denom_dev;
            self.lambda_dev += dl;
            let dirs = [
                particles[i0].position,
                particles[i1_idx].position,
                particles[i2].position,
                particles[i3].position,
            ];
            let centroid = Vec3::new(
                dirs.iter().map(|p| p.x).sum::<Real>() / 4.0,
                dirs.iter().map(|p| p.y).sum::<Real>() / 4.0,
                dirs.iter().map(|p| p.z).sum::<Real>() / 4.0,
            );
            let idxs = [i0, i1_idx, i2, i3];
            for (k, &idx) in idxs.iter().enumerate() {
                let d = particles[idx].position - centroid;
                let d_len = d.norm().max(1e-12);
                particles[idx].position += (d / d_len) * (dl * ws[k]);
            }
        }
        let alpha_vol = self.lambda_lame / (self.rest_volume * (dt_sub * dt_sub) + 1e-30);
        let denom_vol = w_sum + alpha_vol;
        if denom_vol > 1e-14 && c_vol.abs() > 1e-12 {
            let dl = (-c_vol - alpha_vol * self.lambda_vol) / denom_vol;
            self.lambda_vol += dl;
            let idxs = [i0, i1_idx, i2, i3];
            for (k, &idx) in idxs.iter().enumerate() {
                particles[idx].position.y += dl * ws[k] * 0.25;
                let _ = k;
            }
        }
    }
}
