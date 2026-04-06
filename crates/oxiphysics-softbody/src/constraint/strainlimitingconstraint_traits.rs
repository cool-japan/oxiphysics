//! # StrainLimitingConstraint - Trait Implementations
//!
//! This module contains trait implementations for `StrainLimitingConstraint`.
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
use super::types::StrainLimitingConstraint;

impl SoftConstraint for StrainLimitingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let d1 = p1 - p0;
        let d2 = p2 - p0;
        let e1_len = d1.norm();
        if e1_len < 1e-14 {
            return;
        }
        let u = d1 / e1_len;
        let normal = d1.cross(&d2);
        let n_len = normal.norm();
        if n_len < 1e-14 {
            return;
        }
        let v = normal.cross(&d1).normalize();
        let ax = d1.dot(&u);
        let ay = d1.dot(&v);
        let bx = d2.dot(&u);
        let by = d2.dot(&v);
        let f00 = ax * self.inv_rest[0] + bx * self.inv_rest[2];
        let f01 = ax * self.inv_rest[1] + bx * self.inv_rest[3];
        let f10 = ay * self.inv_rest[0] + by * self.inv_rest[2];
        let f11 = ay * self.inv_rest[1] + by * self.inv_rest[3];
        let ftf00 = f00 * f00 + f10 * f10;
        let ftf01 = f00 * f01 + f10 * f11;
        let ftf11 = f01 * f01 + f11 * f11;
        let trace = ftf00 + ftf11;
        let det = ftf00 * ftf11 - ftf01 * ftf01;
        let disc = (trace * trace - 4.0 * det).max(0.0);
        let sqrt_disc = disc.sqrt();
        let s1_sq = (trace + sqrt_disc) / 2.0;
        let s2_sq = (trace - sqrt_disc) / 2.0;
        let s1 = s1_sq.max(0.0).sqrt();
        let s2 = s2_sq.max(0.0).sqrt();
        let needs_correction = s1 > self.max_stretch
            || s2 > self.max_stretch
            || s1 < self.min_compression
            || s2 < self.min_compression;
        if !needs_correction {
            return;
        }
        let s1_clamped = s1.clamp(self.min_compression, self.max_stretch);
        let s2_clamped = s2.clamp(self.min_compression, self.max_stretch);
        let c_val = (s1 - s1_clamped).abs() + (s2 - s2_clamped).abs();
        if c_val < 1e-14 {
            return;
        }
        let centroid = (p0 + p1 + p2) / 3.0;
        let correction_factor = c_val * 0.3;
        let alpha_tilde = self.compliance / (dt_sub * dt_sub);
        let effective_factor = correction_factor / (1.0 + alpha_tilde);
        for &idx in &[i0, i1, i2] {
            if particles[idx].inverse_mass == 0.0 {
                continue;
            }
            let to_centroid = centroid - particles[idx].position;
            particles[idx].position += to_centroid * effective_factor;
        }
    }
}
