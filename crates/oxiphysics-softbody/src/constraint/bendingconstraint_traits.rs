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
use std::f64::consts::PI;

use super::functions::SoftConstraint;
use super::types::BendingConstraint;

/// Constraint violation (rad) below which no correction is applied.
const NEGLIGIBLE_VIOLATION: Real = 1.0e-10;

/// Smallest weighted-gradient-plus-compliance denominator that is divided by.
const MIN_DENOMINATOR: Real = 1.0e-12;

impl SoftConstraint for BendingConstraint {
    /// One XPBD update of the hinge constraint `C = theta - rest_angle`.
    ///
    /// `theta = pi - |phi|` is the dihedral angle between the face normals and
    /// `phi` the signed fold angle of the hinge (see
    /// [`BendingConstraint::rest_angle`]).  Its gradient is `-sign(phi) * d phi / dx`
    /// for all four hinge particles, so the correction restores a fold whichever
    /// way the cloth is bent, moves the hinge endpoints as well as the wings and
    /// conserves linear and angular momentum.  A flat hinge (`phi = 0`) with a
    /// creased rest angle folds towards positive `phi`.
    ///
    /// The constraint does nothing for a non-positive or non-finite `dt_sub`, an
    /// index outside `particles`, or a degenerate hinge (zero-length edge or a
    /// wing vertex on the hinge line).
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let usable = dt_sub.is_finite()
            && dt_sub > 0.0
            && self.indices.iter().all(|&index| index < particles.len());
        if !usable {
            return;
        }
        let [i0, i1, i2, i3] = self.indices;
        let Some(fold) = Self::hinge_fold(
            &particles[i0].position,
            &particles[i1].position,
            &particles[i2].position,
            &particles[i3].position,
        ) else {
            return;
        };
        let c = PI - fold.angle.abs() - self.rest_angle;
        if c.abs() < NEGLIGIBLE_VIOLATION {
            return;
        }
        // theta = pi - |phi|  =>  d theta / dx = -sign(phi) d phi / dx, with the
        // tie at phi = 0 broken towards positive phi.
        let orientation = if fold.angle < 0.0 { 1.0 } else { -1.0 };
        let gradient = fold.gradient.map(|g| g * orientation);
        let weights = self.indices.map(|index| particles[index].inverse_mass);
        let weighted_gradient_norm: Real = (0..4)
            .map(|k| weights[k] * gradient[k].norm_squared())
            .sum();
        let alpha_tilde = self.compliance.max(0.0) / (dt_sub * dt_sub);
        let denominator = weighted_gradient_norm + alpha_tilde;
        if denominator < MIN_DENOMINATOR {
            return;
        }
        let delta_lambda = (-c - alpha_tilde * self.lambda) / denominator;
        if !delta_lambda.is_finite() {
            return;
        }
        self.lambda += delta_lambda;
        for (k, &index) in self.indices.iter().enumerate() {
            particles[index].position += gradient[k] * (weights[k] * delta_lambda);
        }
    }

    fn reset_lambda(&mut self) {
        Self::reset_lambda(self);
    }
}
