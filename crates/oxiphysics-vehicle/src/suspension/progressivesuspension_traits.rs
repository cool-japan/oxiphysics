//! # ProgressiveSuspension - Trait Implementations
//!
//! This module contains trait implementations for `ProgressiveSuspension`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `SuspensionModel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Real;

use super::functions::SuspensionModel;
#[allow(unused_imports)]
use super::functions::*;
use super::types::ProgressiveSuspension;

impl Default for ProgressiveSuspension {
    fn default() -> Self {
        Self {
            progressive_factor: 0.5,
        }
    }
}

impl SuspensionModel for ProgressiveSuspension {
    fn compute_force(
        &self,
        rest_length: Real,
        current_length: Real,
        velocity: Real,
        stiffness: Real,
        damping: Real,
    ) -> Real {
        let displacement = rest_length - current_length;
        if displacement <= 0.0 {
            return 0.0;
        }
        let ratio = if rest_length > 1e-10 {
            displacement / rest_length
        } else {
            0.0
        };
        let nonlinear_stiffness = stiffness * (1.0 + self.progressive_factor * ratio);
        let spring_force = nonlinear_stiffness * displacement;
        let damper_force = -damping * velocity;
        (spring_force + damper_force).max(0.0)
    }
}
