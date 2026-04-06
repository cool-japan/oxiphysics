//! # LinearSuspension - Trait Implementations
//!
//! This module contains trait implementations for `LinearSuspension`.
//!
//! ## Implemented Traits
//!
//! - `SuspensionModel`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Real;

use super::functions::SuspensionModel;
#[allow(unused_imports)]
use super::functions::*;
use super::types::LinearSuspension;

impl SuspensionModel for LinearSuspension {
    fn compute_force(
        &self,
        rest_length: Real,
        current_length: Real,
        velocity: Real,
        stiffness: Real,
        damping: Real,
    ) -> Real {
        let displacement = rest_length - current_length;
        let spring_force = stiffness * displacement;
        let damper_force = -damping * velocity;
        (spring_force + damper_force).max(0.0)
    }
}
