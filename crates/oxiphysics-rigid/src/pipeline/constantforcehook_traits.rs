//! # ConstantForceHook - Trait Implementations
//!
//! This module contains trait implementations for `ConstantForceHook`.
//!
//! ## Implemented Traits
//!
//! - `PreStepHook`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::PreStepHook;
#[allow(unused_imports)]
use super::functions::*;
use super::types::{BodySnapshot, ConstantForceHook};

impl PreStepHook for ConstantForceHook {
    fn pre_step(&self, bodies: &mut [BodySnapshot], _sim_time: f64, dt: f64) {
        for b in bodies.iter_mut() {
            if b.active && b.inv_mass > 0.0 {
                b.velocity[0] += self.force[0] * b.inv_mass * dt;
                b.velocity[1] += self.force[1] * b.inv_mass * dt;
                b.velocity[2] += self.force[2] * b.inv_mass * dt;
            }
        }
    }
}
