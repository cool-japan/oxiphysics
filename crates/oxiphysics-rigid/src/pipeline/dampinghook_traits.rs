//! # DampingHook - Trait Implementations
//!
//! This module contains trait implementations for `DampingHook`.
//!
//! ## Implemented Traits
//!
//! - `PreStepHook`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::PreStepHook;
#[allow(unused_imports)]
use super::functions::*;
use super::types::{BodySnapshot, DampingHook};

impl PreStepHook for DampingHook {
    fn pre_step(&self, bodies: &mut [BodySnapshot], _sim_time: f64, dt: f64) {
        let factor = 1.0 - (self.coeff * dt).min(1.0);
        for b in bodies.iter_mut() {
            if b.active && b.inv_mass > 0.0 {
                b.velocity[0] *= factor;
                b.velocity[1] *= factor;
                b.velocity[2] *= factor;
            }
        }
    }
}
