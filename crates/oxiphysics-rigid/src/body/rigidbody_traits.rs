//! # RigidBody - Trait Implementations
//!
//! This module contains trait implementations for `RigidBody`.
//!
//! ## Implemented Traits
//!
//! - `Steppable`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;
use oxiphysics_core::{Steppable, TimeStep};

use super::rigidbody_type::RigidBody;

impl Steppable for RigidBody {
    fn step(&mut self, time_step: &TimeStep) {
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        self.integrate_forces(time_step.dt, &gravity);
        self.integrate_velocity(time_step.dt);
    }
}
