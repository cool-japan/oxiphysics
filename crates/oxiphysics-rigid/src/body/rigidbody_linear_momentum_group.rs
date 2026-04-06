//! # RigidBody - linear_momentum_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Linear momentum: m * v.
    #[allow(dead_code)]
    pub fn linear_momentum(&self) -> Vec3 {
        self.velocity * self.mass
    }
}
