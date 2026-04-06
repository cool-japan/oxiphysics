//! # RigidBody - accumulated_torque_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Return the total accumulated torque (before clearing).
    #[allow(dead_code)]
    pub fn accumulated_torque(&self) -> Vec3 {
        self.torque_accumulator
    }
}
