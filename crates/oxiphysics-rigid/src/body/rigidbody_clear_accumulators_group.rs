//! # RigidBody - clear_accumulators_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Clear force and torque accumulators without integrating them.
    pub fn clear_accumulators(&mut self) {
        self.force_accumulator = Vec3::zeros();
        self.torque_accumulator = Vec3::zeros();
    }
}
