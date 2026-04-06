//! # RigidBody - accumulated_force_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Return the total accumulated force (before clearing).
    #[allow(dead_code)]
    pub fn accumulated_force(&self) -> Vec3 {
        self.force_accumulator
    }
}
