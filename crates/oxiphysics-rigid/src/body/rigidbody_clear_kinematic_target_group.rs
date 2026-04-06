//! # RigidBody - clear_kinematic_target_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Clear kinematic target.
    #[allow(dead_code)]
    pub fn clear_kinematic_target(&mut self) {
        self.kinematic_target = None;
    }
}
