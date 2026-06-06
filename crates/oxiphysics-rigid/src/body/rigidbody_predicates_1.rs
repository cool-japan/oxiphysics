//! # RigidBody - predicates Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::BodyType;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Is this body kinematic?
    pub fn is_kinematic(&self) -> bool {
        self.body_type == BodyType::Kinematic
    }
}
