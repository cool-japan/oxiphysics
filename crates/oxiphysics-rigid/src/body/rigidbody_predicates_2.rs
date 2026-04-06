//! # RigidBody - predicates Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::BodyType;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Is this body static?
    #[allow(dead_code)]
    pub fn is_static(&self) -> bool {
        self.body_type == BodyType::Static
    }
}
