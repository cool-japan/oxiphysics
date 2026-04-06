//! # RigidBody - wake_up_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::BodyState;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Wake this body if sleeping.
    pub fn wake_up(&mut self) {
        self.state = BodyState::Active;
        self.sleep_timer = 0.0;
    }
}
