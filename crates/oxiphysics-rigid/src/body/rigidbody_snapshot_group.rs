//! # RigidBody - snapshot_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{BodySnapshot, BodyState, BodyType};

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Capture the body state as a serializable snapshot.
    #[allow(dead_code)]
    pub fn snapshot(&self) -> BodySnapshot {
        let q = self.transform.rotation.as_vector();
        let body_type_str = match self.body_type {
            BodyType::Dynamic => "dynamic",
            BodyType::Kinematic => "kinematic",
            BodyType::Static => "static",
        };
        let state_str = match self.state {
            BodyState::Active => "active",
            BodyState::Sleeping => "sleeping",
        };
        BodySnapshot {
            position: [
                self.transform.position.x,
                self.transform.position.y,
                self.transform.position.z,
            ],
            rotation: [q[0], q[1], q[2], q[3]],
            velocity: [self.velocity.x, self.velocity.y, self.velocity.z],
            angular_velocity: [
                self.angular_velocity.x,
                self.angular_velocity.y,
                self.angular_velocity.z,
            ],
            mass: self.mass,
            body_type: body_type_str,
            state: state_str,
            linear_damping: self.linear_damping,
            angular_damping: self.angular_damping,
            gravity_scale: self.gravity_scale,
        }
    }
}
