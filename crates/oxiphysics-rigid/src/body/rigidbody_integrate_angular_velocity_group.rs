//! # RigidBody - integrate_angular_velocity_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::{Quat, Real, Unit};

use super::types::{BodyState, BodyType};

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Integrate angular velocity only (for split-step methods).
    #[allow(dead_code)]
    pub fn integrate_angular_velocity(&mut self, dt: Real) {
        if self.body_type == BodyType::Static || self.state == BodyState::Sleeping {
            return;
        }
        let omega = self.angular_velocity;
        let omega_len = omega.norm();
        if omega_len > 1e-10 {
            let half_angle = omega_len * dt * 0.5;
            let axis = omega / omega_len;
            let delta_rot = Quat::from_axis_angle(&Unit::new_normalize(axis), 2.0 * half_angle);
            self.transform.rotation = delta_rot * self.transform.rotation;
        }
        self.update_world_inertia();
    }
}
