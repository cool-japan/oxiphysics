//! # RigidBody - integrate_velocity_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_core::math::{Quat, Real, Unit};

use super::types::{BodyState, BodyType};

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Integrate velocity -> position (semi-implicit Euler, second half).
    pub fn integrate_velocity(&mut self, dt: Real) {
        if self.body_type == BodyType::Static || self.state == BodyState::Sleeping {
            return;
        }
        if self.body_type == BodyType::Kinematic {
            if let Some(target) = &self.kinematic_target {
                let target_pos = target.target_position;
                let target_rot = target.target_rotation;
                self.velocity = (target_pos - self.transform.position) / dt.max(1e-10);
                self.transform.position = target_pos;
                self.transform.rotation = target_rot;
            } else {
                self.transform.position += self.velocity * dt;
                let omega = self.angular_velocity;
                let omega_len = omega.norm();
                if omega_len > 1e-10 {
                    let half_angle = omega_len * dt * 0.5;
                    let axis = omega / omega_len;
                    let delta_rot =
                        Quat::from_axis_angle(&Unit::new_normalize(axis), 2.0 * half_angle);
                    self.transform.rotation = delta_rot * self.transform.rotation;
                }
            }
            self.step_count += 1;
            return;
        }
        self.transform.position += self.velocity * dt;
        let omega = self.angular_velocity;
        let omega_len = omega.norm();
        if omega_len > 1e-10 {
            let half_angle = omega_len * dt * 0.5;
            let axis = omega / omega_len;
            let delta_rot = Quat::from_axis_angle(&Unit::new_normalize(axis), 2.0 * half_angle);
            self.transform.rotation = delta_rot * self.transform.rotation;
        }
        self.update_world_inertia();
        self.step_count += 1;
    }
}
