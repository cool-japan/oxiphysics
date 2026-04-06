//! # RigidBody - new_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::Transform;
use oxiphysics_core::math::{Mat3, Quat, Real, Unit, Vec3};

use super::rigidbody_type::RigidBody;
use super::types::{BodySnapshot, BodyState, BodyType};

impl RigidBody {
    /// Create a new dynamic rigid body with the given mass.
    pub fn new(mass: Real) -> Self {
        let inverse_mass = if mass > 0.0 { 1.0 / mass } else { 0.0 };
        Self {
            transform: Transform::default(),
            velocity: Vec3::zeros(),
            angular_velocity: Vec3::zeros(),
            mass,
            inverse_mass,
            local_inertia: Mat3::identity() * mass,
            world_inverse_inertia: Mat3::identity() * inverse_mass,
            force_accumulator: Vec3::zeros(),
            torque_accumulator: Vec3::zeros(),
            body_type: BodyType::Dynamic,
            state: BodyState::Active,
            linear_damping: 0.01,
            angular_damping: 0.01,
            gravity_scale: 1.0,
            sleep_timer: 0.0,
            kinematic_target: None,
            step_count: 0,
        }
    }
    /// Create a new static rigid body.
    pub fn new_static() -> Self {
        let mut body = Self::new(0.0);
        body.body_type = BodyType::Static;
        body.inverse_mass = 0.0;
        body.world_inverse_inertia = Mat3::zeros();
        body.local_inertia = Mat3::zeros();
        body
    }
    /// Create a new kinematic rigid body.
    pub fn new_kinematic() -> Self {
        let mut body = Self::new(0.0);
        body.body_type = BodyType::Kinematic;
        body.inverse_mass = 0.0;
        body.world_inverse_inertia = Mat3::zeros();
        body.local_inertia = Mat3::zeros();
        body
    }
    /// Restore body state from a snapshot.
    #[allow(dead_code)]
    pub fn restore_from_snapshot(&mut self, snap: &BodySnapshot) {
        self.transform.position = Vec3::new(snap.position[0], snap.position[1], snap.position[2]);
        let qw = snap.rotation[3];
        let qx = snap.rotation[0];
        let qy = snap.rotation[1];
        let qz = snap.rotation[2];
        let axis_len = (qx * qx + qy * qy + qz * qz).sqrt();
        if axis_len > 1e-12 {
            let angle = 2.0 * axis_len.atan2(qw);
            let axis = Unit::new_normalize(Vec3::new(qx, qy, qz));
            self.transform.rotation = Quat::from_axis_angle(&axis, angle);
        } else {
            self.transform.rotation = Quat::identity();
        }
        self.velocity = Vec3::new(snap.velocity[0], snap.velocity[1], snap.velocity[2]);
        self.angular_velocity = Vec3::new(
            snap.angular_velocity[0],
            snap.angular_velocity[1],
            snap.angular_velocity[2],
        );
        self.mass = snap.mass;
        self.inverse_mass = if snap.mass > 0.0 {
            1.0 / snap.mass
        } else {
            0.0
        };
        self.linear_damping = snap.linear_damping;
        self.angular_damping = snap.angular_damping;
        self.gravity_scale = snap.gravity_scale;
        self.body_type = match snap.body_type {
            "kinematic" => BodyType::Kinematic,
            "static" => BodyType::Static,
            _ => BodyType::Dynamic,
        };
        self.state = match snap.state {
            "sleeping" => BodyState::Sleeping,
            _ => BodyState::Active,
        };
        self.update_world_inertia();
    }
    /// Merge two bodies into one combined body (at center of mass).
    ///
    /// The resulting body has combined mass, weighted center of mass,
    /// and momentum-conserving velocity.
    #[allow(dead_code)]
    pub fn merge(a: &RigidBody, b: &RigidBody) -> RigidBody {
        let total_mass = a.mass + b.mass;
        if total_mass < 1e-30 {
            return RigidBody::new_static();
        }
        let inv_total = 1.0 / total_mass;
        let com = (a.transform.position * a.mass + b.transform.position * b.mass) * inv_total;
        let combined_velocity = (a.velocity * a.mass + b.velocity * b.mass) * inv_total;
        let combined_angular_momentum =
            a.local_inertia * a.angular_velocity + b.local_inertia * b.angular_velocity;
        let combined_inertia = a.local_inertia + b.local_inertia;
        let mut merged = RigidBody::new(total_mass);
        merged.transform.position = com;
        merged.velocity = combined_velocity;
        merged.local_inertia = combined_inertia;
        if let Some(inv_inertia) = combined_inertia.try_inverse() {
            merged.angular_velocity = inv_inertia * combined_angular_momentum;
        }
        merged.update_world_inertia();
        merged
    }
}
