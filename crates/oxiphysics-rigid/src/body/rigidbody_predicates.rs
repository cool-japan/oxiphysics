//! # RigidBody - predicates Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::{Real, Vec3};

use super::types::{BodyState, BodyType};

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Is this body dynamic?
    pub fn is_dynamic(&self) -> bool {
        self.body_type == BodyType::Dynamic
    }
    /// Apply a force at the center of mass (accumulated until next step).
    pub fn apply_force(&mut self, force: Vec3) {
        if self.is_dynamic() {
            self.force_accumulator += force;
        }
    }
    /// Apply a torque (accumulated until next step).
    pub fn apply_torque(&mut self, torque: Vec3) {
        if self.is_dynamic() {
            self.torque_accumulator += torque;
        }
    }
    /// Apply an impulse at the center of mass (instantaneous velocity change).
    pub fn apply_impulse(&mut self, impulse: Vec3) {
        if self.is_dynamic() {
            self.velocity += impulse * self.inverse_mass;
        }
    }
    /// Apply an impulse at a world-space point.
    pub fn apply_impulse_at_point(&mut self, impulse: Vec3, point: Vec3) {
        if self.is_dynamic() {
            self.velocity += impulse * self.inverse_mass;
            let r = point - self.transform.position;
            self.angular_velocity += self.world_inverse_inertia * r.cross(&impulse);
        }
    }
    /// Apply a force at a world-space point.
    pub fn apply_force_at_point(&mut self, force: Vec3, point: Vec3) {
        if self.is_dynamic() {
            self.force_accumulator += force;
            let r = point - self.transform.position;
            self.torque_accumulator += r.cross(&force);
        }
    }
    /// Integrate forces -> velocity (semi-implicit Euler, first half).
    pub fn integrate_forces(&mut self, dt: Real, gravity: &Vec3) {
        if !self.is_dynamic() || self.state == BodyState::Sleeping {
            return;
        }
        self.velocity += *gravity * self.gravity_scale * dt;
        self.velocity += self.force_accumulator * self.inverse_mass * dt;
        self.angular_velocity += self.world_inverse_inertia * self.torque_accumulator * dt;
        self.velocity *= (1.0 - self.linear_damping).max(0.0);
        self.angular_velocity *= (1.0 - self.angular_damping).max(0.0);
        self.force_accumulator = Vec3::zeros();
        self.torque_accumulator = Vec3::zeros();
    }
    /// Check if this body should go to sleep.
    pub fn check_sleep(
        &mut self,
        dt: Real,
        linear_threshold: Real,
        angular_threshold: Real,
        time_before_sleep: Real,
    ) {
        if !self.is_dynamic() {
            return;
        }
        if self.velocity.norm_squared() < linear_threshold * linear_threshold
            && self.angular_velocity.norm_squared() < angular_threshold * angular_threshold
        {
            self.sleep_timer += dt;
            if self.sleep_timer >= time_before_sleep {
                self.state = BodyState::Sleeping;
                self.velocity = Vec3::zeros();
                self.angular_velocity = Vec3::zeros();
            }
        } else {
            self.sleep_timer = 0.0;
            self.state = BodyState::Active;
        }
    }
    /// Apply angular impulse (instantaneous angular velocity change).
    #[allow(dead_code)]
    pub fn apply_angular_impulse(&mut self, impulse: Vec3) {
        if self.is_dynamic() {
            self.angular_velocity += self.world_inverse_inertia * impulse;
        }
    }
}
