//! # MotorJoint - Trait Implementations
//!
//! This module contains trait implementations for `MotorJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_rigid::RigidBodySet;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{apply_angular_impulse, read_body};
use super::types::MotorJoint;

impl Constraint for MotorJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, _dt: f64) {
        if let Some(body) = bodies.get(self.body_a) {
            self.world_axis = (body.transform.rotation * self.local_axis).normalize();
        }
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let omega_a = a.ang_vel.dot(&self.world_axis);
        let omega_b = b.ang_vel.dot(&self.world_axis);
        let current_velocity = omega_a - omega_b;
        let velocity_error = self.target_velocity - current_velocity;
        let k_a = self.world_axis.dot(&(a.inv_inertia * self.world_axis));
        let k_b = self.world_axis.dot(&(b.inv_inertia * self.world_axis));
        let k = k_a + k_b;
        if k < 1e-12 {
            return;
        }
        let max_impulse = self.max_torque * dt;
        let lambda = (velocity_error / k).clamp(-max_impulse, max_impulse);
        let impulse = self.world_axis * lambda;
        apply_angular_impulse(bodies, self.body_a, self.body_b, impulse);
    }
    fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
