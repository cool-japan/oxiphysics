//! # UniversalJoint - Trait Implementations
//!
//! This module contains trait implementations for `UniversalJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_rigid::RigidBodySet;

use super::functions::{
    JOINT_BAUMGARTE, apply_angular_impulse, apply_pair_impulse, read_body, small_rotation_quat,
};
use super::types::UniversalJoint;

impl Constraint for UniversalJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, _dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        self.r_a = a.rotation * self.local_anchor_a;
        self.r_b = b.rotation * self.local_anchor_b;
        self.world_axis_a = (a.rotation * self.local_axis_a).normalize();
        self.world_axis_b = (b.rotation * self.local_axis_b).normalize();
        self.constrained_axis = self.world_axis_a.cross(&self.world_axis_b).normalize();
        let _ = b;
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let va = a.vel + a.ang_vel.cross(&self.r_a);
        let vb = b.vel + b.ang_vel.cross(&self.r_b);
        let dv = va - vb;
        let k = a.inv_mass + b.inv_mass;
        if k > 1e-12 {
            let impulse = -dv / k;
            apply_pair_impulse(
                bodies,
                self.body_a,
                self.body_b,
                impulse,
                &self.r_a,
                &self.r_b,
            );
        }
        let a2 = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b2 = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let d_omega = a2.ang_vel - b2.ang_vel;
        let constrained_component = d_omega.dot(&self.constrained_axis) * self.constrained_axis;
        if constrained_component.norm_squared() > 1e-20 {
            let ang_k = a2.inv_inertia + b2.inv_inertia;
            if let Some(ang_k_inv) = ang_k.try_inverse() {
                let ang_impulse = -(ang_k_inv * constrained_component);
                apply_angular_impulse(bodies, self.body_a, self.body_b, ang_impulse);
            }
        }
    }
    fn solve_position(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let world_a = a.position + a.rotation * self.local_anchor_a;
        let world_b = b.position + b.rotation * self.local_anchor_b;
        let error = world_a - world_b;
        let k = a.inv_mass + b.inv_mass;
        if k > 1e-12 {
            let correction = error * JOINT_BAUMGARTE;
            if let Some(body) = bodies.get_mut(self.body_a) {
                body.transform.position -= correction * (a.inv_mass / k);
            }
            if let Some(body) = bodies.get_mut(self.body_b) {
                body.transform.position += correction * (b.inv_mass / k);
            }
        }
        let angular_err = self.angular_error();
        if angular_err.abs() > 1e-6 {
            let ang_correction = self.constrained_axis * (angular_err * JOINT_BAUMGARTE);
            let ang_k = a.inv_inertia + b.inv_inertia;
            if let Some(ang_k_inv) = ang_k.try_inverse() {
                let delta = ang_k_inv * ang_correction;
                if let Some(body) = bodies.get_mut(self.body_a) {
                    let dq = small_rotation_quat(&(-delta * 0.5));
                    body.transform.rotation = dq * body.transform.rotation;
                }
                if let Some(body) = bodies.get_mut(self.body_b) {
                    let dq = small_rotation_quat(&(delta * 0.5));
                    body.transform.rotation = dq * body.transform.rotation;
                }
            }
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
