//! # PrismaticJoint - Trait Implementations
//!
//! This module contains trait implementations for `PrismaticJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBodySet;

use super::functions::JOINT_BAUMGARTE;
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{apply_angular_impulse, apply_pair_impulse, read_body};
use super::types::PrismaticJoint;

impl Constraint for PrismaticJoint {
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
        self.world_axis = (a.rotation * self.local_axis).normalize();
        let reference = if self.world_axis.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        self.perp1 = self.world_axis.cross(&reference).normalize();
        self.perp2 = self.world_axis.cross(&self.perp1);
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
        if k < 1e-12 {
            return;
        }
        let free_v = dv.dot(&self.world_axis) * self.world_axis;
        let constrained_v = dv - free_v;
        if constrained_v.norm_squared() > 1e-20 {
            let impulse = -constrained_v / k;
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
        if d_omega.norm_squared() > 1e-20 {
            let ang_k = a2.inv_inertia + b2.inv_inertia;
            if let Some(ang_k_inv) = ang_k.try_inverse() {
                let ang_impulse = -(ang_k_inv * d_omega);
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
        let along = error.dot(&self.world_axis) * self.world_axis;
        let perp_error = error - along;
        let k = a.inv_mass + b.inv_mass;
        if k < 1e-12 {
            return;
        }
        let correction = perp_error * JOINT_BAUMGARTE;
        if let Some(body) = bodies.get_mut(self.body_a) {
            body.transform.position -= correction * (a.inv_mass / k);
        }
        if let Some(body) = bodies.get_mut(self.body_b) {
            body.transform.position += correction * (b.inv_mass / k);
        }
        if let (Some(lower), Some(upper)) = (self.lower_limit, self.upper_limit) {
            let dist = along.dot(&self.world_axis);
            if dist < lower {
                let limit_correction = self.world_axis * ((lower - dist) * JOINT_BAUMGARTE);
                if let Some(body) = bodies.get_mut(self.body_a) {
                    body.transform.position += limit_correction * (a.inv_mass / k);
                }
                if let Some(body) = bodies.get_mut(self.body_b) {
                    body.transform.position -= limit_correction * (b.inv_mass / k);
                }
            } else if dist > upper {
                let limit_correction = self.world_axis * ((upper - dist) * JOINT_BAUMGARTE);
                if let Some(body) = bodies.get_mut(self.body_a) {
                    body.transform.position += limit_correction * (a.inv_mass / k);
                }
                if let Some(body) = bodies.get_mut(self.body_b) {
                    body.transform.position -= limit_correction * (b.inv_mass / k);
                }
            }
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
