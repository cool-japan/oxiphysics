//! # SphericalJoint - Trait Implementations
//!
//! This module contains trait implementations for `SphericalJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_rigid::RigidBodySet;

use super::functions::JOINT_BAUMGARTE;
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{apply_pair_impulse, read_body};
use super::types::SphericalJoint;

impl Constraint for SphericalJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, _dt: f64) {
        if let Some(body) = bodies.get(self.body_a) {
            self.r_a = body.transform.rotation * self.local_anchor_a;
        }
        if let Some(body) = bodies.get(self.body_b) {
            self.r_b = body.transform.rotation * self.local_anchor_b;
        }
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
        let identity = oxiphysics_core::math::Mat3::identity();
        let k = (a.inv_mass + b.inv_mass) * identity + a.inv_inertia + b.inv_inertia;
        if let Some(k_inv) = k.try_inverse() {
            let impulse = k_inv * (-dv);
            apply_pair_impulse(
                bodies,
                self.body_a,
                self.body_b,
                impulse,
                &self.r_a,
                &self.r_b,
            );
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
        if k < 1e-12 {
            return;
        }
        let correction = error * JOINT_BAUMGARTE;
        if let Some(body) = bodies.get_mut(self.body_a) {
            body.transform.position -= correction * (a.inv_mass / k);
        }
        if let Some(body) = bodies.get_mut(self.body_b) {
            body.transform.position += correction * (b.inv_mass / k);
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
