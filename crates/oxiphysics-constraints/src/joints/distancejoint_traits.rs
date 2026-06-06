//! # DistanceJoint - Trait Implementations
//!
//! This module contains trait implementations for `DistanceJoint`.
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

use super::functions::{JOINT_BAUMGARTE, apply_pair_impulse, linear_effective_mass, read_body};
use super::types::DistanceJoint;

impl Constraint for DistanceJoint {
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
        let world_a = a.position + self.r_a;
        let world_b = b.position + self.r_b;
        let diff = world_b - world_a;
        let dist = diff.norm();
        self.constraint_dir = if dist > 1e-10 {
            diff / dist
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        self.eff_mass = linear_effective_mass(
            a.inv_mass,
            b.inv_mass,
            &a.inv_inertia,
            &b.inv_inertia,
            &self.r_a,
            &self.r_b,
            &self.constraint_dir,
        );
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
        let rel_v = (vb - va).dot(&self.constraint_dir);
        let world_a = a.position + self.r_a;
        let world_b = b.position + self.r_b;
        let dist = (world_b - world_a).norm();
        let active = dist < self.min_distance || dist > self.max_distance;
        let target_v = 0.0_f64;
        if active && self.eff_mass > 1e-12 {
            let lambda = (target_v - rel_v) * self.eff_mass;
            let impulse = self.constraint_dir * lambda;
            apply_pair_impulse(
                bodies,
                self.body_a,
                self.body_b,
                impulse,
                &self.r_a,
                &self.r_b,
            );
            self.accumulated_lambda += lambda;
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
        let diff = world_b - world_a;
        let dist = diff.norm();
        let error;
        if dist < self.min_distance {
            error = dist - self.min_distance;
        } else if dist > self.max_distance {
            error = dist - self.max_distance;
        } else {
            return;
        }
        let dir = if dist > 1e-10 {
            diff / dist
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        let correction = dir * (error * JOINT_BAUMGARTE);
        let k = a.inv_mass + b.inv_mass;
        if k > 1e-12 {
            if let Some(body) = bodies.get_mut(self.body_a) {
                body.transform.position += correction * (a.inv_mass / k);
            }
            if let Some(body) = bodies.get_mut(self.body_b) {
                body.transform.position -= correction * (b.inv_mass / k);
            }
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
