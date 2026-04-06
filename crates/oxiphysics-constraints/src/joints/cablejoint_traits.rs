//! # CableJoint - Trait Implementations
//!
//! This module contains trait implementations for `CableJoint`.
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
use super::functions::{linear_effective_mass, read_body};
use super::types::CableJoint;

impl Constraint for CableJoint {
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
        self.taut = dist >= self.max_length - 1e-6;
        self.cable_dir = if dist > 1e-10 {
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
            &self.cable_dir,
        );
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        self.apply_tension_constraint(bodies, dt);
    }
    fn solve_position(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        if !self.taut {
            return;
        }
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
        if dist <= self.max_length {
            return;
        }
        let error = dist - self.max_length;
        let dir = if dist > 1e-10 {
            diff / dist
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        let correction = dir * (error * JOINT_BAUMGARTE);
        let k = a.inv_mass + b.inv_mass;
        if k < 1e-12 {
            return;
        }
        if let Some(body) = bodies.get_mut(self.body_a) {
            body.transform.position += correction * (a.inv_mass / k);
        }
        if let Some(body) = bodies.get_mut(self.body_b) {
            body.transform.position -= correction * (b.inv_mass / k);
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
