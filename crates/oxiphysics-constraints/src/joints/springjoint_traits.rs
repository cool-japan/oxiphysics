//! # SpringJoint - Trait Implementations
//!
//! This module contains trait implementations for `SpringJoint`.
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

use super::functions::{apply_pair_impulse, read_body};
use super::types::SpringJoint;

impl Constraint for SpringJoint {
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
        let delta = world_b - world_a;
        self.current_distance = delta.norm();
        if self.current_distance > 1e-10 {
            self.direction = delta / self.current_distance;
        } else {
            self.direction = Vec3::new(1.0, 0.0, 0.0);
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
        let stretch = self.current_distance - self.rest_length;
        let spring_force = self.stiffness * stretch;
        let va = a.vel + a.ang_vel.cross(&self.r_a);
        let vb = b.vel + b.ang_vel.cross(&self.r_b);
        let rel_v = (vb - va).dot(&self.direction);
        let damping_force = self.damping * rel_v;
        let total_force = spring_force + damping_force;
        let impulse = self.direction * (total_force * dt);
        apply_pair_impulse(
            bodies,
            self.body_a,
            self.body_b,
            impulse,
            &self.r_a,
            &self.r_b,
        );
    }
    fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
