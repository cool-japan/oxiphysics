//! # BreakableJoint - Trait Implementations
//!
//! This module contains trait implementations for `BreakableJoint`.
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

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{apply_pair_impulse, read_body};
use super::types::BreakableJoint;

impl Constraint for BreakableJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, _dt: f64) {
        if self.broken {
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
        self.r_a = a.rotation * Vec3::zeros();
        self.r_b = b.rotation * Vec3::zeros();
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        if self.broken {
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
        let world_a = a.position;
        let world_b = b.position;
        let diff = world_b - world_a;
        let dist = diff.norm();
        let dir = if dist > 1e-10 { diff / dist } else { return };
        let extension = dist - self.natural_length;
        let rel_v = (b.vel - a.vel).dot(&dir);
        let spring_force = self.stiffness * extension + self.damping * rel_v;
        let lambda = spring_force * dt;
        self.accumulated_force += lambda.abs();
        if self.accumulated_force > self.break_threshold {
            self.broken = true;
            return;
        }
        let impulse = dir * lambda;
        let k = a.inv_mass + b.inv_mass;
        if k > 1e-12 {
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
    fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
