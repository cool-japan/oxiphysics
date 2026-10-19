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

use super::constraint_math::{SoftRow, spring_softness};
use super::functions::read_body;
use super::types::SpringJoint;

impl Constraint for SpringJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
        self.spring_row = SoftRow::inactive();
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
        // Implicit spring-damper: one soft row whose accumulated impulse is the
        // impulse of the whole step, however many iterations the solver runs.
        let stretch = self.current_distance - self.rest_length;
        if let Some((gamma, bias_scale)) = spring_softness(self.stiffness, self.damping, dt) {
            self.spring_row = SoftRow::new(
                &a,
                &b,
                &self.r_a,
                &self.r_b,
                self.direction,
                gamma,
                stretch * bias_scale,
            );
        }
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        self.spring_row
            .iterate(bodies, self.body_a, self.body_b, &self.r_a, &self.r_b, 1.0);
    }
    fn solve_velocity_relaxed(&mut self, bodies: &mut RigidBodySet, _dt: f64, omega: f64) -> f64 {
        self.spring_row.iterate(
            bodies,
            self.body_a,
            self.body_b,
            &self.r_a,
            &self.r_b,
            omega,
        )
    }
    fn residual(&self) -> f64 {
        self.spring_row.last_update()
    }
    fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
