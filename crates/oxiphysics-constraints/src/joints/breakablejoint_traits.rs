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

use super::constraint_math::{SoftRow, spring_softness};
use super::functions::read_body;
use super::types::BreakableJoint;

impl BreakableJoint {
    /// One (possibly relaxed) Gauss-Seidel update of the spring row, followed
    /// by the break test on the impulse transmitted so far in this step.
    fn iterate_spring(&mut self, bodies: &mut RigidBodySet, omega: f64) -> f64 {
        if self.broken {
            return 0.0;
        }
        let update = self.spring_row.iterate(
            bodies,
            self.body_a,
            self.body_b,
            &self.r_a,
            &self.r_b,
            omega,
        );
        self.accumulated_force = self.spring_row.lambda().abs();
        if self.accumulated_force > self.break_threshold {
            self.broken = true;
        }
        update
    }
}

impl Constraint for BreakableJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
        self.accumulated_force = 0.0;
        self.spring_row = SoftRow::inactive();
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
        self.r_a = Vec3::zeros();
        self.r_b = Vec3::zeros();
        let diff = b.position - a.position;
        let dist = diff.norm();
        if dist <= 1e-10 {
            return;
        }
        // Implicit spring-damper between the body centres: the row's
        // accumulated impulse is the impulse transmitted during this step.
        let stretch = dist - self.natural_length;
        if let Some((gamma, bias_scale)) = spring_softness(self.stiffness, self.damping, dt) {
            self.spring_row = SoftRow::new(
                &a,
                &b,
                &self.r_a,
                &self.r_b,
                diff / dist,
                gamma,
                stretch * bias_scale,
            );
        }
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        self.iterate_spring(bodies, 1.0);
    }
    fn solve_velocity_relaxed(&mut self, bodies: &mut RigidBodySet, _dt: f64, omega: f64) -> f64 {
        self.iterate_spring(bodies, omega)
    }
    fn residual(&self) -> f64 {
        self.spring_row.last_update()
    }
    fn is_active(&self) -> bool {
        !self.broken
    }
    fn solve_position(&mut self, _bodies: &mut RigidBodySet, _dt: f64) {}
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
