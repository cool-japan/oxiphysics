//! # CollisionConstraint - Trait Implementations
//!
//! This module contains trait implementations for `CollisionConstraint`.
//!
//! ## Implemented Traits
//!
//! - `SoftConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::particle::SoftParticle;
use oxiphysics_core::math::Real;

use super::functions::SoftConstraint;
#[allow(unused_imports)]
use super::functions::*;
use super::types::CollisionConstraint;

impl SoftConstraint for CollisionConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], _dt_sub: Real) {
        let p = &mut particles[self.particle_index];
        if p.is_static() {
            return;
        }
        let d = self.normal.dot(&(p.position - self.point_on_plane));
        if d < 0.0 {
            p.position -= self.normal * d;
        }
    }
}
