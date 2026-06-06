//! # IslandBodyIter - Trait Implementations
//!
//! This module contains trait implementations for `IslandBodyIter`.
//!
//! ## Implemented Traits
//!
//! - `Iterator`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::BodyHandle;

use super::types::IslandBodyIter;

impl<'a> Iterator for IslandBodyIter<'a> {
    type Item = (BodyHandle, f64, f64);
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let &handle = self.handles.next()?;
            if let Some(body) = self.bodies.get(handle) {
                let lin = body.velocity.norm();
                let ang = body.angular_velocity.norm();
                return Some((handle, lin, ang));
            }
        }
    }
}
