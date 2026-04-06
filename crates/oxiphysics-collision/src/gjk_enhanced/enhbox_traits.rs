//! # EnhBox - Trait Implementations
//!
//! This module contains trait implementations for `EnhBox`.
//!
//! ## Implemented Traits
//!
//! - `ConvexShape`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ConvexShape;
#[allow(unused_imports)]
use super::functions::*;
use super::types::EnhBox;

impl ConvexShape for EnhBox {
    fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        [
            self.centre[0]
                + if dir[0] >= 0.0 {
                    self.half[0]
                } else {
                    -self.half[0]
                },
            self.centre[1]
                + if dir[1] >= 0.0 {
                    self.half[1]
                } else {
                    -self.half[1]
                },
            self.centre[2]
                + if dir[2] >= 0.0 {
                    self.half[2]
                } else {
                    -self.half[2]
                },
        ]
    }
    fn centre(&self) -> [f64; 3] {
        self.centre
    }
}
