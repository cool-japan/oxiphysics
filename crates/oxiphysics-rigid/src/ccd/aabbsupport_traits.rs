//! # AabbSupport - Trait Implementations
//!
//! This module contains trait implementations for `AabbSupport`.
//!
//! ## Implemented Traits
//!
//! - `ConvexSupport`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::ConvexSupport;
#[allow(unused_imports)]
use super::functions::*;
use super::types::AabbSupport;

impl ConvexSupport for AabbSupport {
    fn support(&self, d: [f64; 3]) -> [f64; 3] {
        [
            if d[0] >= 0.0 {
                self.max[0]
            } else {
                self.min[0]
            },
            if d[1] >= 0.0 {
                self.max[1]
            } else {
                self.min[1]
            },
            if d[2] >= 0.0 {
                self.max[2]
            } else {
                self.min[2]
            },
        ]
    }
}
