//! # SdfLineSegment - Trait Implementations
//!
//! This module contains trait implementations for `SdfLineSegment`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfLineSegment;

impl Sdf for SdfLineSegment {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let pa = sub(p, self.a);
        let ba = sub(self.b, self.a);
        let h = clamp01(dot(pa, ba) / dot(ba, ba));
        len(sub(pa, scale(ba, h))) - self.radius
    }
}
