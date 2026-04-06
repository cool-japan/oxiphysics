//! # SdfBox - Trait Implementations
//!
//! This module contains trait implementations for `SdfBox`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions::*;
use super::types::SdfBox;

impl Sdf for SdfBox {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let q = [
            p[0].abs() - self.half_extents[0],
            p[1].abs() - self.half_extents[1],
            p[2].abs() - self.half_extents[2],
        ];
        len([q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)]) + q[0].max(q[1]).max(q[2]).min(0.0)
    }
}
