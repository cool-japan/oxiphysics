//! # SdfGyroid - Trait Implementations
//!
//! This module contains trait implementations for `SdfGyroid`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::*;
use super::types::SdfGyroid;

impl Sdf for SdfGyroid {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let q = scale(p, self.scale * 2.0 * PI);
        let v = q[0].sin() * q[1].cos() + q[1].sin() * q[2].cos() + q[2].sin() * q[0].cos();
        v.abs() / self.scale - self.thickness
    }
}
