//! # SdfExtrude - Trait Implementations
//!
//! This module contains trait implementations for `SdfExtrude`.
//!
//! ## Implemented Traits
//!
//! - `Sdf`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use super::types::SdfExtrude;

impl<F: Fn([f64; 2]) -> f64 + Send + Sync> Sdf for SdfExtrude<F> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let d2 = (self.profile)([p[0], p[2]]);
        let dy = p[1].abs() - self.half_height;
        let out_lateral = d2.max(0.0);
        let out_cap = dy.max(0.0);
        let inside = d2.min(0.0).max(dy.min(0.0));
        len2([out_lateral, out_cap]) + inside
    }
}
