//! # SdfRepeat - Trait Implementations
//!
//! This module contains trait implementations for `SdfRepeat`.
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
use super::types::SdfRepeat;

impl<S: Sdf> Sdf for SdfRepeat<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        fn rmod(v: f64, period: f64) -> f64 {
            if period < 1e-15 {
                return v;
            }
            v - period * (v / period).round()
        }
        let q = [
            rmod(p[0], self.period_x),
            rmod(p[1], self.period_y),
            rmod(p[2], self.period_z),
        ];
        self.inner.dist(q)
    }
}
