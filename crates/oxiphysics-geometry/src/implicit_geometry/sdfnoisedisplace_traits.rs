//! # SdfNoiseDisplace - Trait Implementations
//!
//! This module contains trait implementations for `SdfNoiseDisplace`.
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
use super::types::SdfNoiseDisplace;

impl<S: Sdf> Sdf for SdfNoiseDisplace<S> {
    fn dist(&self, p: [f64; 3]) -> f64 {
        let n = fbm_noise(scale(p, self.noise_scale), self.octaves, 2.0, 0.5);
        self.inner.dist(p) + n * self.amplitude
    }
}
