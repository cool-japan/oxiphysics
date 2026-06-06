//! # QuatAlgebra - Trait Implementations
//!
//! This module contains trait implementations for `QuatAlgebra`.
//!
//! ## Implemented Traits
//!
//! - `Mul`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::ops::Mul;

use super::types::QuatAlgebra;

impl Mul for QuatAlgebra {
    type Output = Self;
    /// Hamilton product.
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
            self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
            self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
            self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w,
        )
    }
}
