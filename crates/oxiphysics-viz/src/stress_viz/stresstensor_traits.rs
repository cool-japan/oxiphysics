//! # StressTensor - Trait Implementations
//!
//! This module contains trait implementations for `StressTensor`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Add`
//! - `Sub`
//! - `Mul`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::StressTensor;

impl Default for StressTensor {
    fn default() -> Self {
        Self::zero()
    }
}

impl std::ops::Add for StressTensor {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut v = [0.0f64; 6];
        for (v_i, (a, b)) in v.iter_mut().zip(self.voigt.iter().zip(rhs.voigt.iter())) {
            *v_i = a + b;
        }
        Self::new(v)
    }
}

impl std::ops::Sub for StressTensor {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut v = [0.0f64; 6];
        for (v_i, (a, b)) in v.iter_mut().zip(self.voigt.iter().zip(rhs.voigt.iter())) {
            *v_i = a - b;
        }
        Self::new(v)
    }
}

impl std::ops::Mul<f64> for StressTensor {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        let mut v = self.voigt;
        for x in v.iter_mut() {
            *x *= s;
        }
        Self::new(v)
    }
}
