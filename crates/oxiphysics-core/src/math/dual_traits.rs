//! # Dual - Trait Implementations
//!
//! This module contains trait implementations for `Dual`.
//!
//! ## Implemented Traits
//!
//! - `Add`
//! - `Sub`
//! - `Mul`
//! - `Div`
//! - `Neg`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::Dual;

impl std::ops::Add for Dual {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            du: self.du + rhs.du,
        }
    }
}

impl std::ops::Sub for Dual {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            du: self.du - rhs.du,
        }
    }
}

impl std::ops::Mul for Dual {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re,
            du: self.re * rhs.du + self.du * rhs.re,
        }
    }
}

impl std::ops::Div for Dual {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let re = self.re / rhs.re;
        let du = (self.du * rhs.re - self.re * rhs.du) / (rhs.re * rhs.re);
        Self { re, du }
    }
}

impl std::ops::Neg for Dual {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            du: -self.du,
        }
    }
}
