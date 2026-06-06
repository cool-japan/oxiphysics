//! # Interval - Trait Implementations
//!
//! This module contains trait implementations for `Interval`.
//!
//! ## Implemented Traits
//!
//! - `Neg`
//! - `Add`
//! - `Sub`
//! - `Mul`
//! - `Div`
//! - `Mul`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::ops::{Add, Div, Mul, Neg, Sub};

use super::types::Interval;

impl Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
}

impl Add for Interval {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            lo: self.lo + rhs.lo,
            hi: self.hi + rhs.hi,
        }
    }
}

impl Sub for Interval {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            lo: self.lo - rhs.hi,
            hi: self.hi - rhs.lo,
        }
    }
}

impl Mul for Interval {
    type Output = Self;
    /// Interval multiplication using the four-corners rule.
    fn mul(self, rhs: Self) -> Self {
        let c1 = self.lo * rhs.lo;
        let c2 = self.lo * rhs.hi;
        let c3 = self.hi * rhs.lo;
        let c4 = self.hi * rhs.hi;
        Self {
            lo: c1.min(c2).min(c3).min(c4),
            hi: c1.max(c2).max(c3).max(c4),
        }
    }
}

impl Div for Interval {
    type Output = Self;
    /// Interval division. Panics if the divisor contains zero.
    ///
    /// Computes `[a,b] / [c,d]` directly via the four-corners rule on
    /// `[a,b] * [1/d, 1/c]`, equivalent to multiplication by the reciprocal.
    fn div(self, rhs: Self) -> Self {
        let rec = rhs.reciprocal();
        let c1 = self.lo * rec.lo;
        let c2 = self.lo * rec.hi;
        let c3 = self.hi * rec.lo;
        let c4 = self.hi * rec.hi;
        Self {
            lo: c1.min(c2).min(c3).min(c4),
            hi: c1.max(c2).max(c3).max(c4),
        }
    }
}

/// Multiply an interval by a scalar.
impl Mul<f64> for Interval {
    type Output = Interval;
    fn mul(self, s: f64) -> Interval {
        if s >= 0.0 {
            Interval {
                lo: self.lo * s,
                hi: self.hi * s,
            }
        } else {
            Interval {
                lo: self.hi * s,
                hi: self.lo * s,
            }
        }
    }
}
