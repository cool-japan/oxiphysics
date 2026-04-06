//! # Polynomial - Trait Implementations
//!
//! This module contains trait implementations for `Polynomial`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::Polynomial;
use std::fmt;

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for (i, &c) in self.coeffs.iter().enumerate().rev() {
            if c.abs() < 1e-15 {
                continue;
            }
            if !first && c > 0.0 {
                write!(f, " + ")?;
            } else if !first && c < 0.0 {
                write!(f, " - ")?;
            }
            let ac = if first { c } else { c.abs() };
            first = false;
            if i == 0 {
                write!(f, "{:.6}", ac)?;
            } else if i == 1 {
                if (ac - 1.0).abs() < 1e-15 {
                    write!(f, "x")?;
                } else {
                    write!(f, "{:.6}*x", ac)?;
                }
            } else if (ac - 1.0).abs() < 1e-15 {
                write!(f, "x^{}", i)?;
            } else {
                write!(f, "{:.6}*x^{}", ac, i)?;
            }
        }
        if first {
            write!(f, "0")?;
        }
        Ok(())
    }
}
