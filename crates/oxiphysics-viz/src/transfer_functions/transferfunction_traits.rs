//! # TransferFunction - Trait Implementations
//!
//! This module contains trait implementations for `TransferFunction`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{Rgba, TransferFunction};

impl Default for TransferFunction {
    fn default() -> Self {
        let mut tf = Self { points: Vec::new() };
        tf.add_point(0.0, Rgba::transparent());
        tf.add_point(1.0, Rgba::white());
        tf
    }
}
