//! # AnifHeader - Trait Implementations
//!
//! This module contains trait implementations for `AnifHeader`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{AnifEndian, AnifHeader};

impl Default for AnifHeader {
    fn default() -> Self {
        Self {
            version: (1, 0),
            num_points: 0,
            x_start: 0.0,
            x_delta: 1.0,
            endian: AnifEndian::Little,
            title: String::new(),
            instrument: String::new(),
            date: String::new(),
        }
    }
}
