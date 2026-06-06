//! # Joint - Trait Implementations
//!
//! This module contains trait implementations for `Joint`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::Joint;

impl Default for Joint {
    fn default() -> Self {
        let mut ibm = [0.0f32; 16];
        ibm[0] = 1.0;
        ibm[5] = 1.0;
        ibm[10] = 1.0;
        ibm[15] = 1.0;
        Joint {
            name: String::new(),
            translation: [0.0, 0.0, 0.0],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0, 1.0, 1.0],
            children: Vec::new(),
            inverse_bind_matrix: ibm,
        }
    }
}
