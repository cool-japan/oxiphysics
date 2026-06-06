//! # GltfNode - Trait Implementations
//!
//! This module contains trait implementations for `GltfNode`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::GltfNode;

impl Default for GltfNode {
    fn default() -> Self {
        GltfNode {
            name: String::new(),
            mesh: None,
            translation: [0.0, 0.0, 0.0],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0, 1.0, 1.0],
            children: Vec::new(),
        }
    }
}
