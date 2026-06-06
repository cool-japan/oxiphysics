// # WheelContact - Trait Implementations
//
// This module contains trait implementations for `WheelContact`.
//
// ## Implemented Traits
//
// - `Default`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::WheelContact;

impl Default for WheelContact {
    fn default() -> Self {
        Self {
            hit: false,
            position: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
            penetration: 0.0,
            body_idx: None,
        }
    }
}
