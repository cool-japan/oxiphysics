//! # IncrementalExportConfig - Trait Implementations
//!
//! This module contains trait implementations for `IncrementalExportConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::IncrementalExportConfig;

impl Default for IncrementalExportConfig {
    fn default() -> Self {
        Self {
            min_speed_threshold: 0.0,
            max_batch_size: 1000,
            include_sleeping: true,
        }
    }
}
