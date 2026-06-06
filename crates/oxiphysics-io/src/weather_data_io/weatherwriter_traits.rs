//! # WeatherWriter - Trait Implementations
//!
//! This module contains trait implementations for `WeatherWriter`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::WeatherWriter;

impl Default for WeatherWriter {
    fn default() -> Self {
        Self {
            include_header: true,
            separator: ',',
            missing_value: "NaN".to_string(),
            precision: 2,
        }
    }
}
