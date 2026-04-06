//! # CsvWriterConfig - Trait Implementations
//!
//! This module contains trait implementations for `CsvWriterConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::CsvWriterConfig;

impl Default for CsvWriterConfig {
    fn default() -> Self {
        Self {
            delimiter: ',',
            line_ending: "\n".to_string(),
            quote_all: false,
            precision: 6,
        }
    }
}
