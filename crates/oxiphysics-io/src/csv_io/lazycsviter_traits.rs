//! # LazyCsvIter - Trait Implementations
//!
//! This module contains trait implementations for `LazyCsvIter`.
//!
//! ## Implemented Traits
//!
//! - `Iterator`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::LazyCsvIter;

impl<'a> Iterator for LazyCsvIter<'a> {
    type Item = Vec<String>;
    fn next(&mut self) -> Option<Self::Item> {
        let line = self.lines.next()?;
        Some(
            line.split(self.delimiter)
                .map(str::trim)
                .map(String::from)
                .collect(),
        )
    }
}
