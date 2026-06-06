//! # MolVizError - Trait Implementations
//!
//! This module contains trait implementations for `MolVizError`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//! - `Error`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::MolVizError;
use std::fmt;

impl fmt::Display for MolVizError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ParseError { message } => write!(f, "parse error: {message}"),
            Self::MissingSection { section } => write!(f, "missing section: {section}"),
            Self::DimensionMismatch { message } => {
                write!(f, "dimension mismatch: {message}")
            }
            Self::UnknownFormat { detected } => write!(f, "unknown format: {detected}"),
            Self::IoError { message } => write!(f, "I/O error: {message}"),
            Self::OutOfBounds { message } => write!(f, "out of bounds: {message}"),
        }
    }
}

impl std::error::Error for MolVizError {}
