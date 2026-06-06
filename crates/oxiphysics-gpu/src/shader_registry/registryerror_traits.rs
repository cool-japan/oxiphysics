//! # RegistryError - Trait Implementations
//!
//! This module contains trait implementations for `RegistryError`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//! - `Error`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::RegistryError;

impl std::fmt::Display for RegistryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RegistryError::UnknownShader(n) => write!(f, "unknown shader: {n}"),
            RegistryError::MissingDefine { shader, define } => {
                write!(f, "shader '{shader}' requires define '{define}'")
            }
            RegistryError::SourceTooLarge {
                shader,
                size,
                limit,
            } => {
                write!(f, "shader '{shader}' is {size} bytes (limit {limit})")
            }
        }
    }
}

impl std::error::Error for RegistryError {}
