//! # GpuError - Trait Implementations
//!
//! This module contains trait implementations for `GpuError`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::GpuError;

impl std::fmt::Display for GpuError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GpuError::InvalidBuffer(id) => write!(f, "invalid buffer id: {}", id.0),
            GpuError::SizeMismatch { expected, got } => {
                write!(f, "size mismatch: expected {expected}, got {got}")
            }
            GpuError::EmptyBuffer => write!(f, "reduction on empty buffer"),
            GpuError::NotFound(name) => write!(f, "not found: {name}"),
        }
    }
}
