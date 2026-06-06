//! # VoxelVolume - Trait Implementations
//!
//! This module contains trait implementations for `VoxelVolume`.
//!
//! ## Implemented Traits
//!
//! - `Display`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::VoxelVolume;
use std::fmt;

impl fmt::Display for VoxelVolume {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "VoxelVolume({}x{}x{}, spacing=[{:.3},{:.3},{:.3}])",
            self.dims[0],
            self.dims[1],
            self.dims[2],
            self.spacing[0],
            self.spacing[1],
            self.spacing[2],
        )
    }
}
