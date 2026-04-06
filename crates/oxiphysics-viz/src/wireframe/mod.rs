// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Debug wireframe generation for bounding volumes and contact visualization.

pub mod extras;
pub mod mesh_ops;
pub mod primitives;
pub mod renderer;

// Re-export everything for backwards compatibility
pub use extras::*;
pub use mesh_ops::*;
pub use primitives::*;
pub use renderer::*;
