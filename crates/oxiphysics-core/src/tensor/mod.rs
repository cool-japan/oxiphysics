// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tensor algebra operations for continuum mechanics.
//!
//! Provides second-order tensors (`Tensor2`), fourth-order tensors (`Tensor4`),
//! and a Voigt notation helper for symmetric second-order tensors.

mod decomposition;
mod operations;
mod types;

// Re-export all public items
pub use decomposition::*;
pub use operations::*;
pub use types::*;
