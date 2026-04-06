// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core math types and operations for the OxiPhysics engine.
//!
//! **All other OxiPhysics crates should import math types from
//! `oxiphysics_core::math` rather than depending on nalgebra directly.**

mod dual_traits;
mod geometry;
mod linear_algebra;
mod types;

// Re-export all public items
pub use geometry::*;
pub use linear_algebra::*;
pub use types::*;
