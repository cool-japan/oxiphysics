// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Computational homogenization for multiscale materials.
//!
//! Implements RVE-based homogenization and asymptotic expansion methods
//! including Voigt-Reuss-Hill bounds, Mori-Tanaka micromechanics,
//! Hill-Mandel condition verification, effective property extraction,
//! periodic RVE analysis, and multi-scale FE coupling.

pub mod advanced;
pub mod bounds;
pub mod cell;
pub mod extended;
pub mod matrix_utils;
pub mod schemes;

// Re-export everything for backwards compatibility
pub use advanced::*;
pub use bounds::*;
pub use cell::*;
pub use extended::*;
pub use matrix_utils::*;
pub use schemes::*;
