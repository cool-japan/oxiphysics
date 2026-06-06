// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Corotational FEM soft body simulation.
//!
//! Each tetrahedral element uses a corotational formulation:
//!
//! 1. Compute the deformation gradient **F**.
//! 2. Polar-decompose F = R * S to extract the rotation **R**.
//! 3. Compute elastic forces: **f = R * Ke * (R^T * x - X0)**.

mod math_helpers;

pub mod corot_functions;
pub mod corot_sim;
pub mod corotational;
pub mod integrators;
pub mod materials;
pub mod stress;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_extended;

// Re-export all public items for backwards compatibility
pub use corot_functions::*;
pub use corot_sim::*;
pub use corotational::*;
pub use integrators::*;
pub use materials::*;
pub use stress::*;
