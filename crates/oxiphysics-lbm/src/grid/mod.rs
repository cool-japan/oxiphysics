//! Grid and field storage for LBM simulations.
//!
//! Provides 2D and 3D grids that hold distribution functions and
//! macroscopic fields (density, velocity).

pub mod functions;
pub mod grid_extended;
pub mod types;

// Test-only module
#[cfg(test)]
mod functions_2;

// Re-export all types
pub use functions::*;
pub use grid_extended::*;
pub use types::*;
