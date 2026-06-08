//! Numerical algorithms for physics simulation.
//!
//! Provides root finding, quadrature, finite differences, and special functions
//! used throughout the OxiPhysics engine.

pub mod functions;
pub mod functions_2;
pub mod types;

// Re-export all types
pub use functions::*;
pub use functions_2::*;
pub use types::*;

#[cfg(test)]
mod tests;
