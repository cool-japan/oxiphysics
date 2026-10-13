//! Auto-generated module structure

/// Essential boundary conditions for MLS approximations.
mod boundary;
#[cfg(test)]
mod boundary_tests;
/// Functions and utilities.
mod functions;
/// Type definitions.
mod types;

// Re-export all types
pub use boundary::*;
pub use functions::*;
pub use types::*;
