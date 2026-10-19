//! Auto-generated module structure

/// Functions and utilities.
mod functions;
/// Regression tests of the persistent contact manifolds.
#[cfg(test)]
mod manifold_tests;
/// Type definitions.
mod types;

// Re-export all types
pub use functions::*;
pub use types::*;
