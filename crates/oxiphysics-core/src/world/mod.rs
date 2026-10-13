//! Auto-generated module structure

/// Trait implementations for body_traits.
mod body_traits;
/// Position-based (XPBD) constraint solver of the world step.
mod constraint_solver;
/// Functions and utilities.
mod functions;
/// Trait implementations for physicsworld_traits.
mod physicsworld_traits;
/// Type definitions.
mod types;

#[cfg(test)]
mod accumulator_tests;
#[cfg(test)]
mod constraint_tests;
#[cfg(test)]
mod rotation_tests;

// Re-export all types
pub use functions::*;
pub use types::*;
