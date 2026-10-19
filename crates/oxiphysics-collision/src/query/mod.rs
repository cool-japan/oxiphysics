//! Auto-generated module structure

/// Trait implementations for compoundquery_traits.
mod compoundquery_traits;
/// Functions and utilities.
mod functions;
/// Closest-point projection onto the convex hull of a vertex set.
mod hull_projection;
/// Trait implementations for querybox_traits.
mod querybox_traits;
/// Trait implementations for querycapsule_traits.
mod querycapsule_traits;
/// Trait implementations for querysphere_traits.
mod querysphere_traits;
/// Type definitions.
mod types;

// Re-export all types
pub use functions::*;
pub use hull_projection::project_point_convex;
pub use types::*;
