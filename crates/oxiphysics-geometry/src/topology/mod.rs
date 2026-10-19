//! Auto-generated module structure

#[cfg(test)]
mod fixtures;
/// Functions and utilities.
mod functions;
/// Trait implementations for halfedgemesh_traits.
mod halfedgemesh_traits;
/// Simplicial homology over `Z/2Z` (exact Betti numbers).
mod homology;
/// Voronoi-pole approximation of the medial axis of a surface sample.
mod medial_poles;
/// Trait implementations for meshtopology_traits.
mod meshtopology_traits;
/// Trait implementations for quadedgemesh_traits.
mod quadedgemesh_traits;
/// Validated index-form triangle surface used by the surgery operators.
mod surface;
/// Topological surgery operators (edge split/collapse, handle attach/delete).
mod surgery;
#[cfg(test)]
mod surgery_tests;
/// Type definitions.
mod types;
/// Trait implementations for wingededgemesh_traits.
mod wingededgemesh_traits;

// Re-export all types
pub use functions::*;
pub use homology::*;
pub use medial_poles::*;
pub use surgery::*;
pub use types::*;
