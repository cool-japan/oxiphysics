//! Auto-generated module structure

pub mod batchnarrowphase_traits;
pub mod contactfilter_traits;
pub mod dispatch;
pub mod epa;
pub mod functions;
pub mod functions_2;
pub mod gjk;
pub mod gjk_cache;
pub mod manifold;
pub mod obb_sat;
pub mod primitives;
pub mod simple;
pub mod specialized;
pub mod types;

// Re-export all types
// NOTE: contactfilter_traits re-exports specialized::*, which conflicts with
// simple::sphere_sphere. We skip the contactfilter_traits glob and re-export
// specialized explicitly via its own module.
pub use dispatch::*;
pub use epa::*;
pub use functions::*;
pub use functions_2::*;
pub use gjk::*;
pub use gjk_cache::*;
pub use manifold::*;
pub use obb_sat::*;
// NOTE: primitives::ClosestFeature conflicts with gjk::ClosestFeature.
// We skip the glob re-export for primitives; use narrowphase::primitives:: explicitly.
pub use types::*;
