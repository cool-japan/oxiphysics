//! Auto-generated module structure

mod box_manifold;
pub mod compoundshape_traits;
mod convex_manifold;
pub mod dispatchconfig_traits;
pub mod dispatchqueue_traits;
mod functions;
pub mod narrowphasedispatcher_traits;
pub mod shapefeaturecache_traits;
pub mod speculativeconfig_traits;
mod types;

// Re-export all types
pub use box_manifold::*;
pub use convex_manifold::*;
pub use functions::*;
pub use types::*;
