//! Auto-generated module structure

pub mod compoundshape_traits;
pub mod dispatchconfig_traits;
pub mod dispatchqueue_traits;
mod box_manifold;
mod functions;
pub mod narrowphasedispatcher_traits;
pub mod shapefeaturecache_traits;
pub mod speculativeconfig_traits;
mod types;

// Re-export all types
pub use functions::*;
pub use box_manifold::*;
pub use types::*;
