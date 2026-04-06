//! Medical visualization for the OxiPhysics viz crate.
//!
//! Provides volume rendering (ray casting, MIP), isosurface extraction, slice planes,
//! color transfer functions, anatomy atlases, measurement tools, heatmap overlays,
//! vessel tree rendering, and electrophysiology activation maps.

pub mod colortransferfunc_traits;
pub mod functions;
pub(crate) mod helpers;
pub mod incisionpath_traits;
pub mod measurementtool_traits;
pub mod types_core;
pub mod types_ext;

// Re-export all public items
pub use functions::*;
pub use types_core::*;
pub use types_ext::*;
