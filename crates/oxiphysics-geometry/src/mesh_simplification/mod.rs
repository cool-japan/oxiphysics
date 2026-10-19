//! Triangle mesh simplification.
//!
//! - [`qem_simplify`] / [`qem_simplify_with_seams`]: quadric error metric edge
//!   collapse (Garland & Heckbert) with a manifold-preserving link condition,
//!   normal-flip and normal-deviation checks, boundary and texture seam handling
//!   ([`QemConfig`], [`TextureSeamFlags`]).
//! - [`check_normal_deviation`]: the normal test of a collapse on a
//!   [`HalfEdgeMesh`].
//! - [`vertex_clustering`], [`edge_length_decimation`]: cheaper simplifiers.
//! - Quality and topology helpers: [`validate_topology`], [`compute_mesh_stats`],
//!   [`hausdorff_symmetric`], [`euler_characteristic`] and others.

pub mod collapsecandidate_traits;
pub mod functions;
pub mod functions_2;
mod qem;
#[cfg(test)]
mod qem_tests;
pub mod textureseamflags_traits;
pub mod types;

// Re-export all types
pub use functions::*;
pub use types::*;
