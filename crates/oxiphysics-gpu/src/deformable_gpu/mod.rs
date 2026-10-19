//! Deformable GPU simulation module.

mod collisionresponseparams_traits;
mod corotational;
mod deformablegpumesh_traits;
pub mod functions;
mod physicspipelineconfig_traits;
mod stvk;
pub mod types;

// Re-export all public items
pub use types::*;

#[cfg(test)]
mod fem_kernel_tests;
