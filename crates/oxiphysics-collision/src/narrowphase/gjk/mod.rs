//! Auto-generated module structure

pub(crate) mod distance;
mod functions;
mod functions_2;
pub mod gjksolver_traits;
pub mod mpr;
mod ray_cast;
pub mod simplex_traits;
mod types;
pub mod warmstartgjk_traits;

// Re-export all types
pub use functions::*;
pub use mpr::{mpr_contact, mpr_full};
pub use ray_cast::gjk_ray_cast;
pub use types::*;
