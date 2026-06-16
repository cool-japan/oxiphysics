//! Auto-generated module structure

pub mod constraintparams_traits;
pub mod cubic_toi;
pub mod deformcollisionconfig_traits;
pub mod functions;
pub mod penaltyparams_traits;
pub mod types;

// Re-export all types
pub use cubic_toi::{cubic_toi_edge_edge, cubic_toi_vertex_face, solve_cubic_in_01};
pub use functions::*;
pub use types::*;
