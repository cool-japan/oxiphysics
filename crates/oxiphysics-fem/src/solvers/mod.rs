//! Auto-generated module structure

pub mod amg;
pub mod assembly_coloring;
pub mod functions;
pub mod functions_2;
pub mod solvererror_traits;
pub mod types;

// Re-export all types
pub use amg::{
    AmgClassical, AmgHierarchy, AmgLevel, AmgPreconditioner, CycleKind, Preconditioner,
    SmoothedAggregationAmg,
};
pub use functions::*;
pub use functions_2::*;
pub use types::*;
