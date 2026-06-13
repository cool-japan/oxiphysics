//! Auto-generated module structure

pub mod amg;
pub mod assembly_coloring;
pub mod block_schur_pc;
pub mod functions;
pub mod functions_2;
pub mod matrix_free;
pub mod solvererror_traits;
pub mod types;

// Re-export all types
pub use amg::{
    AmgClassical, AmgHierarchy, AmgLevel, AmgPreconditioner, CycleKind, Preconditioner,
    SmoothedAggregationAmg, chebyshev_smoother,
};
pub use block_schur_pc::{
    BlockSchurPreconditioner, block_schur_gmres, extract_velocity_block, gmres_left_preconditioned,
    lump_pressure_mass,
};
pub use functions::*;
pub use functions_2::*;
pub use matrix_free::{MatrixFreeError, MatrixFreeOperator, MatrixFreePcg, SumFactPoisson};
pub use types::*;
