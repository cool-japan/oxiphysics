//! Auto-generated module structure

pub mod areaconstraint_traits;
pub mod balloon_volume;
pub mod balloonpressureconstraint_traits;
mod bending_hinge;
pub mod bendingconstraint_traits;
pub mod collisionconstraint_traits;
pub mod distanceconstraint_traits;
pub mod functions;
pub mod isometricbendingconstraint_traits;
pub mod isometricbendingconstraintv2_traits;
pub mod lraconstraint_traits;
pub mod motorconstraint_traits;
pub mod neohookeanconstraint_traits;
pub mod rigidbindingconstraint_traits;
pub mod shapematchingconstraint_traits;
pub mod strainlimitingconstraint_traits;
pub mod surfacetensionconstraint_traits;
pub mod types;
pub mod volumeconstraint_traits;

#[cfg(test)]
mod bending_tests;
#[cfg(test)]
mod neohookean_tests;
#[cfg(test)]
mod reset_lambda_tests;

// Re-export all types
pub use balloon_volume::BalloonVolumeConstraint;
pub use functions::*;
pub use types::*;
