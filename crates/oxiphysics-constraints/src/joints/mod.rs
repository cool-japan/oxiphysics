// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Joint constraints: fixed, revolute, prismatic, ball, and spring joints.

mod balljoint_traits;
mod breakablejoint_traits;
mod cablejoint_traits;
mod conejoint_traits;
mod constraint_math;
mod cylindricaljoint_traits;
mod distancejoint_traits;
mod fixedjoint_traits;
pub(crate) mod functions;
mod jointforcemeter_traits;
mod motorjoint_traits;
mod prismaticjoint_traits;
mod revolutejoint_traits;
mod sphericaljoint_traits;
mod springjoint_traits;
#[cfg(test)]
mod test_support;
#[cfg(test)]
mod tests_point_and_orientation;
#[cfg(test)]
mod tests_revolute;
#[cfg(test)]
mod tests_soft_joints;
pub mod types;
mod universaljoint_traits;
mod weldjoint_traits;
mod xpbdjoint_traits;

// Re-export all public types
pub use types::*;
