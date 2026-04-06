// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! High-level molecular dynamics simulation driver.
//!
//! Provides both a generic [`MdSimulation`] (integrates with the full trait
//! system) and a self-contained `MdSim` / [`MdState`] / [`MdConfig`] trio
//! that operates purely on plain `Vec<[f64;3]>` data and has no external
//! crate dependencies.

pub mod advanced;
pub mod analysis;
pub mod checkpoint;
pub mod config;
pub mod core_sim;
pub mod forces;
pub mod generic_sim;
pub mod plain_sim;
pub mod tracking;

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_advanced;

// Re-export all public items for backwards compatibility
pub use advanced::*;
pub use analysis::*;
pub use checkpoint::*;
pub use config::*;
pub use core_sim::*;
pub use forces::*;
pub use generic_sim::*;
pub use plain_sim::*;
pub use tracking::*;
