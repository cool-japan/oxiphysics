// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Algebraic Multigrid (AMG) solver module.

pub mod aggregation;
pub mod classical;
pub mod cycle;
pub mod galerkin;
pub mod graph;
pub mod near_null_space;
pub mod preconditioner;
pub mod smoothed_aggregation;
pub mod smoothers;

pub use classical::AmgClassical;
pub use cycle::{AmgHierarchy, AmgLevel, CycleKind};
pub use preconditioner::{AmgPreconditioner, Preconditioner};
pub use smoothed_aggregation::SmoothedAggregationAmg;
