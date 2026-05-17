// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fluid FEM and stabilized methods.
//!
//! Re-exports all types from the two sub-modules:
//! - [`fluid_fem_stabilized`]: SUPG/PSPG/VMS stabilized methods, ALE, level-set
//! - [`fluid_fem_stokes`]: Taylor-Hood Stokes FEM and 1-D Navier-Stokes FEM

pub mod fluid_fem_stabilized;
pub mod fluid_fem_stokes;

pub use fluid_fem_stabilized::*;
pub use fluid_fem_stokes::*;
