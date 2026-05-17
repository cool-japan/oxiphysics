// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrostatics module including simplified PME and reaction-field methods.
//!
//! Additional features:
//! - Generalized Born (GB) solvation model
//! - Linearized Poisson-Boltzmann model
//! - Charge group handling for cutoff-based methods
//! - Wolf summation for long-range electrostatics
//!
//! Physical units:
//! - energy  : kJ mol^-1
//! - distance: angstrom
//! - charge  : electron units (e)

mod coulomb;
mod multipole;
pub mod pme;
mod polarization;
mod reaction_field;
mod solvation;

pub use coulomb::*;
pub use multipole::*;
pub use pme::{
    ewald_real_space_energy, pme_reciprocal_energy, pme_reciprocal_forces, pme_self_energy,
};
pub use polarization::*;
pub use reaction_field::*;
pub use solvation::*;
