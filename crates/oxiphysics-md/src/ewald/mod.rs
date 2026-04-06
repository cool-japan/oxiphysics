// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle-Mesh Ewald (PME) electrostatics.
//!
//! Implements a simplified PME approach where only the real-space (direct-sum)
//! part with Ewald damping is computed.  The reciprocal-space contribution is
//! accounted for analytically via a self-energy correction term.
//!
//! Additional features:
//! - Ewald parameter optimization (optimal alpha, error estimation)
//! - Tinfoil and vacuum boundary conditions (dipole correction)
//! - Smooth PME with B-spline charge spreading
//! - Force and energy decomposition utilities
//!
//! Physical units throughout:
//! - energy  : kJ mol^-1
//! - distance: angstrom
//! - charge  : electron units (e)

mod params;
mod real_space;
mod reciprocal;
mod summation;

pub use params::*;
pub use real_space::*;
pub use reciprocal::*;
pub use summation::*;
