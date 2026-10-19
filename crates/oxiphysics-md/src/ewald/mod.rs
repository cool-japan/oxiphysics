// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ewald summation for periodic electrostatics.
//!
//! The Ewald method splits the Coulomb energy of a charge-neutral periodic
//! system into an `erfc`-screened real-space sum, a reciprocal-space sum over
//! the k-vectors of the box and a self-energy correction.
//!
//! - [`EwaldSummation`]: the complete sum for an orthorhombic box.  It returns
//!   energies and analytic forces (the negative gradients of the energies),
//!   supports tinfoil and vacuum boundary conditions (dipole correction) and
//!   provides energy and force decompositions.
//! - [`EwaldSum`]: a compact cubic-box variant with Madelung-constant helpers.
//! - [`EwaldForceIntegrator`]: velocity-Verlet MD driven by the full Ewald
//!   forces.
//! - [`PmeElectrostatics`]: a deliberately reduced force field that evaluates
//!   only the screened real-space part and the self-energy correction, with
//!   *no* reciprocal-space term.  The smooth particle-mesh Ewald method
//!   (B-spline charge spreading and FFT) lives in [`crate::electrostatics::pme`].
//! - Parameter selection and error control: [`EwaldParams::optimize`],
//!   [`AlphaOptimizer`], [`EwaldAccuracyBudget`] and the truncation-error
//!   estimates of [`EwaldSummation`].
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
