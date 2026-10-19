// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Curved boundary treatment for the Lattice Boltzmann Method.
//!
//! Every scheme in this module works *link-wise*: a boundary link connects a
//! fluid node `x_f` to a solid node `x_s = x_f + c_i` and is cut by the wall
//! at `x_w = x_f + q c_i` with `q ∈ (0, 1]`.  After streaming, the population
//! `f_ī(x_f, t+1)` travelling back from the solid node (`ī` is the direction
//! opposite to `i`) is unknown; the schemes reconstruct it from post-collision
//! populations `f̃` sampled along the link line:
//!
//! - [`bouzidi_bc`] — linear Bouzidi–Firdaouss–Lallemand (BFL) interpolated
//!   bounce-back (Bouzidi et al., *Phys. Fluids* 13, 2001).
//! - [`bouzidi_quadratic_bc`] — quadratic BFL interpolation.
//! - [`interpolated_bounce_back`] — the unified single-formula linear scheme of
//!   Yu, Mei, Luo & Shyy (*Prog. Aerosp. Sci.* 39, 2003).
//! - [`filippova_hanel`] — Filippova–Hänel fictitious-equilibrium scheme with
//!   the Mei–Luo–Shyy (*J. Comput. Phys.* 155, 1999) stability improvement.
//! - [`moving_wall_bc`] — half-way bounce-back with the Ladd momentum
//!   correction `−2 w_i ρ (c_i·u_w)/c_s²`.
//! - [`CurvedBoundarySet`] — a collection of [`CurvedWallNode`] links that
//!   applies any [`CurvedBcMethod`] to a whole distribution field, including a
//!   moment-based reconstruction that enforces mass and momentum constraints
//!   on the unknown populations of each boundary node.
//! - [`RotatingCylinder`] — a 2-D rotating cylinder whose boundary links carry
//!   the local tangential wall velocity.

mod cylinder;
mod kernels;
mod lattice;
mod set;

#[cfg(test)]
mod channel_tests;
#[cfg(test)]
mod tests;

pub use cylinder::RotatingCylinder;
pub use kernels::{
    LinkPopulations, WallLink, bouzidi_bc, bouzidi_quadratic_bc, bouzidi_wall_correction, box_sdf,
    d2q9_equilibrium, direct_forcing, fictitious_equilibrium, filippova_hanel,
    half_way_bounce_back, interpolated_bounce_back, mei_luo_shyy_chi, moving_wall_bc,
    quadratic_interpolation_bc, sphere_sdf,
};
pub use lattice::CurvedLattice;
pub use set::{CurvedBoundarySet, CurvedWallNode, LatticeDomain};

// ---------------------------------------------------------------------------
// Enums
// ---------------------------------------------------------------------------

/// Method used for curved boundary treatment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvedBcMethod {
    /// Filippova–Hänel fictitious-equilibrium bounce-back with the Mei–Luo–Shyy
    /// choice of `u_bf` and `χ` (second-order accurate).
    FilippovaHanel,
    /// Linear Bouzidi–Firdaouss–Lallemand interpolated bounce-back.
    Bouzidi,
    /// Quadratic Bouzidi–Firdaouss–Lallemand interpolated bounce-back using
    /// three nodes along the link line.
    QuadraticInterpolation,
    /// Moment-based reconstruction: linear BFL estimates of the unknown
    /// populations are corrected by the minimal weighted perturbation that
    /// keeps the node's mass and imposes the wall-interpolated momentum.
    MomentBased,
}
