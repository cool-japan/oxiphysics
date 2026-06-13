// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Direct (pivoting) solvers for the Linear Complementarity Problem (LCP).
//!
//! Given a square matrix `M` and a vector `q`, the LCP asks for vectors
//! `w` and `z` satisfying:
//!
//! ```text
//! w - M z = q
//! w >= 0,  z >= 0,  wᵀ z = 0
//! ```
//!
//! These are the conditions that arise from velocity-level rigid-body
//! contact resolution (Signorini–Fichera): each row `i` is a contact whose
//! normal force `z_i` and post-solve relative normal velocity `w_i` are
//! non-negative and complementary.
//!
//! Two classical pivoting methods are provided, both exact for the
//! positive-semi-definite (PSD) matrices produced by contact assembly:
//!
//! * [`lemke_solve`] — Lemke's complementary pivoting (1965) with a
//!   lexicographic minimum-ratio test for anti-cycling on degenerate
//!   problems.
//! * [`dantzig_solve`] — Dantzig / Cottle principal pivoting (the
//!   incremental "one contact at a time" scheme used by ODE's `dSolveLCP`).
//!
//! The implementation is deliberately **nalgebra-free**: a small dense
//! tableau over plain `f64` arrays ([`LcpTableau`]) backs the pivoting.
//! For the `≤ 64`-contact systems these solvers target, a hand-rolled dense
//! tableau is both simpler and faster than pulling in an external linear
//! algebra crate, and it keeps the crate's "plain `f64` arrays" house style.
//!
//! A [`harness`] of three canonical contact scenes (box stack, wedge,
//! high-mass-ratio) is included; it doubles as the engine for the solver
//! conformance matrix and reports the complementarity residual achieved by
//! each method.

pub mod dantzig;
pub mod dense;
pub mod harness;
pub mod lemke;

pub use dantzig::dantzig_solve;
pub use dense::LcpTableau;
pub use harness::{
    HarnessResult, complementarity_residual, run_harness, scene_box_stack, scene_high_mass_ratio,
    scene_wedge,
};
pub use lemke::lemke_solve;

use thiserror::Error;

/// Errors produced by the direct LCP solvers.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum LcpError {
    /// Lemke's method hit a secondary ray: no leaving variable exists for
    /// the entering column. The LCP is infeasible or unbounded (for a
    /// copositive `M` this means no solution exists).
    #[error("Lemke: ray termination — problem is infeasible or unbounded")]
    RayTermination,

    /// The complementary pivoting loop exceeded its iteration budget without
    /// reaching a complementary basic feasible solution.
    #[error("Lemke: exceeded {max_pivots} pivots without terminating")]
    MaxPivots {
        /// The pivot budget that was exhausted.
        max_pivots: usize,
    },

    /// Dantzig's principal pivoting could not restore feasibility for the
    /// active subproblem at the given contact index.
    #[error("Dantzig: infeasible contact subproblem at contact {i}")]
    DantzigInfeasible {
        /// Index of the contact whose subproblem failed.
        i: usize,
    },

    /// `M` is not square or its dimension disagrees with `q`.
    #[error("LCP dimension mismatch: M is {m}×{m} but q has length {n}")]
    DimensionMismatch {
        /// The side length of the (square) `M` that was supplied.
        m: usize,
        /// The length of the `q` vector.
        n: usize,
    },
}
