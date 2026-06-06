// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reactive transport for the LBM framework.
//!
//! This module provides:
//!
//! - **ConcentrationField**: 2D concentration field with explicit diffusion and periodic BC.
//! - **first_order_reaction**: C *= exp(-k*dt) applied pointwise.
//! - **bimolecular_reaction**: delta_A = -rate*A*B*dt (explicit Euler).
//! - **LbmPassiveScalar**: LBM-based passive scalar advection-diffusion on D2Q9.
//! - **MultiStepReaction**: sequential A->B->C multi-step reaction chains.
//! - **CatalyticSurface**: catalytic surface with adsorption/desorption.
//! - **SpeciesBoundaryCondition**: inlet/outlet/wall BCs for species transport.
//! - **ReactionFrontTracker**: tracks reaction front position in a domain.
//! - **HeatRelease**: coupling between reaction rate and temperature field.

pub mod catalytic;
pub mod combustion;
pub mod concentration;
pub mod extended;
pub mod lbm_reactive;
pub mod species;

// ---------------------------------------------------------------------------
// Re-exports: concentration
// ---------------------------------------------------------------------------
pub use concentration::{
    ConcentrationField, bimolecular_reaction, bimolecular_reaction_with_product,
    first_order_reaction, reversible_reaction,
};

// ---------------------------------------------------------------------------
// Re-exports: species
// ---------------------------------------------------------------------------
#[cfg(test)]
pub(crate) use species::R_GAS;
pub use species::{ArrheniusRate, BimolecularRate, MultiStepReaction, Reaction, Species};

// ---------------------------------------------------------------------------
// Re-exports: catalytic
// ---------------------------------------------------------------------------
pub use catalytic::{
    BoundaryEdge, CatalyticSurface, HeatRelease, ReactionFrontTracker, SpeciesBcType,
    SpeciesBoundaryCondition,
};

// ---------------------------------------------------------------------------
// Re-exports: lbm_reactive
// ---------------------------------------------------------------------------
pub use lbm_reactive::{
    LbmPassiveScalar, ReactionRate, ReactiveLattice, ReactiveLbm, SpeciesTransport, arrhenius_rate,
    compute_reaction_source, damkoehler_number, peclet_number,
};

// ---------------------------------------------------------------------------
// Re-exports: combustion
// ---------------------------------------------------------------------------
pub use combustion::{
    COMBUSTION_STOICH, ChainBranching, CombustionReactor, DELTA_H_COMBUSTION,
    ElementaryReactionLbm, FlameFrontTracker, IgnitionDelayModel, LaminarFlameSpeed,
    MultiSpeciesDiffusion, SemenovExplosion, SpeciesGrid, ZeldovichMechanism, ZeroDReactor,
    apply_heat_source, combustion_indices, species_source_terms,
};

// ---------------------------------------------------------------------------
// Re-exports: extended
// ---------------------------------------------------------------------------
pub use extended::{
    ArrheniusKinetics, GrayScott, MultiComponentField, ReactiveFlowSolver, SpeciesLbm2D,
    TwoSpeciesReaction, equivalence_ratio, lewis_number, mixture_fraction, progress_variable,
    schmidt_number, zeldovich_number,
};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "tests.rs"]
mod tests;

#[cfg(test)]
#[path = "tests_extended.rs"]
mod tests_extended;
