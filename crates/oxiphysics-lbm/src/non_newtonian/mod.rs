// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Non-Newtonian fluid models for the Lattice Boltzmann Method.
//!
//! Provides power-law, Bingham, Cross, Casson, Carreau, and Herschel-Bulkley
//! fluid rheology, together with a helper struct for maintaining a
//! spatially-varying relaxation time field.
//!
//! Extended features:
//! - Cross model
//! - Casson model
//! - Regularized Bingham model (Papanastasiou)
//! - Apparent viscosity iteration (Picard)
//! - Viscosity field output
//!
//! In LBM the local relaxation time is linked to the local effective kinematic
//! viscosity via:
//!
//! ```text
//! tau = 0.5 + nu_eff / cs^2   where cs^2 = 1/3
//! ```
//!
//! so  `tau = 0.5 + 3 * nu_eff`.

mod bingham;
mod carreau;
mod herschel_bulkley;
mod lbm_integration;
mod power_law;
mod viscoelastic;

// ---------------------------------------------------------------------------
// Module-wide constants
// ---------------------------------------------------------------------------

/// Speed of sound squared: cs^2 = 1/3.
const CS2: f64 = 1.0 / 3.0;

/// Minimum effective viscosity (prevents tau from going below 0.5 + epsilon).
const MU_MIN: f64 = 1e-6;
/// Maximum effective viscosity (prevents tau from diverging).
const MU_MAX: f64 = 1e6;
/// Very large viscosity returned for the rigid (un-yielded) Bingham zone.
const RIGID_MU: f64 = 1e8;

// ---------------------------------------------------------------------------
// Trait: NonNewtonianFluid
// ---------------------------------------------------------------------------

/// Common interface for non-Newtonian fluid models.
pub trait NonNewtonianFluid {
    /// Return the effective viscosity for the given scalar shear rate.
    fn effective_viscosity(&self, shear_rate: f64) -> f64;

    /// Return the local LBM relaxation time (lattice units, rho = 1).
    fn local_tau(&self, shear_rate: f64) -> f64 {
        0.5 + self.effective_viscosity(shear_rate) / CS2
    }
}

// ---------------------------------------------------------------------------
// Trait: LocalViscosityModel
// ---------------------------------------------------------------------------

/// A trait for computing local effective viscosity from a scalar shear rate.
///
/// Implement this for any rheology model that can be evaluated per-cell
/// in an LBM simulation.
pub trait LocalViscosityModel {
    /// Return the effective (kinematic) viscosity at the given scalar shear rate
    /// `gamma_dot` (1/s in physical units, or 1/timestep in lattice units).
    fn viscosity(&self, gamma_dot: f64) -> f64;
}

// ---------------------------------------------------------------------------
// Trait: NonNewtonianModel
// ---------------------------------------------------------------------------

/// Simplified trait for non-Newtonian fluid models used in LBM collision.
///
/// Distinct from `NonNewtonianFluid` -- this trait uses `viscosity()` as the
/// method name to align with the short-named structs `PowerLaw`, `Carreau`, etc.
pub trait NonNewtonianModel {
    /// Return the effective viscosity at the given scalar shear rate.
    fn viscosity(&self, shear_rate: f64) -> f64;
}

// ---------------------------------------------------------------------------
// Re-exports: power_law
// ---------------------------------------------------------------------------
pub use power_law::{PowerLaw, PowerLawFluid};

// ---------------------------------------------------------------------------
// Re-exports: bingham
// ---------------------------------------------------------------------------
pub use bingham::{
    BinghamFluid, BinghamPipeFlow, BinghamPlastic, RegularizedBingham, bingham_yield_correction,
    is_yielded_2d, von_mises_stress_2d,
};

// ---------------------------------------------------------------------------
// Re-exports: herschel_bulkley
// ---------------------------------------------------------------------------
pub use herschel_bulkley::{HerschelBulkley, HerschelBulkleyFluid};

// ---------------------------------------------------------------------------
// Re-exports: carreau
// ---------------------------------------------------------------------------
pub use carreau::{Carreau, CarreauFluid, CarreauYasudaFluid, CassonFluid, CrossFluid};

// ---------------------------------------------------------------------------
// Re-exports: viscoelastic
// ---------------------------------------------------------------------------
pub use viscoelastic::{
    GiesekusFluid, JohnsonSegalman, MaxwellFluid, OldroydB, PhanThienTanner, RheologyLookupTable,
    RoliePolyModel, ViscoelasticRelaxation, YieldCriterion,
};

// ---------------------------------------------------------------------------
// Re-exports: lbm_integration
// ---------------------------------------------------------------------------
#[cfg(test)]
pub(crate) use lbm_integration::{C9_NN, W9_NN};
pub use lbm_integration::{
    LocalTauLattice, NonNewtonianLBM, PapanastasiouViscoplastic, ThixotropicFluid,
    ThixotropicModel, ViscosityField, apply_non_newtonian_collision, deborah_number,
    effective_relaxation_frequency, effective_relaxation_time, effective_tau,
    generalized_reynolds_power_law, iterate_apparent_viscosity, kolmogorov_scale,
    metzner_reed_reynolds, non_newtonian_bgk_step_2d, oldroyd_b_viscosity,
    shear_rate_from_strain_tensor, smagorinsky_turbulent_viscosity, trouton_ratio,
    turbulent_effective_viscosity, update_relaxation_field, update_tau_from_shear_rates,
    van_driest_mixing_length, viscosity_index, weissenberg_number,
};

// ---------------------------------------------------------------------------
// RheologyModel enum (dispatch)
// ---------------------------------------------------------------------------

/// Rheology model variants for dispatch.
///
/// Enables selecting a fluid model at runtime without dynamic dispatch,
/// covering all models provided in this module.
#[derive(Debug, Clone, Copy)]
pub enum RheologyModel {
    /// Power-law (Ostwald-de Waele) fluid.
    PowerLaw(PowerLawFluid),
    /// Bingham plastic fluid.
    Bingham(BinghamFluid),
    /// Carreau fluid (a=2 Carreau-Yasuda special case).
    Carreau(CarreauFluid),
    /// Carreau-Yasuda fluid.
    CarreauYasuda(CarreauYasudaFluid),
    /// Cross fluid.
    Cross(CrossFluid),
    /// Casson fluid.
    Casson(CassonFluid),
    /// Herschel-Bulkley fluid.
    HerschelBulkley(HerschelBulkley),
    /// Regularized Bingham (Papanastasiou).
    RegularizedBingham(RegularizedBingham),
}

/// Dispatch function: compute the effective viscosity for any `RheologyModel` variant.
///
/// This is the single-entry-point function for non-Newtonian viscosity in the
/// LBM collision kernel; it avoids requiring dynamic dispatch via trait objects
/// in hot inner loops.
pub fn effective_viscosity(model: RheologyModel, shear_rate: f64) -> f64 {
    match model {
        RheologyModel::PowerLaw(m) => m.effective_viscosity(shear_rate),
        RheologyModel::Bingham(m) => m.effective_viscosity(shear_rate),
        RheologyModel::Carreau(m) => m.viscosity(shear_rate),
        RheologyModel::CarreauYasuda(m) => m.viscosity(shear_rate),
        RheologyModel::Cross(m) => m.viscosity(shear_rate),
        RheologyModel::Casson(m) => m.viscosity(shear_rate),
        RheologyModel::HerschelBulkley(m) => m.viscosity(shear_rate),
        RheologyModel::RegularizedBingham(m) => m.viscosity(shear_rate),
    }
}

/// Compute the LBM relaxation time `tau = 0.5 + nu_eff / cs^2` for any rheology model.
pub fn relaxation_time(model: RheologyModel, shear_rate: f64) -> f64 {
    0.5 + effective_viscosity(model, shear_rate) / CS2
}

/// Compute the LBM relaxation frequency `omega = 1/tau` for any rheology model.
pub fn relaxation_frequency(model: RheologyModel, shear_rate: f64) -> f64 {
    1.0 / relaxation_time(model, shear_rate)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "tests.rs"]
mod tests;

#[cfg(test)]
#[path = "tests_extended.rs"]
mod tests_extended;
