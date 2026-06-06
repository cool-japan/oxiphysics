// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Power-law fluid models for non-Newtonian LBM.

use super::{CS2, LocalViscosityModel, MU_MAX, MU_MIN, NonNewtonianFluid, NonNewtonianModel};

// ---------------------------------------------------------------------------
// PowerLawFluid
// ---------------------------------------------------------------------------

/// Power-law (Ostwald-de Waele) fluid model.
///
/// The effective viscosity follows:
///
/// ```text
/// mu_eff(gamma) = K * |gamma|^(n-1)
/// ```
///
/// clamped to `[mu_min, mu_max]`.
///
/// * `n < 1` -- shear-thinning (pseudo-plastic)
/// * `n = 1` -- Newtonian (reduces to mu = K)
/// * `n > 1` -- shear-thickening (dilatant)
#[derive(Debug, Clone, Copy)]
pub struct PowerLawFluid {
    /// Flow behaviour index.
    pub n: f64,
    /// Consistency coefficient K (units: Pa*s^n in physical space).
    pub consistency_k: f64,
}

impl PowerLawFluid {
    /// Create a new power-law fluid with behaviour index `n` and consistency `k`.
    pub fn new(n: f64, consistency_k: f64) -> Self {
        Self { n, consistency_k }
    }

    /// Effective viscosity at the given scalar shear rate `|gamma|`.
    ///
    /// Returns `K * |gamma|^(n-1)`, clamped to `[mu_min, mu_max]`.
    pub fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        // Avoid division by zero when n < 1 and gamma -> 0.
        let mu = if gamma < 1e-15 {
            if self.n < 1.0 {
                MU_MAX
            } else {
                self.consistency_k
            }
        } else {
            self.consistency_k * gamma.powf(self.n - 1.0)
        };
        mu.clamp(MU_MIN, MU_MAX)
    }

    /// Local relaxation time `tau` from the shear rate.
    ///
    /// Uses `tau = 0.5 + mu_eff / (rho * cs^2)` with `rho = 1` in lattice units.
    pub fn local_tau(&self, shear_rate: f64) -> f64 {
        let mu_eff = self.effective_viscosity(shear_rate);
        0.5 + mu_eff / CS2
    }

    /// Compute the shear stress at the given shear rate.
    ///
    /// `stress = K * |gamma|^n`
    pub fn stress(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        self.consistency_k * gamma.powf(self.n)
    }

    /// Alias: viscosity at the given shear rate (same as `effective_viscosity`).
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        self.effective_viscosity(shear_rate)
    }

    /// Returns `true` if the fluid is shear-thinning (n < 1).
    pub fn is_shear_thinning(&self) -> bool {
        self.n < 1.0
    }
}

impl NonNewtonianFluid for PowerLawFluid {
    fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        self.effective_viscosity(shear_rate)
    }
}

impl LocalViscosityModel for PowerLawFluid {
    fn viscosity(&self, gamma_dot: f64) -> f64 {
        self.effective_viscosity(gamma_dot)
    }
}

// ---------------------------------------------------------------------------
// PowerLaw struct (short-named variant)
// ---------------------------------------------------------------------------

/// Power-law (Ostwald-de Waele) fluid -- short-named variant for ergonomic use.
///
/// `mu_eff(gamma) = k * |gamma|^(n-1)`, clamped to `[mu_min, mu_max]`.
///
/// * `n < 1` -- shear-thinning
/// * `n = 1` -- Newtonian (reduces to mu = k)
/// * `n > 1` -- shear-thickening
#[derive(Debug, Clone, Copy)]
pub struct PowerLaw {
    /// Consistency coefficient k (Pa*s^n in physical units).
    pub k: f64,
    /// Flow behaviour index n.
    pub n: f64,
}

impl PowerLaw {
    /// Create a new power-law fluid.
    pub fn new(k: f64, n: f64) -> Self {
        Self { k, n }
    }

    /// Effective viscosity at scalar shear rate `|gamma|`.
    pub fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let mu = if gamma < 1e-15 {
            if self.n < 1.0 { MU_MAX } else { self.k }
        } else {
            self.k * gamma.powf(self.n - 1.0)
        };
        mu.clamp(MU_MIN, MU_MAX)
    }

    /// LBM relaxation time from shear rate: `tau = 0.5 + mu_eff / cs^2`.
    pub fn local_tau(&self, shear_rate: f64) -> f64 {
        0.5 + self.effective_viscosity(shear_rate) / CS2
    }

    /// Returns `true` if shear-thinning (n < 1).
    pub fn is_shear_thinning(&self) -> bool {
        self.n < 1.0
    }
}

impl NonNewtonianModel for PowerLaw {
    fn viscosity(&self, shear_rate: f64) -> f64 {
        self.effective_viscosity(shear_rate)
    }
}
