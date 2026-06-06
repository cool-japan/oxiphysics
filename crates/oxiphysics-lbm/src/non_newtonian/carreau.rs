// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Carreau, Cross, Casson, and Carreau-Yasuda fluid models.

use super::{CS2, LocalViscosityModel, MU_MAX, MU_MIN, NonNewtonianFluid, NonNewtonianModel};

// ---------------------------------------------------------------------------
// CarreauFluid
// ---------------------------------------------------------------------------

/// Carreau fluid model.
///
/// ```text
/// mu(gamma) = eta_inf + (eta_0 - eta_inf) * (1 + (lambda * gamma)^2)^((n-1)/2)
/// ```
#[derive(Debug, Clone, Copy)]
pub struct CarreauFluid {
    /// Zero-shear-rate viscosity.
    pub eta_0: f64,
    /// Infinite-shear-rate viscosity.
    pub eta_inf: f64,
    /// Relaxation time lambda.
    pub lambda: f64,
    /// Power-law index.
    pub n: f64,
}

impl CarreauFluid {
    /// Create a new Carreau fluid.
    pub fn new(eta_0: f64, eta_inf: f64, lambda: f64, n: f64) -> Self {
        Self {
            eta_0,
            eta_inf,
            lambda,
            n,
        }
    }

    /// Effective viscosity at the given shear rate.
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let lg = self.lambda * gamma;
        self.eta_inf + (self.eta_0 - self.eta_inf) * (1.0 + lg * lg).powf((self.n - 1.0) / 2.0)
    }
}

impl LocalViscosityModel for CarreauFluid {
    fn viscosity(&self, gamma_dot: f64) -> f64 {
        CarreauFluid::viscosity(self, gamma_dot)
    }
}

// ---------------------------------------------------------------------------
// CrossFluid
// ---------------------------------------------------------------------------

/// Cross fluid model.
///
/// ```text
/// mu(gamma) = eta_inf + (eta_0 - eta_inf) / (1 + (K * gamma)^m)
/// ```
///
/// The Cross model is an alternative to the Carreau model that provides
/// a better fit for some polymer solutions.
#[derive(Debug, Clone, Copy)]
pub struct CrossFluid {
    /// Zero-shear-rate viscosity.
    pub eta_0: f64,
    /// Infinite-shear-rate viscosity.
    pub eta_inf: f64,
    /// Cross time constant K.
    pub k: f64,
    /// Cross rate constant m (dimensionless exponent).
    pub m: f64,
}

impl CrossFluid {
    /// Create a new Cross fluid.
    pub fn new(eta_0: f64, eta_inf: f64, k: f64, m: f64) -> Self {
        Self {
            eta_0,
            eta_inf,
            k,
            m,
        }
    }

    /// Effective viscosity at the given shear rate.
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let kg = self.k * gamma;
        self.eta_inf + (self.eta_0 - self.eta_inf) / (1.0 + kg.powf(self.m))
    }
}

// ---------------------------------------------------------------------------
// CassonFluid
// ---------------------------------------------------------------------------

/// Casson fluid model.
///
/// ```text
/// sqrt(stress) = sqrt(tau_y) + sqrt(mu_inf * gamma)
/// stress = (sqrt(tau_y) + sqrt(mu_inf * gamma))^2
/// viscosity = stress / gamma
/// ```
///
/// Used for blood and some food products.
#[derive(Debug, Clone, Copy)]
pub struct CassonFluid {
    /// Yield stress tau_y.
    pub tau_y: f64,
    /// Infinite-shear-rate viscosity mu_inf.
    pub mu_inf: f64,
}

impl CassonFluid {
    /// Create a new Casson fluid.
    pub fn new(tau_y: f64, mu_inf: f64) -> Self {
        Self { tau_y, mu_inf }
    }

    /// Shear stress at the given shear rate.
    pub fn stress(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let s = self.tau_y.sqrt() + (self.mu_inf * gamma).sqrt();
        s * s
    }

    /// Effective viscosity at the given shear rate.
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        if gamma < 1e-15 {
            return MU_MAX;
        }
        (self.stress(shear_rate) / gamma).clamp(MU_MIN, MU_MAX)
    }
}

// ---------------------------------------------------------------------------
// Carreau struct (short-named variant)
// ---------------------------------------------------------------------------

/// Carreau shear-thinning fluid -- short-named variant.
///
/// `mu(gamma) = mu_inf + (mu_0 - mu_inf) * (1 + (lambda*gamma)^2)^((n-1)/2)`
#[derive(Debug, Clone, Copy)]
pub struct Carreau {
    /// Zero-shear-rate viscosity mu_0.
    pub mu_0: f64,
    /// Infinite-shear-rate viscosity mu_inf.
    pub mu_inf: f64,
    /// Relaxation time constant lambda (s).
    pub lambda: f64,
    /// Power-law index n (< 1 for shear-thinning).
    pub n: f64,
}

impl Carreau {
    /// Create a new Carreau fluid.
    pub fn new(mu_0: f64, mu_inf: f64, lambda: f64, n: f64) -> Self {
        Self {
            mu_0,
            mu_inf,
            lambda,
            n,
        }
    }

    /// Effective viscosity at scalar shear rate `|gamma|`.
    pub fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let lg = self.lambda * gamma;
        self.mu_inf + (self.mu_0 - self.mu_inf) * (1.0 + lg * lg).powf((self.n - 1.0) / 2.0)
    }

    /// LBM relaxation time from shear rate: `tau = 0.5 + mu_eff / cs^2`.
    pub fn local_tau(&self, shear_rate: f64) -> f64 {
        0.5 + self.effective_viscosity(shear_rate) / CS2
    }
}

impl NonNewtonianModel for Carreau {
    fn viscosity(&self, shear_rate: f64) -> f64 {
        self.effective_viscosity(shear_rate)
    }
}

// ---------------------------------------------------------------------------
// CarreauYasudaFluid
// ---------------------------------------------------------------------------

/// Carreau-Yasuda fluid model.
///
/// Generalization of the Carreau model with an extra parameter `a` (the Yasuda
/// exponent) controlling the width of the transition region:
///
/// ```text
/// eta = eta_inf + (eta_0 - eta_inf) * (1 + (lambda * gamma)^a)^((n-1)/a)
/// ```
///
/// Setting `a = 2` recovers the standard Carreau model.
/// The Carreau-Yasuda model is widely used for polymer melts and blood.
#[derive(Debug, Clone, Copy)]
pub struct CarreauYasudaFluid {
    /// Zero-shear-rate viscosity eta_0.
    pub eta_0: f64,
    /// Infinite-shear-rate viscosity eta_inf.
    pub eta_inf: f64,
    /// Relaxation/time constant lambda (s).
    pub lambda: f64,
    /// Yasuda exponent `a` (dimensionless, typically 2).
    pub a: f64,
    /// Power-law index `n` (dimensionless).
    pub n: f64,
}

impl CarreauYasudaFluid {
    /// Create a new Carreau-Yasuda fluid.
    pub fn new(eta_0: f64, eta_inf: f64, lambda: f64, a: f64, n: f64) -> Self {
        Self {
            eta_0,
            eta_inf,
            lambda,
            a,
            n,
        }
    }

    /// Effective viscosity at scalar shear rate `gamma`.
    ///
    /// `eta = eta_inf + (eta_0 - eta_inf) * (1 + (lambda*gamma)^a)^((n-1)/a)`
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        let lg = self.lambda * gamma;
        // Use powf; safe because lg >= 0 and a > 0.
        let inner = 1.0 + lg.powf(self.a);
        self.eta_inf + (self.eta_0 - self.eta_inf) * inner.powf((self.n - 1.0) / self.a)
    }

    /// Limiting viscosity at zero shear rate: should return eta_0.
    pub fn zero_shear_viscosity(&self) -> f64 {
        self.eta_0
    }

    /// Limiting viscosity at infinite shear rate: approaches eta_inf.
    pub fn infinite_shear_viscosity(&self) -> f64 {
        self.eta_inf
    }

    /// LBM relaxation time from the local shear rate.
    pub fn local_tau(&self, shear_rate: f64) -> f64 {
        0.5 + self.viscosity(shear_rate) / CS2
    }
}

impl NonNewtonianFluid for CarreauYasudaFluid {
    fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        self.viscosity(shear_rate)
    }
}

impl LocalViscosityModel for CarreauYasudaFluid {
    fn viscosity(&self, gamma_dot: f64) -> f64 {
        CarreauYasudaFluid::viscosity(self, gamma_dot)
    }
}
