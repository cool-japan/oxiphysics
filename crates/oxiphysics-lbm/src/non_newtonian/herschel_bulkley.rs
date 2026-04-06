// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Herschel-Bulkley fluid models for non-Newtonian LBM.

use super::{CS2, LocalViscosityModel, MU_MAX, MU_MIN, NonNewtonianFluid, NonNewtonianModel};

// ---------------------------------------------------------------------------
// HerschelBulkley
// ---------------------------------------------------------------------------

/// Herschel-Bulkley fluid model.
///
/// ```text
/// stress = tau_y + K * |gamma|^n
/// viscosity = stress / |gamma|  (for |gamma| > 0)
/// ```
///
/// When `tau_y = 0`, this reduces to a power-law fluid.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct HerschelBulkley {
    /// Consistency coefficient K.
    pub consistency_k: f64,
    /// Flow behaviour index n.
    pub n: f64,
    /// Yield stress.
    pub tau_y: f64,
}

#[allow(dead_code)]
impl HerschelBulkley {
    /// Create a new Herschel-Bulkley fluid.
    pub fn new(consistency_k: f64, n: f64, tau_y: f64) -> Self {
        Self {
            consistency_k,
            n,
            tau_y,
        }
    }

    /// Shear stress at the given shear rate.
    pub fn stress(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        self.tau_y + self.consistency_k * gamma.powf(self.n)
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

impl LocalViscosityModel for HerschelBulkley {
    fn viscosity(&self, gamma_dot: f64) -> f64 {
        HerschelBulkley::viscosity(self, gamma_dot)
    }
}

impl NonNewtonianModel for HerschelBulkley {
    fn viscosity(&self, shear_rate: f64) -> f64 {
        HerschelBulkley::viscosity(self, shear_rate)
    }
}

// ---------------------------------------------------------------------------
// HerschelBulkleyFluid -- extended version with critical shear rate
// ---------------------------------------------------------------------------

/// Extended Herschel-Bulkley model with additional analysis methods.
///
/// The Herschel-Bulkley model combines yield stress with power-law viscoplasticity:
///
/// ```text
/// tau = tau_y + K |gamma|^n     (for tau > tau_y)
/// gamma = 0                     (for tau <= tau_y)
/// ```
///
/// When `n = 1`, reduces to Bingham plastic.
/// When `tau_y = 0`, reduces to power-law fluid.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct HerschelBulkleyFluid {
    /// Yield stress tau_y (Pa in physical units).
    pub tau_y: f64,
    /// Consistency index K (Pa*s^n).
    pub consistency_k: f64,
    /// Power-law exponent n.
    pub n: f64,
}

#[allow(dead_code)]
impl HerschelBulkleyFluid {
    /// Create a new Herschel-Bulkley fluid.
    pub fn new(tau_y: f64, consistency_k: f64, n: f64) -> Self {
        Self {
            tau_y,
            consistency_k,
            n,
        }
    }

    /// Shear stress at the given shear rate (yielded region).
    pub fn stress(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        self.tau_y + self.consistency_k * gamma.powf(self.n)
    }

    /// Effective viscosity `mu_eff = tau / |gamma|` (yielded region).
    pub fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        if gamma < 1e-15 {
            return MU_MAX;
        }
        (self.stress(shear_rate) / gamma).clamp(MU_MIN, MU_MAX)
    }

    /// Critical shear rate: the shear rate at which the fluid transitions
    /// from rigid to yielded behaviour.  For HB this is `gamma_c = (tau_y / K)^(1/n)`.
    ///
    /// Returns `0.0` if `tau_y = 0` (no yield stress).
    pub fn critical_shear_rate(&self) -> f64 {
        if self.tau_y < 1e-30 || self.consistency_k < 1e-30 {
            return 0.0;
        }
        (self.tau_y / self.consistency_k).powf(1.0 / self.n)
    }

    /// Check whether the material has yielded at the given shear rate.
    pub fn is_yielded(&self, shear_rate: f64) -> bool {
        shear_rate.abs() > self.critical_shear_rate()
    }

    /// Bingham number `Bn = tau_y / (K * gamma_ref^n)` -- ratio of yield stress to
    /// viscoplastic stress at a reference shear rate.
    pub fn bingham_number(&self, gamma_ref: f64) -> f64 {
        if gamma_ref < 1e-30 || self.consistency_k < 1e-30 {
            return f64::INFINITY;
        }
        self.tau_y / (self.consistency_k * gamma_ref.powf(self.n))
    }

    /// Effective relaxation time for LBM: `tau = 0.5 + mu_eff / cs^2`.
    pub fn lbm_relaxation_time(&self, shear_rate: f64) -> f64 {
        0.5 + self.effective_viscosity(shear_rate) / CS2
    }

    /// Generalized (apparent) viscosity at the given shear rate via
    /// the Metzner-Reed definition: `mu_app = tau / |gamma|`.
    pub fn apparent_viscosity(&self, shear_rate: f64) -> f64 {
        self.effective_viscosity(shear_rate)
    }

    /// Pipe-flow mean velocity for a Herschel-Bulkley fluid in a tube of
    /// radius `R` under pressure gradient `G = dP / L`.
    ///
    /// Uses a simplified numerical integration approach.
    pub fn pipe_mean_velocity(&self, radius: f64, pressure_gradient: f64) -> f64 {
        let r_p = if pressure_gradient < 1e-30 {
            return 0.0;
        } else {
            2.0 * self.tau_y / pressure_gradient
        };
        if r_p >= radius {
            return 0.0; // no flow -- pressure gradient insufficient to yield
        }
        let n = self.n;
        let k = self.consistency_k;
        let g = pressure_gradient;
        // Numerical integration using 200 trapezoid intervals
        let n_pts = 200usize;
        let dr = radius / n_pts as f64;
        let mut integral = 0.0_f64;
        for j in 0..n_pts {
            let r1 = j as f64 * dr;
            let r2 = r1 + dr;
            let u1 = if r1 < r_p {
                let factor = (g / (2.0 * k)).powf(1.0 / n);
                let integrate_annulus = |r_inner: f64| -> f64 {
                    if r_inner <= r_p {
                        return (radius - r_p).powf((n + 1.0) / n) * n / (n + 1.0);
                    }
                    let upper =
                        (radius - r_p).powf((n + 1.0) / n) - (r_inner - r_p).powf((n + 1.0) / n);
                    factor * upper * n / (n + 1.0)
                };
                r1 * 2.0 * std::f64::consts::PI * factor * integrate_annulus(r_p)
            } else {
                let factor = (g / (2.0 * k)).powf(1.0 / n);
                let u = factor
                    * (n / (n + 1.0))
                    * ((radius - r_p).powf((n + 1.0) / n) - (r1 - r_p).powf((n + 1.0) / n));
                r1 * 2.0 * std::f64::consts::PI * u
            };
            let u2 = if r2 < r_p {
                let factor = (g / (2.0 * k)).powf(1.0 / n);
                let integrate_annulus = |r_inner: f64| -> f64 {
                    if r_inner <= r_p {
                        return (radius - r_p).powf((n + 1.0) / n) * n / (n + 1.0);
                    }
                    let upper =
                        (radius - r_p).powf((n + 1.0) / n) - (r_inner - r_p).powf((n + 1.0) / n);
                    factor * upper * n / (n + 1.0)
                };
                r2 * 2.0 * std::f64::consts::PI * factor * integrate_annulus(r_p)
            } else {
                let factor = (g / (2.0 * k)).powf(1.0 / n);
                let u = factor
                    * (n / (n + 1.0))
                    * ((radius - r_p).powf((n + 1.0) / n) - (r2 - r_p).powf((n + 1.0) / n));
                r2 * 2.0 * std::f64::consts::PI * u.max(0.0)
            };
            integral += 0.5 * (u1 + u2) * dr;
        }
        // Divide by pipe area to get mean velocity
        let area = std::f64::consts::PI * radius * radius;
        (integral / area).max(0.0)
    }
}

impl NonNewtonianFluid for HerschelBulkleyFluid {
    fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        self.effective_viscosity(shear_rate)
    }
}

impl LocalViscosityModel for HerschelBulkleyFluid {
    fn viscosity(&self, gamma_dot: f64) -> f64 {
        self.effective_viscosity(gamma_dot)
    }
}
