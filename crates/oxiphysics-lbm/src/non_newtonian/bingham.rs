// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Bingham plastic models and yield-surface helpers.

use super::{MU_MAX, MU_MIN, NonNewtonianFluid, NonNewtonianModel, RIGID_MU};

// ---------------------------------------------------------------------------
// BinghamFluid
// ---------------------------------------------------------------------------

/// Bingham plastic fluid model.
///
/// For a Bingham fluid the constitutive relation is:
///
/// ```text
/// if |tau_shear| < tau_yield   ->  rigid (effectively infinite viscosity)
/// else  mu_eff = mu_p + tau_yield / |gamma|
/// ```
///
/// The regularised form used here returns `RIGID_MU` for the un-yielded zone
/// instead of true infinity, which keeps tau finite and numerically tractable.
#[derive(Debug, Clone, Copy)]
pub struct BinghamFluid {
    /// Yield stress tau_y (Pa in physical units; dimensionless in lattice units).
    pub yield_stress: f64,
    /// Plastic viscosity mu_p.
    pub plastic_viscosity: f64,
}

impl BinghamFluid {
    /// Create a new Bingham fluid.
    pub fn new(yield_stress: f64, plastic_viscosity: f64) -> Self {
        Self {
            yield_stress,
            plastic_viscosity,
        }
    }

    /// Effective viscosity at scalar shear rate `|gamma|`.
    ///
    /// * Below yield: returns `RIGID_MU` (effectively rigid).
    /// * Above yield: returns `mu_p + tau_yield / |gamma|`.
    pub fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        if gamma < 1e-15 {
            return RIGID_MU;
        }
        // Check whether the material has yielded.
        if self.plastic_viscosity * gamma < self.yield_stress {
            RIGID_MU
        } else {
            (self.plastic_viscosity + self.yield_stress / gamma).clamp(MU_MIN, MU_MAX)
        }
    }
}

impl NonNewtonianFluid for BinghamFluid {
    fn effective_viscosity(&self, shear_rate: f64) -> f64 {
        self.effective_viscosity(shear_rate)
    }
}

impl super::LocalViscosityModel for BinghamFluid {
    fn viscosity(&self, gamma_dot: f64) -> f64 {
        self.effective_viscosity(gamma_dot)
    }
}

// ---------------------------------------------------------------------------
// BinghamPlastic (alternative naming, wraps BinghamFluid)
// ---------------------------------------------------------------------------

/// Bingham plastic model with named fields `tau_y` and `mu_p`.
///
/// This is an alias-like struct for ergonomic API usage alongside the
/// existing `BinghamFluid`.
#[derive(Debug, Clone, Copy)]
pub struct BinghamPlastic {
    /// Yield stress.
    pub tau_y: f64,
    /// Plastic viscosity.
    pub mu_p: f64,
}

impl BinghamPlastic {
    /// Create a new Bingham plastic.
    pub fn new(tau_y: f64, mu_p: f64) -> Self {
        Self { tau_y, mu_p }
    }

    /// Effective viscosity at the given shear rate.
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        let inner = BinghamFluid::new(self.tau_y, self.mu_p);
        inner.effective_viscosity(shear_rate)
    }

    /// Shear stress: `tau_y + mu_p * |gamma|` (yielded region).
    pub fn stress(&self, shear_rate: f64) -> f64 {
        self.tau_y + self.mu_p * shear_rate.abs()
    }

    /// Returns `true` if the material is flowing (shear stress exceeds yield).
    pub fn is_flowing(&self, shear_rate: f64, tol: f64) -> bool {
        self.mu_p * shear_rate.abs() >= self.tau_y - tol
    }
}

impl NonNewtonianModel for BinghamPlastic {
    fn viscosity(&self, shear_rate: f64) -> f64 {
        BinghamPlastic::viscosity(self, shear_rate)
    }
}

// ---------------------------------------------------------------------------
// RegularizedBingham (Papanastasiou)
// ---------------------------------------------------------------------------

/// Regularized Bingham model using the Papanastasiou exponential regularization.
///
/// ```text
/// mu_eff = mu_p + tau_y * (1 - exp(-m * |gamma|)) / |gamma|
/// ```
///
/// where `m` is the regularization parameter (large m -> sharp transition).
/// This avoids the discontinuity at the yield point, making the model
/// suitable for iterative LBM solvers.
#[derive(Debug, Clone, Copy)]
pub struct RegularizedBingham {
    /// Yield stress tau_y.
    pub tau_y: f64,
    /// Plastic viscosity mu_p.
    pub mu_p: f64,
    /// Regularization parameter m (dimensionless, typically 100-10000).
    pub m: f64,
}

impl RegularizedBingham {
    /// Create a new regularized Bingham model.
    pub fn new(tau_y: f64, mu_p: f64, m: f64) -> Self {
        Self { tau_y, mu_p, m }
    }

    /// Effective viscosity at the given shear rate.
    pub fn viscosity(&self, shear_rate: f64) -> f64 {
        let gamma = shear_rate.abs();
        if gamma < 1e-15 {
            // Limit as gamma -> 0: tau_y * m + mu_p
            return (self.tau_y * self.m + self.mu_p).clamp(MU_MIN, MU_MAX);
        }
        let mu = self.mu_p + self.tau_y * (1.0 - (-self.m * gamma).exp()) / gamma;
        mu.clamp(MU_MIN, MU_MAX)
    }

    /// Shear stress at the given shear rate.
    pub fn stress(&self, shear_rate: f64) -> f64 {
        self.viscosity(shear_rate) * shear_rate.abs()
    }
}

// ---------------------------------------------------------------------------
// BinghamPipeFlow (Buckingham equation)
// ---------------------------------------------------------------------------

/// Bingham plastic pipe flow analysis using the Buckingham equation.
///
/// For pressure-driven flow in a circular pipe of radius `R` with a Bingham
/// fluid, there is a central plug region of radius `r_p = 2 * tau_y * L / dP`
/// where the fluid moves as a rigid body.
///
/// The Buckingham-Reiner equation relates flow rate Q to pressure drop dP:
///
/// ```text
/// Q = (pi R^4 dP) / (8 mu_p L) * [1 - 4/3 * (tau_y R / dP L) + 1/3 * (tau_y R / dP L)^4]
/// ```
///
/// (for the yielded regime where dP / (2L) > tau_y / R).
#[derive(Debug, Clone, Copy)]
pub struct BinghamPipeFlow {
    /// Yield stress tau_y (Pa).
    pub tau_y: f64,
    /// Plastic viscosity mu_p (Pa*s).
    pub mu_p: f64,
    /// Pipe radius R (m).
    pub radius: f64,
    /// Pipe length L (m).
    pub length: f64,
}

impl BinghamPipeFlow {
    /// Create a new Bingham pipe flow instance.
    pub fn new(tau_y: f64, mu_p: f64, radius: f64, length: f64) -> Self {
        Self {
            tau_y,
            mu_p,
            radius,
            length,
        }
    }

    /// Wall shear stress: `tau_wall = R * dP / (2L)`.
    pub fn wall_shear_stress(&self, delta_p: f64) -> f64 {
        self.radius * delta_p / (2.0 * self.length)
    }

    /// Plug-flow radius: `r_p = tau_y / (dP / (2L)) = 2 * tau_y * L / dP`.
    ///
    /// Returns `None` if the wall shear stress is below the yield stress
    /// (no flow).
    pub fn plug_radius(&self, delta_p: f64) -> Option<f64> {
        let tau_wall = self.wall_shear_stress(delta_p);
        if tau_wall <= self.tau_y {
            return None; // no flow
        }
        Some(self.tau_y * self.radius / tau_wall)
    }

    /// Velocity profile `u(r)` at radial position `r` (m from centre).
    ///
    /// Returns `None` if `r > R` or if the wall shear stress is below yield.
    ///
    /// Inside the plug (`r < r_p`): `u = u_plug`.
    /// Outside the plug (`r_p < r <= R`):
    ///
    /// ```text
    /// u(r) = 1/(2 mu_p) * (tau_wall - tau_y) * (R - r) - tau_y/(2 mu_p) * (r - r_p)
    /// ```
    pub fn velocity_profile(&self, r: f64, delta_p: f64) -> Option<f64> {
        if r > self.radius {
            return None;
        }
        let tau_wall = self.wall_shear_stress(delta_p);
        if tau_wall <= self.tau_y {
            return Some(0.0); // no flow
        }
        let r_p = self.tau_y * self.radius / tau_wall;
        if r <= r_p {
            // plug region -- same velocity as plug centre
            let u_plug = (delta_p / (4.0 * self.mu_p * self.length)) * (self.radius - r_p).powi(2);
            Some(u_plug)
        } else {
            // yielded annulus
            let factor = delta_p / (4.0 * self.mu_p * self.length);
            let u = factor * ((self.radius.powi(2) - r.powi(2)) - 2.0 * r_p * (self.radius - r));
            Some(u.max(0.0))
        }
    }

    /// Volume-averaged flow velocity (mean velocity in the pipe).
    ///
    /// Obtained by integrating the velocity profile over the cross-section.
    pub fn mean_velocity(&self, delta_p: f64) -> f64 {
        let tau_wall = self.wall_shear_stress(delta_p);
        if tau_wall <= self.tau_y {
            return 0.0;
        }
        let xi = self.tau_y / tau_wall; // dimensionless yield point
        // Buckingham-Reiner formula: u_mean = (R^2 * (dP/2muL)) / 8 * f(xi)
        // f(xi) = 1 - 4/3 * xi + 1/3 * xi^4
        let f_xi = 1.0 - (4.0 / 3.0) * xi + (1.0 / 3.0) * xi.powi(4);
        let pressure_factor = delta_p * self.radius.powi(2) / (8.0 * self.mu_p * self.length);
        pressure_factor * f_xi
    }

    /// Volume flow rate Q from the Buckingham-Reiner equation.
    ///
    /// `Q = pi R^2 * U_mean`
    pub fn flow_rate(&self, delta_p: f64) -> f64 {
        std::f64::consts::PI * self.radius.powi(2) * self.mean_velocity(delta_p)
    }

    /// Bingham number `Bn = tau_y * L / (mu_p * U_ref)`.
    ///
    /// Ratio of yield stress to viscous stress; Bn >> 1 means yield stress dominates.
    pub fn bingham_number(&self, u_ref: f64) -> f64 {
        if u_ref < 1e-30 || self.mu_p < 1e-30 {
            return f64::INFINITY;
        }
        self.tau_y * self.length / (self.mu_p * u_ref)
    }

    /// Check whether flow is possible for a given pressure drop.
    pub fn is_flowing(&self, delta_p: f64) -> bool {
        self.wall_shear_stress(delta_p) > self.tau_y
    }

    /// Friction factor `f = 16/Re_eff` for the Bingham plastic in pipe flow
    /// (generalized Fanning friction factor).
    pub fn fanning_friction_factor(&self, u_mean: f64, rho: f64) -> f64 {
        if u_mean < 1e-30 || rho < 1e-30 {
            return f64::INFINITY;
        }
        let re = rho * u_mean * 2.0 * self.radius / self.mu_p;
        16.0 / re
    }
}

// ---------------------------------------------------------------------------
// Viscoplastic yield-surface check
// ---------------------------------------------------------------------------

/// Check whether a point in a 2-D flow is in the yielded or un-yielded region
/// for a viscoplastic fluid with yield stress `tau_y`.
///
/// The von-Mises yield criterion in 2-D is:
///
/// ```text
/// sigma_vm = sqrt(s_xx^2 + s_yy^2 + 2 * s_xy^2)
/// ```
///
/// where `s_ij` are the deviatoric stress components.  If `sigma_vm >= tau_y`
/// the material has yielded; otherwise it behaves as a rigid solid.
///
/// # Returns
/// `true` if the material is yielded (flowing), `false` if rigid.
pub fn is_yielded_2d(s_xx: f64, s_yy: f64, s_xy: f64, tau_y: f64) -> bool {
    let vm = (s_xx * s_xx + s_yy * s_yy + 2.0 * s_xy * s_xy).sqrt();
    vm >= tau_y
}

/// Compute the von-Mises stress invariant for a symmetric 2-D stress tensor.
///
/// `sigma_vm = sqrt(s_xx^2 + s_yy^2 + 2*s_xy^2)`
pub fn von_mises_stress_2d(s_xx: f64, s_yy: f64, s_xy: f64) -> f64 {
    (s_xx * s_xx + s_yy * s_yy + 2.0 * s_xy * s_xy).sqrt()
}

/// Apply Bingham yield-surface correction to a stress-tensor component.
///
/// In the Papanastasiou / Augmented Lagrangian framework the deviatoric stress
/// is returned as-is if `|s| >= tau_y`, and clamped to `|s|` direction scaled
/// to `tau_y` otherwise (the material does not flow below yield).
///
/// This scalar version applies to the magnitude `|s|` of the tensor, returning
/// the corrected magnitude.
pub fn bingham_yield_correction(stress_magnitude: f64, tau_y: f64) -> f64 {
    if stress_magnitude >= tau_y {
        stress_magnitude
    } else {
        0.0 // un-yielded: return zero effective stress
    }
}
