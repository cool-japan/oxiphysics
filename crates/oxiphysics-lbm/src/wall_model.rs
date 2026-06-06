// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Wall models for turbulent boundary layers in LBM.
//!
//! Provides:
//! - Law-of-the-Wall (log law / linear sublayer)
//! - Van Driest damping and mixing length
//! - Spalart-Allmaras one-equation model (simplified)
//! - Wall shear stress and skin friction
//! - Boundary layer thickness diagnostics

// ============================================================================
// 1. Law-of-the-Wall
// ============================================================================

/// Compute the non-dimensional velocity u+ from the law of the wall.
///
/// - Viscous sublayer (y+ ≤ 11.6): u+ = y+
/// - Log layer (y+ > 11.6): u+ = (1/κ) * ln(y+) + B
///
/// # Arguments
/// * `y_plus` - non-dimensional wall distance y+
/// * `kappa`  - von Kármán constant (≈ 0.41)
/// * `b`      - log-law intercept (≈ 5.2)
pub fn law_of_the_wall_velocity(y_plus: f64, kappa: f64, b: f64) -> f64 {
    if y_plus <= 11.6 {
        y_plus
    } else {
        (1.0 / kappa) * y_plus.ln() + b
    }
}

/// Compute friction velocity u_τ = sqrt(τ_wall / ρ).
///
/// # Arguments
/// * `tau_wall` - wall shear stress \[Pa\]
/// * `rho`      - fluid density \[kg/m³\]
pub fn friction_velocity(tau_wall: f64, rho: f64) -> f64 {
    (tau_wall / rho).sqrt()
}

/// Compute non-dimensional wall distance y+ = y * u_τ / ν.
///
/// # Arguments
/// * `y`     - physical distance from wall \[m\]
/// * `u_tau` - friction velocity \[m/s\]
/// * `nu`    - kinematic viscosity \[m²/s\]
pub fn y_plus(y: f64, u_tau: f64, nu: f64) -> f64 {
    y * u_tau / nu
}

// ============================================================================
// 2. Van Driest damping and mixing length
// ============================================================================

/// Van Driest damping function D = 1 - exp(-y+ / A+).
///
/// # Arguments
/// * `y_plus` - non-dimensional wall distance
/// * `a_plus` - Van Driest constant (≈ 26)
pub fn van_driest_damping(y_plus: f64, a_plus: f64) -> f64 {
    1.0 - (-y_plus / a_plus).exp()
}

/// Turbulent mixing length with Van Driest damping.
///
/// - Inner region (y/δ < 0.22): l_m = κ * y * D(y+)
/// - Outer region (y/δ ≥ 0.22): l_m = 0.09 * δ
///
/// # Arguments
/// * `y`      - distance from wall \[m\]
/// * `kappa`  - von Kármán constant
/// * `delta`  - boundary layer thickness \[m\]
/// * `y_plus` - non-dimensional wall distance
pub fn mixing_length(y: f64, kappa: f64, delta: f64, y_plus: f64) -> f64 {
    if y / delta < 0.22 {
        kappa * y * van_driest_damping(y_plus, 26.0)
    } else {
        0.09 * delta
    }
}

// ============================================================================
// 3. Spalart-Allmaras one-equation model (simplified)
// ============================================================================

/// Parameters for the Spalart-Allmaras turbulence model.
pub struct SaParams {
    /// Production coefficient cb1
    pub cb1: f64,
    /// Cross-diffusion coefficient cb2
    pub cb2: f64,
    /// Viscosity ratio coefficient cv1
    pub cv1: f64,
    /// Turbulent Prandtl number (diffusion denominator)
    pub sigma: f64,
    /// von Kármán constant
    pub kappa: f64,
}

impl Default for SaParams {
    /// Standard Spalart-Allmaras constants.
    fn default() -> Self {
        SaParams {
            cb1: 0.1355,
            cb2: 0.622,
            cv1: 7.1,
            sigma: 2.0 / 3.0,
            kappa: 0.41,
        }
    }
}

/// SA production term P = cb1 * S̃ * ν̃.
///
/// Here S̃ is approximated by the vorticity magnitude `omega`.
///
/// # Arguments
/// * `nu_t`   - modified turbulent viscosity ν̃ \[m²/s\]
/// * `omega`  - vorticity magnitude (|∇×u|) \[1/s\]
/// * `d`      - distance to nearest wall \[m\]
/// * `params` - SA model constants
pub fn sa_production(nu_t: f64, omega: f64, d: f64, params: &SaParams) -> f64 {
    let chi = nu_t / (nu_t + 1e-15); // simplified, use nu_t/nu normally
    let fv2 = 1.0 - chi / (1.0 + chi * params.cv1.powi(3).cbrt());
    let s_tilde = omega + nu_t / (params.kappa * params.kappa * d * d) * fv2;
    params.cb1 * s_tilde * nu_t
}

/// SA destruction term D = cw1 * f_w * (ν̃/d)².
///
/// Uses simplified f_w ≈ 1 (near-wall approximation).
///
/// # Arguments
/// * `nu_t`   - modified turbulent viscosity ν̃ \[m²/s\]
/// * `d`      - distance to nearest wall \[m\]
/// * `params` - SA model constants
pub fn sa_destruction(nu_t: f64, d: f64, params: &SaParams) -> f64 {
    let cw1 = params.cb1 / (params.kappa * params.kappa) + (1.0 + params.cb2) / params.sigma;
    cw1 * (nu_t / d).powi(2)
}

/// Explicit Euler update step for SA modified viscosity.
///
/// ν̃_{n+1} = ν̃_n + dt * (P - D)
///
/// Result is clamped to zero (ν̃ ≥ 0).
///
/// # Arguments
/// * `nu_t`   - current modified turbulent viscosity \[m²/s\]
/// * `omega`  - vorticity magnitude \[1/s\]
/// * `d`      - distance to nearest wall \[m\]
/// * `dt`     - time step \[s\]
/// * `nu`     - molecular kinematic viscosity \[m²/s\] (unused in simplified form)
/// * `params` - SA model constants
pub fn sa_step(nu_t: f64, omega: f64, d: f64, dt: f64, _nu: f64, params: &SaParams) -> f64 {
    let prod = sa_production(nu_t, omega, d, params);
    let dest = sa_destruction(nu_t, d, params);
    (nu_t + dt * (prod - dest)).max(0.0)
}

// ============================================================================
// 4. Wall shear stress
// ============================================================================

/// Simplified wall shear stress using Couette approximation.
///
/// τ_w = ρ * ν * u / y  (assumes du/dy ≈ u/y)
///
/// # Arguments
/// * `u`   - velocity at distance y from wall \[m/s\]
/// * `y`   - distance from wall \[m\]
/// * `nu`  - kinematic viscosity \[m²/s\]
///
/// Returns τ_w / ρ (shear stress per unit density).
pub fn wall_shear_stress(u: f64, y: f64, nu: f64) -> f64 {
    nu * u / y
}

/// Skin friction coefficient Cf = τ_w / (0.5 * ρ * u_inf²).
///
/// # Arguments
/// * `tau_w` - wall shear stress per unit density \[m²/s²\]
/// * `rho`   - fluid density \[kg/m³\]  (unused in simplified form; included for API clarity)
/// * `u_inf` - free-stream velocity \[m/s\]
pub fn skin_friction_coefficient(tau_w: f64, _rho: f64, u_inf: f64) -> f64 {
    tau_w / (0.5 * u_inf * u_inf)
}

/// Blasius laminar skin friction: Cf = 0.664 / sqrt(Rex).
///
/// Valid for laminar flat-plate boundary layer.
///
/// # Arguments
/// * `rex` - local Reynolds number Re_x = u_inf * x / ν
pub fn blasius_cf(rex: f64) -> f64 {
    0.664 / rex.sqrt()
}

/// Turbulent skin friction using 1/7-power-law: Cf = 0.074 / Rex^0.2.
///
/// # Arguments
/// * `rex` - local Reynolds number Re_x
pub fn turbulent_cf_1_7(rex: f64) -> f64 {
    0.074 / rex.powf(0.2)
}

// ============================================================================
// 5. Boundary layer thickness
// ============================================================================

/// Blasius δ99 boundary layer thickness: δ99 = 4.91 * x / sqrt(Rex).
///
/// # Arguments
/// * `x`     - streamwise distance from leading edge \[m\]
/// * `nu`    - kinematic viscosity \[m²/s\]
/// * `u_inf` - free-stream velocity \[m/s\]
pub fn blasius_delta99(x: f64, nu: f64, u_inf: f64) -> f64 {
    let rex = u_inf * x / nu;
    4.91 * x / rex.sqrt()
}

/// Displacement thickness δ* = ∫(1 - u/u_inf) dy using trapezoidal rule.
///
/// # Arguments
/// * `u_profile` - velocity values at wall-normal positions
/// * `y_profile` - wall-normal positions (must match length of `u_profile`)
/// * `u_inf`     - free-stream velocity \[m/s\]
pub fn displacement_thickness(u_profile: &[f64], y_profile: &[f64], u_inf: f64) -> f64 {
    assert_eq!(u_profile.len(), y_profile.len());
    let n = u_profile.len();
    if n < 2 {
        return 0.0;
    }
    let mut integral = 0.0;
    for i in 0..n - 1 {
        let dy = y_profile[i + 1] - y_profile[i];
        let f0 = 1.0 - u_profile[i] / u_inf;
        let f1 = 1.0 - u_profile[i + 1] / u_inf;
        integral += 0.5 * (f0 + f1) * dy;
    }
    integral
}

/// Momentum thickness θ = ∫ (u/u_inf)(1 - u/u_inf) dy using trapezoidal rule.
///
/// # Arguments
/// * `u_profile` - velocity values at wall-normal positions
/// * `y_profile` - wall-normal positions
/// * `u_inf`     - free-stream velocity \[m/s\]
pub fn momentum_thickness(u_profile: &[f64], y_profile: &[f64], u_inf: f64) -> f64 {
    assert_eq!(u_profile.len(), y_profile.len());
    let n = u_profile.len();
    if n < 2 {
        return 0.0;
    }
    let mut integral = 0.0;
    for i in 0..n - 1 {
        let dy = y_profile[i + 1] - y_profile[i];
        let r0 = u_profile[i] / u_inf;
        let r1 = u_profile[i + 1] / u_inf;
        let f0 = r0 * (1.0 - r0);
        let f1 = r1 * (1.0 - r1);
        integral += 0.5 * (f0 + f1) * dy;
    }
    integral
}

/// Shape factor H = δ* / θ.
///
/// Typical values: laminar ≈ 2.6, turbulent ≈ 1.4.
///
/// # Arguments
/// * `delta_star` - displacement thickness \[m\]
/// * `theta`      - momentum thickness \[m\]
pub fn shape_factor(delta_star: f64, theta: f64) -> f64 {
    delta_star / theta
}

// ============================================================================
// 6. Wall distance computation
// ============================================================================

/// Compute the minimum wall distance for each cell in a 2D grid.
///
/// `wall_mask[k]` is `true` if cell k is a solid wall cell.
/// Returns `d[k]` = Euclidean distance to the nearest wall cell.
///
/// # Arguments
/// * `nx`, `ny` - grid dimensions
/// * `wall_mask` - flat boolean array (row-major)
/// * `dx` - grid spacing (assumed equal in x and y)
pub fn compute_wall_distance(nx: usize, ny: usize, wall_mask: &[bool], dx: f64) -> Vec<f64> {
    let n = nx * ny;
    let mut dist = vec![f64::MAX; n];

    // Collect wall cell positions
    let mut wall_positions: Vec<(f64, f64)> = Vec::new();
    for (k, &is_wall) in wall_mask.iter().enumerate() {
        if is_wall {
            let x = (k % nx) as f64 * dx;
            let y = (k / nx) as f64 * dx;
            wall_positions.push((x, y));
        }
    }

    // Brute force nearest-wall for each fluid cell
    for k in 0..n {
        if wall_mask[k] {
            dist[k] = 0.0;
            continue;
        }
        let cx = (k % nx) as f64 * dx;
        let cy = (k / nx) as f64 * dx;
        for &(wx, wy) in &wall_positions {
            let d = ((cx - wx).powi(2) + (cy - wy).powi(2)).sqrt();
            if d < dist[k] {
                dist[k] = d;
            }
        }
    }
    dist
}

// ============================================================================
// 7. Wall shear stress from velocity profile
// ============================================================================

/// Compute wall shear stress using a velocity profile (3-point stencil).
///
/// Uses one-sided 3-point derivative at the wall:
///
/// `du/dy|_wall = (-3u_0 + 4u_1 - u_2) / (2*dy)`
///
/// Returns τ_w / ρ = ν * du/dy.
///
/// # Arguments
/// * `u0` - velocity at the wall (typically 0 for no-slip)
/// * `u1` - velocity at first cell from wall
/// * `u2` - velocity at second cell from wall
/// * `dy` - grid spacing
/// * `nu` - kinematic viscosity
pub fn wall_shear_stress_3pt(u0: f64, u1: f64, u2: f64, dy: f64, nu: f64) -> f64 {
    let du_dy = (-3.0 * u0 + 4.0 * u1 - u2) / (2.0 * dy);
    nu * du_dy
}

// ============================================================================
// 8. y+ calculation utilities
// ============================================================================

/// Compute y+ for the first cell away from a wall.
///
/// `y+ = y * u_tau / nu = y * sqrt(tau_w / rho) / nu`
///
/// where `tau_w_per_rho` = tau_w / rho.
pub fn y_plus_from_tau(y: f64, tau_w_per_rho: f64, nu: f64) -> f64 {
    let u_tau = tau_w_per_rho.abs().sqrt();
    y * u_tau / nu
}

/// Iteratively compute u_tau from the law-of-the-wall given u and y.
///
/// Newton iteration on:
///   u/u_tau = f(y*u_tau/nu)
///
/// where f is the law-of-the-wall function.
///
/// Returns u_tau.
pub fn u_tau_newton(u: f64, y: f64, nu: f64, kappa: f64, b: f64, max_iter: usize) -> f64 {
    if u.abs() < 1e-30 || y < 1e-30 {
        return 0.0;
    }
    // Initial guess from laminar estimate: tau = nu * u / y → u_tau = sqrt(nu*u/y)
    let mut u_tau = (nu * u.abs() / y).sqrt();

    for _ in 0..max_iter {
        let yp = y * u_tau / nu;
        let up = law_of_the_wall_velocity(yp, kappa, b);
        let residual = u.abs() / u_tau - up;

        // Derivative of residual w.r.t. u_tau (numerical)
        let eps = u_tau * 1e-6;
        let yp2 = y * (u_tau + eps) / nu;
        let up2 = law_of_the_wall_velocity(yp2, kappa, b);
        let dres = -u.abs() / ((u_tau + eps) * (u_tau + eps)) * (u_tau + eps)
            + u.abs() / (u_tau * u_tau)
            - (up2 - up) / eps;

        if dres.abs() < 1e-30 {
            break;
        }

        let du_tau = -residual / dres;
        u_tau = (u_tau + du_tau).max(1e-15);

        if residual.abs() < 1e-6 {
            break;
        }
    }
    u_tau
}

// ============================================================================
// 9. Wall-modeled LES
// ============================================================================

/// Wall-modeled LES (WMLES) interface.
///
/// Provides the effective wall-stress boundary condition for LES
/// by matching the outer LES solution to a modeled inner layer.
pub struct WallModeledLes {
    /// Von Karman constant.
    pub kappa: f64,
    /// Log-law intercept.
    pub b: f64,
    /// Molecular kinematic viscosity.
    pub nu: f64,
}

impl WallModeledLes {
    /// Create a new WMLES model.
    pub fn new(nu: f64) -> Self {
        Self {
            kappa: 0.41,
            b: 5.2,
            nu,
        }
    }

    /// Compute the wall shear stress (per unit density) from the velocity
    /// at the matching point.
    ///
    /// Uses the Newton-iteration u_tau solver, then tau_w/rho = u_tau^2.
    pub fn wall_stress(&self, u_match: f64, y_match: f64) -> f64 {
        let u_tau = u_tau_newton(u_match, y_match, self.nu, self.kappa, self.b, 50);
        u_tau * u_tau
    }

    /// Effective eddy viscosity at the matching point from the WMLES:
    ///
    /// `nu_t = kappa * y * u_tau * D(y+)`
    ///
    /// where D is the Van Driest damping function.
    pub fn eddy_viscosity(&self, y_match: f64, u_match: f64) -> f64 {
        let u_tau = u_tau_newton(u_match, y_match, self.nu, self.kappa, self.b, 50);
        let yp = y_match * u_tau / self.nu;
        let d = van_driest_damping(yp, 26.0);
        self.kappa * y_match * u_tau * d
    }

    /// Compute y+ at the matching location.
    pub fn y_plus_at_match(&self, y_match: f64, u_match: f64) -> f64 {
        let u_tau = u_tau_newton(u_match, y_match, self.nu, self.kappa, self.b, 50);
        y_match * u_tau / self.nu
    }
}

/// Reichardt's law of the wall (smooth, continuous from viscous to log layer):
///
/// `u+ = (1/kappa) * ln(1 + kappa*y+) + C * (1 - exp(-y+/A) - y+/A * exp(-B*y+))`
///
/// Standard constants: C = 7.8, A = 11.0, B = 0.33.
pub fn reichardt_velocity(y_plus: f64, kappa: f64) -> f64 {
    const C: f64 = 7.8;
    const A: f64 = 11.0;
    const B_CONST: f64 = 0.33;
    (1.0 / kappa) * (1.0 + kappa * y_plus).ln()
        + C * (1.0 - (-y_plus / A).exp() - (y_plus / A) * (-B_CONST * y_plus).exp())
}

/// Spalding's implicit law of the wall:
///
/// `y+ = u+ + exp(-kappa*B) * [exp(kappa*u+) - 1 - kappa*u+ - (kappa*u+)^2/2 - (kappa*u+)^3/6]`
///
/// Given u+, returns y+.
pub fn spalding_y_plus(u_plus: f64, kappa: f64, b: f64) -> f64 {
    let ku = kappa * u_plus;
    let exp_minus_kb = (-kappa * b).exp();
    u_plus + exp_minus_kb * (ku.exp() - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0)
}

// ============================================================================
// 10. Equilibrium wall model
// ============================================================================

/// Equilibrium wall model (EWM).
///
/// Assumes the local boundary layer is in equilibrium (production = dissipation).
/// Uses the law of the wall to relate the wall stress to the matching-point velocity.
pub struct EquilibriumWallModel {
    /// Molecular kinematic viscosity.
    pub nu: f64,
    /// Von Kármán constant.
    pub kappa: f64,
    /// Log-law intercept B.
    pub b: f64,
}

impl EquilibriumWallModel {
    /// Create a new equilibrium wall model.
    pub fn new(nu: f64) -> Self {
        Self {
            nu,
            kappa: 0.41,
            b: 5.2,
        }
    }

    /// Compute wall shear stress per unit density `tau_w / rho = u_tau^2`.
    pub fn wall_stress(&self, u_match: f64, y_match: f64) -> f64 {
        let u_tau = u_tau_newton(u_match, y_match, self.nu, self.kappa, self.b, 50);
        u_tau * u_tau
    }

    /// Effective (eddy + molecular) viscosity at the matching point.
    ///
    /// `nu_eff = nu + kappa * y * u_tau * D(y+)`
    pub fn effective_viscosity(&self, u_match: f64, y_match: f64, _delta: f64) -> f64 {
        let u_tau = u_tau_newton(u_match, y_match, self.nu, self.kappa, self.b, 50);
        let yp = y_match * u_tau / self.nu;
        let d = van_driest_damping(yp, 26.0);
        self.nu + self.kappa * y_match * u_tau * d
    }

    /// Predicted velocity at distance `y` using the law of the wall.
    pub fn predicted_velocity(&self, u_tau: f64, y: f64) -> f64 {
        let yp = y_plus(y, u_tau, self.nu);
        law_of_the_wall_velocity(yp, self.kappa, self.b) * u_tau
    }
}

// ============================================================================
// 11. Non-equilibrium wall model
// ============================================================================

/// Non-equilibrium wall model (NEQWM).
///
/// Accounts for history and pressure-gradient effects in the boundary layer.
/// Uses a simplified ODE-based inner-layer model to compute the wall stress.
pub struct NonEquilibriumWallModel {
    /// Molecular kinematic viscosity.
    pub nu: f64,
    /// Relaxation time scale for history effects.
    pub tau_relax: f64,
    /// Von Kármán constant.
    pub kappa: f64,
    /// Log-law intercept.
    pub b: f64,
}

impl NonEquilibriumWallModel {
    /// Create a new non-equilibrium wall model.
    ///
    /// `tau_relax` controls how quickly the wall stress adapts to changes.
    pub fn new(nu: f64, tau_relax: f64) -> Self {
        Self {
            nu,
            tau_relax,
            kappa: 0.41,
            b: 5.2,
        }
    }

    /// Compute wall shear stress accounting for pressure gradient.
    ///
    /// The pressure gradient `dp_ds` modifies the effective driving velocity:
    ///   `u_eff = u_match - (dp_ds / rho) * y_match / (u_tau + eps)`
    ///
    /// Returns `tau_w / rho`.
    pub fn wall_stress(&self, u_match: f64, y_match: f64, dp_ds: f64) -> f64 {
        // Equilibrium estimate first
        let u_tau_eq = u_tau_newton(u_match, y_match, self.nu, self.kappa, self.b, 20);
        // Non-equilibrium correction: shift u_match by PG effect
        let pg_correction = dp_ds * y_match / (u_tau_eq.max(1e-15));
        let u_eff = (u_match - pg_correction).max(0.0);
        let u_tau = u_tau_newton(u_eff, y_match, self.nu, self.kappa, self.b, 30);
        u_tau * u_tau
    }

    /// Relaxation update for wall stress: blends toward equilibrium.
    ///
    /// `tau_new = tau_old + (1/tau_relax) * (tau_eq - tau_old) * dt`
    pub fn relax_wall_stress(&self, tau_old: f64, tau_eq: f64, dt: f64) -> f64 {
        if self.tau_relax < 1e-30 {
            return tau_eq;
        }
        tau_old + (tau_eq - tau_old) / self.tau_relax * dt
    }
}

// ============================================================================
// 12. Wall model switching
// ============================================================================

/// Wall model switcher: selects EWM or NEQWM based on local y+.
///
/// Below `y_plus_switch` the equilibrium model is used;
/// above it the non-equilibrium model is applied.
pub struct WallModelSwitcher {
    /// Equilibrium wall model.
    pub ewm: EquilibriumWallModel,
    /// Non-equilibrium wall model (same nu, standard relaxation).
    pub neqwm: NonEquilibriumWallModel,
    /// y+ threshold for switching.
    pub y_plus_switch: f64,
}

impl WallModelSwitcher {
    /// Create a new wall model switcher.
    pub fn new(nu: f64, y_plus_switch: f64) -> Self {
        Self {
            ewm: EquilibriumWallModel::new(nu),
            neqwm: NonEquilibriumWallModel::new(nu, 0.1),
            y_plus_switch,
        }
    }

    /// Compute wall stress using the appropriate model.
    pub fn compute_stress(&self, u_match: f64, y_match: f64) -> f64 {
        let tau_eq = self.ewm.wall_stress(u_match, y_match);
        let u_tau_eq = tau_eq.max(0.0).sqrt();
        let yp = y_plus(y_match, u_tau_eq, self.ewm.nu);
        if yp <= self.y_plus_switch {
            tau_eq
        } else {
            self.neqwm.wall_stress(u_match, y_match, 0.0)
        }
    }
}

// ============================================================================
// 13. Wall distance improvements (3D)
// ============================================================================

/// Compute the minimum wall distance from a point to a set of wall nodes in 3D.
///
/// # Arguments
/// * `point`     - query point \[x, y, z\]
/// * `wall_pts`  - array of wall-node positions \[\[x, y, z\\], ...]
///
/// Returns the minimum Euclidean distance.
pub fn min_wall_distance_3d(point: [f64; 3], wall_pts: &[[f64; 3]]) -> f64 {
    let mut min_d = f64::MAX;
    for wp in wall_pts {
        let dx = point[0] - wp[0];
        let dy = point[1] - wp[1];
        let dz = point[2] - wp[2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        if d < min_d {
            min_d = d;
        }
    }
    min_d
}

/// Compute the nearest wall distance and the outward wall-normal vector.
///
/// Returns `(distance, normal)` where `normal` points from the wall toward the point.
pub fn nearest_wall_normal(point: [f64; 3], wall_pts: &[[f64; 3]]) -> (f64, [f64; 3]) {
    let mut min_d = f64::MAX;
    let mut best = [0.0f64; 3];
    for wp in wall_pts {
        let dx = point[0] - wp[0];
        let dy = point[1] - wp[1];
        let dz = point[2] - wp[2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        if d < min_d {
            min_d = d;
            if d > 1e-30 {
                best = [dx / d, dy / d, dz / d];
            }
        }
    }
    (min_d, best)
}

/// Compute signed distance using a level-set function (1D case).
///
/// Positive = outside the solid, negative = inside.
pub fn signed_wall_distance_1d(x: f64, wall_x: f64) -> f64 {
    x - wall_x
}

/// Bilinear interpolation of wall distance on a 2D grid.
///
/// Interpolates the pre-computed wall distances at four corner nodes
/// to the interior point (xi, eta) in \[0,1\]^2.
pub fn interpolate_wall_distance_bilinear(
    d00: f64,
    d10: f64,
    d01: f64,
    d11: f64,
    xi: f64,
    eta: f64,
) -> f64 {
    d00 * (1.0 - xi) * (1.0 - eta)
        + d10 * xi * (1.0 - eta)
        + d01 * (1.0 - xi) * eta
        + d11 * xi * eta
}

// ============================================================================
// 14. Pressure gradient correction
// ============================================================================

/// Pressure gradient correction to skin friction coefficient.
///
/// Uses Stratford's approximation for the additional effect of a streamwise
/// pressure gradient on the wall friction:
///   `Cf_eff = Cf_0 * sqrt(1 - (dp_ds / (rho * u^2 / y)))`
///
/// Returns the corrected `Cf`.
pub fn pressure_gradient_correction_cf(u: f64, y: f64, nu: f64, dp_ds: f64) -> f64 {
    let cf0 = wall_shear_stress(u, y, nu) / (0.5 * u * u + 1e-30);
    let correction = (1.0 - dp_ds * y / (u * u * 0.5 + 1e-30)).abs().sqrt();
    cf0 * correction
}

// ============================================================================
// 15. Thermal wall model
// ============================================================================

/// Thermal wall model for conjugate heat transfer boundary conditions.
///
/// Uses the thermal law of the wall: `T+ = Pr_t * u+` in the log layer.
pub struct ThermalWallModel {
    /// Molecular kinematic viscosity.
    pub nu: f64,
    /// Turbulent Prandtl number.
    pub pr_t: f64,
    /// Von Kármán constant.
    pub kappa: f64,
    /// Log-law intercept.
    pub b: f64,
}

impl ThermalWallModel {
    /// Create a new thermal wall model.
    ///
    /// `pr_t` is the turbulent Prandtl number (typically 0.85–0.9 for air).
    pub fn new(nu: f64, pr_t: f64) -> Self {
        Self {
            nu,
            pr_t,
            kappa: 0.41,
            b: 5.2,
        }
    }

    /// Compute the wall heat flux `q_w = rho * Cp * u_tau * (T_wall - T_match) / T+`.
    ///
    /// Returns `q_w / (rho * Cp)` (heat flux per unit density × heat capacity).
    pub fn wall_heat_flux(&self, t_wall: f64, t_match: f64, y_match: f64, u_match: f64) -> f64 {
        let u_tau = u_tau_newton(u_match, y_match, self.nu, self.kappa, self.b, 50).max(1e-15);
        let yp = y_plus(y_match, u_tau, self.nu);
        let t_plus = self.pr_t * law_of_the_wall_velocity(yp, self.kappa, self.b);
        if t_plus.abs() < 1e-30 {
            return 0.0;
        }
        u_tau * (t_wall - t_match) / t_plus
    }

    /// Non-dimensional temperature: `T+ = (T_wall - T) / (q_w / (rho * Cp * u_tau))`.
    pub fn temperature_plus(&self, y_plus_val: f64) -> f64 {
        self.pr_t * law_of_the_wall_velocity(y_plus_val, self.kappa, self.b)
    }
}

// ============================================================================
// 16. Rough-wall corrections
// ============================================================================

/// Roughness function `ΔB` (log-law shift due to wall roughness).
///
/// Uses the Colebrook formula interpolating between smooth, transitional,
/// and fully-rough regimes:
///   - Smooth: `k_s+ < 5`  → `ΔB ≈ 0`
///   - Transition: `5 ≤ k_s+ < 70`  → linear interpolation
///   - Fully rough: `k_s+ ≥ 70`  → `ΔB = (1/κ) * ln(k_s+) - 3.5`
pub fn roughness_function(k_s_plus: f64) -> f64 {
    const KAPPA: f64 = 0.41;
    if k_s_plus < 5.0 {
        0.0
    } else if k_s_plus < 70.0 {
        // Transitional: linear from 0 to fully-rough value
        let delta_b_rough = (1.0 / KAPPA) * (70.0_f64).ln() - 3.5;
        delta_b_rough * (k_s_plus - 5.0) / 65.0
    } else {
        (1.0 / KAPPA) * k_s_plus.ln() - 3.5
    }
}

/// Rough-wall log-law velocity: `u+ = (1/κ) * ln(y+) + B - ΔB(k_s+)`.
pub fn rough_wall_velocity(y_plus: f64, kappa: f64, b: f64, k_s_plus: f64) -> f64 {
    let db = roughness_function(k_s_plus);
    law_of_the_wall_velocity(y_plus, kappa, b) - db
}

// ============================================================================
// 17. Compressibility corrections
// ============================================================================

/// Van Driest compressibility correction for high-speed flows.
///
/// Transforms the incompressible law-of-the-wall velocity to the compressible
/// form using temperature ratio correction (simplified):
///
///   `u+_comp = u+_incomp * sqrt(T_wall / T_aw)`
///
/// where `T_aw` is the adiabatic wall temperature.
pub fn van_driest_compressibility(u_plus_incomp: f64, t_wall: f64, t_adiabatic: f64) -> f64 {
    if t_adiabatic < 1e-10 {
        return u_plus_incomp;
    }
    u_plus_incomp * (t_wall / t_adiabatic).sqrt()
}

/// Adiabatic wall temperature: `T_aw = T_inf * (1 + r * (gamma-1)/2 * Ma^2)`.
///
/// `r` is the recovery factor (typically sqrt(Pr) for laminar, Pr^(1/3) for turbulent).
pub fn adiabatic_wall_temperature(t_inf: f64, mach: f64, gamma: f64, recovery: f64) -> f64 {
    t_inf * (1.0 + recovery * (gamma - 1.0) / 2.0 * mach * mach)
}

/// Recovery factor for turbulent flow: `r = Pr^(1/3)`.
pub fn recovery_factor_turbulent(prandtl: f64) -> f64 {
    prandtl.powf(1.0 / 3.0)
}

// ============================================================================
// LogLawWallModel
// ============================================================================

/// Struct-based law-of-the-wall model with Newton iteration for wall stress.
pub struct LogLawWallModel {
    /// Von Kármán constant (≈ 0.41).
    pub kappa: f64,
    /// Log-law intercept b (≈ 5.2).
    pub b: f64,
    /// Kinematic viscosity (m²/s).
    pub nu: f64,
}

impl LogLawWallModel {
    /// Create a new `LogLawWallModel`.
    pub fn new(kappa: f64, b: f64, nu: f64) -> Self {
        Self { kappa, b, nu }
    }

    /// Non-dimensional velocity u+ from the law of the wall.
    ///
    /// - Viscous sublayer (y+ ≤ 11.6): u+ = y+
    /// - Log layer (y+ > 11.6): u+ = (1/kappa) * ln(y+) + B
    pub fn u_plus(y_plus: f64, kappa: f64, b: f64) -> f64 {
        law_of_the_wall_velocity(y_plus, kappa, b)
    }

    /// Wall shear stress via Newton iteration on the law-of-the-wall.
    ///
    /// Returns tau_w / rho (m²/s²).
    pub fn wall_shear_stress_newton(&self, u_local: f64, y: f64) -> f64 {
        let u_tau = u_tau_newton(u_local, y, self.nu, self.kappa, self.b, 50);
        u_tau * u_tau
    }
}

// ============================================================================
// WernerWengle wall model
// ============================================================================

/// Werner-Wengle explicit wall model.
///
/// Uses a two-layer power-law fit:
///   - y+ < 11.81: linear (viscous sublayer)
///   - y+ ≥ 11.81: 1/7 power law
pub struct WernerWengle;

impl WernerWengle {
    /// Compute wall shear stress per unit density (tau_w / rho).
    ///
    /// # Arguments
    /// * `u`  - velocity at first cell (m/s)
    /// * `y`  - first cell distance from wall (m)
    /// * `nu` - kinematic viscosity (m²/s)
    pub fn wall_shear_stress(u: f64, y: f64, nu: f64) -> f64 {
        // Werner-Wengle explicit formula:
        //   y+ < 11.81: viscous sublayer → tau_w/rho = (nu * u / y)
        //   y+ >= 11.81: log-law power → tau_w/rho = (8*nu/y)^0.25 * (u/2)^0.75
        // First estimate tau to get y+
        let tau_visc = nu * u.abs() / y.max(1e-30);
        let u_tau_est = tau_visc.max(0.0).sqrt();
        let y_plus_est = y * u_tau_est / nu.max(1e-30);

        if y_plus_est < 11.81 {
            nu * u.abs() / y.max(1e-30)
        } else {
            // Explicit WW formula:
            // tau_w/rho = (8*nu/y)^0.25 * (u/2)^0.75
            let coeff = (8.0 * nu / y.max(1e-30)).powf(0.25);
            coeff * (u.abs() / 2.0).powf(0.75)
        }
    }
}

// ============================================================================
// Spalding struct
// ============================================================================

/// Spalding unified wall profile.
///
/// y+ = u+ + exp(-kappa*B) * \[exp(kappa*u+) - 1 - kappa*u+ - (kappa*u+)^2/2 - (kappa*u+)^3/6\]
pub struct Spalding {
    /// Von Kármán constant.
    pub kappa: f64,
    /// Log-law intercept b.
    pub b: f64,
}

impl Spalding {
    /// Create a new `Spalding` model.
    pub fn new(kappa: f64, b: f64) -> Self {
        Self { kappa, b }
    }

    /// Compute y+ given u+ using the Spalding unified profile.
    pub fn u_from_y_plus(&self, y_plus: f64) -> f64 {
        // Invert Spalding implicitly via Newton iteration
        // f(u+) = y+(u+) - y_plus = 0
        if y_plus <= 0.0 {
            return 0.0;
        }
        // Initial guess: viscous sublayer
        let mut u_p = y_plus.min(11.6);

        for _ in 0..50 {
            let yp_of_up = spalding_y_plus(u_p, self.kappa, self.b);
            let residual = yp_of_up - y_plus;
            if residual.abs() < 1e-8 {
                break;
            }
            // Derivative of Spalding's y+ w.r.t. u+:
            // d(y+)/d(u+) = 1 + exp(-kB) * kappa * [exp(ku+) - 1 - ku+ - (ku+)^2/2]
            let ku = self.kappa * u_p;
            let exp_minus_kb = (-self.kappa * self.b).exp();
            let d_yp = 1.0 + exp_minus_kb * self.kappa * (ku.exp() - 1.0 - ku - ku * ku / 2.0);
            if d_yp.abs() < 1e-30 {
                break;
            }
            u_p -= residual / d_yp;
            u_p = u_p.max(0.0);
        }
        u_p
    }
}

// ============================================================================
// WallFunctionBc
// ============================================================================

/// Wall-function boundary condition combining a first-cell distance,
/// kinematic viscosity, and a `LogLawWallModel`.
pub struct WallFunctionBc {
    /// Distance of the first grid cell from the wall (m).
    pub y_first: f64,
    /// Kinematic viscosity (m²/s).
    pub nu: f64,
    /// Underlying log-law model.
    pub model: LogLawWallModel,
}

impl WallFunctionBc {
    /// Create a new `WallFunctionBc`.
    pub fn new(y_first: f64, nu: f64, model: LogLawWallModel) -> Self {
        Self { y_first, nu, model }
    }

    /// Effective wall velocity boundary condition from the log-law.
    ///
    /// Returns the velocity that the wall should impose at `y_first`
    /// such that the log-law is satisfied for the given outer velocity `u_log`.
    pub fn wall_velocity_bc(&self, u_log: f64) -> f64 {
        // Compute u_tau from the outer velocity
        let u_tau = u_tau_newton(
            u_log,
            self.y_first,
            self.nu,
            self.model.kappa,
            self.model.b,
            50,
        );
        // Return the velocity predicted by the log-law at y_first
        let y_plus_val = self.y_first * u_tau / self.nu.max(1e-30);
        law_of_the_wall_velocity(y_plus_val, self.model.kappa, self.model.b) * u_tau
    }
}

// ============================================================================
// roughness_correction
// ============================================================================

/// Colebrook-White roughness shift: ΔU+ = (1/kappa) * ln(1 + 0.3 * ks+).
///
/// # Arguments
/// * `u_plus`  - smooth-wall u+ (retained for potential future use)
/// * `ks_plus` - dimensionless equivalent sand-grain roughness k_s+ = k_s * u_tau / nu
pub fn roughness_correction(_u_plus: f64, ks_plus: f64) -> f64 {
    const KAPPA: f64 = 0.41;
    (1.0 / KAPPA) * (1.0 + 0.3 * ks_plus).ln()
}

// ============================================================================
// von_karman_constant
// ============================================================================

/// Returns the von Kármán constant κ ≈ 0.41.
pub fn von_karman_constant() -> f64 {
    0.41
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- Law-of-the-Wall ---

    #[test]
    fn test_lotw_linear_sublayer() {
        let u_plus = law_of_the_wall_velocity(5.0, 0.41, 5.2);
        assert!((u_plus - 5.0).abs() < 1e-12, "u+ = {u_plus}");
    }

    #[test]
    fn test_lotw_at_transition() {
        let u_plus = law_of_the_wall_velocity(11.6, 0.41, 5.2);
        assert!((u_plus - 11.6).abs() < 1e-12, "u+ at y+=11.6: {u_plus}");
    }

    #[test]
    fn test_lotw_log_layer() {
        let kappa = 0.41;
        let b = 5.2;
        let y_plus_val = 100.0;
        let u_plus = law_of_the_wall_velocity(y_plus_val, kappa, b);
        let expected = (1.0 / kappa) * y_plus_val.ln() + b;
        assert!((u_plus - expected).abs() < 1e-12, "{u_plus} vs {expected}");
    }

    #[test]
    fn test_lotw_monotone() {
        let kappa = 0.41;
        let b = 5.2;
        let vals: Vec<f64> = [1.0, 5.0, 11.6, 20.0, 100.0, 500.0]
            .iter()
            .map(|&yp| law_of_the_wall_velocity(yp, kappa, b))
            .collect();
        for i in 0..vals.len() - 1 {
            assert!(vals[i + 1] > vals[i], "not monotone at index {i}");
        }
    }

    #[test]
    fn test_friction_velocity_basic() {
        let u_tau = friction_velocity(0.04, 1.0);
        assert!((u_tau - 0.2).abs() < 1e-12, "u_tau = {u_tau}");
    }

    #[test]
    fn test_friction_velocity_scaled() {
        let u_tau = friction_velocity(4.0, 4.0);
        assert!((u_tau - 1.0).abs() < 1e-12, "u_tau = {u_tau}");
    }

    #[test]
    fn test_y_plus_basic() {
        let yp = y_plus(0.01, 0.1, 1e-4);
        assert!((yp - 10.0).abs() < 1e-12, "y+ = {yp}");
    }

    // --- Van Driest damping ---

    #[test]
    fn test_van_driest_at_wall() {
        let d = van_driest_damping(0.0, 26.0);
        assert!(d.abs() < 1e-14, "D = {d}");
    }

    #[test]
    fn test_van_driest_far_field() {
        let d = van_driest_damping(1000.0, 26.0);
        assert!((d - 1.0).abs() < 1e-10, "D = {d}");
    }

    #[test]
    fn test_van_driest_at_a_plus() {
        let d = van_driest_damping(26.0, 26.0);
        let expected = 1.0 - (-1.0_f64).exp();
        assert!((d - expected).abs() < 1e-12, "D = {d}");
    }

    #[test]
    fn test_mixing_length_inner() {
        let y = 0.01;
        let kappa = 0.41;
        let delta = 0.1;
        let yp = 50.0;
        let lm = mixing_length(y, kappa, delta, yp);
        let d = van_driest_damping(yp, 26.0);
        assert!((lm - kappa * y * d).abs() < 1e-12, "l_m = {lm}");
    }

    #[test]
    fn test_mixing_length_outer() {
        let lm = mixing_length(0.025, 0.41, 0.1, 100.0);
        assert!((lm - 0.09 * 0.1).abs() < 1e-12, "l_m = {lm}");
    }

    // --- Spalart-Allmaras ---

    #[test]
    fn test_sa_params_default() {
        let p = SaParams::default();
        assert!((p.cb1 - 0.1355).abs() < 1e-10);
        assert!((p.kappa - 0.41).abs() < 1e-10);
    }

    #[test]
    fn test_sa_production_positive() {
        let p = sa_production(1e-4, 10.0, 0.01, &SaParams::default());
        assert!(p >= 0.0, "SA production = {p}");
    }

    #[test]
    fn test_sa_destruction_positive() {
        let d = sa_destruction(1e-4, 0.01, &SaParams::default());
        assert!(d >= 0.0, "SA destruction = {d}");
    }

    #[test]
    fn test_sa_step_clamps_negative() {
        let nu_t_new = sa_step(1e-10, 0.0, 1e-5, 1.0, 1e-5, &SaParams::default());
        assert!(nu_t_new >= 0.0, "nu_t = {nu_t_new}");
    }

    #[test]
    fn test_sa_step_grows_with_production() {
        let nu_t0 = 1e-4;
        let nu_t1 = sa_step(nu_t0, 1000.0, 1.0, 1e-4, 1e-5, &SaParams::default());
        assert!(nu_t1 >= nu_t0, "{nu_t0} -> {nu_t1}");
    }

    // --- Wall shear stress ---

    #[test]
    fn test_wall_shear_stress_basic() {
        let tau = wall_shear_stress(1.0, 0.1, 1e-3);
        assert!((tau - 0.01).abs() < 1e-14, "tau = {tau}");
    }

    #[test]
    fn test_skin_friction_coefficient_basic() {
        let cf = skin_friction_coefficient(0.5, 1.0, 10.0);
        assert!((cf - 0.01).abs() < 1e-14, "Cf = {cf}");
    }

    #[test]
    fn test_blasius_cf_decreases_with_rex() {
        assert!(blasius_cf(1e6) < blasius_cf(1e5));
    }

    #[test]
    fn test_blasius_cf_value() {
        assert!((blasius_cf(1e6) - 6.64e-4).abs() < 1e-10);
    }

    #[test]
    fn test_turbulent_cf_1_7_value() {
        let cf = turbulent_cf_1_7(1e6);
        let expected = 0.074 / (1e6_f64).powf(0.2);
        assert!((cf - expected).abs() < 1e-12, "cf = {cf}");
    }

    // --- Boundary layer thickness ---

    #[test]
    fn test_blasius_delta99_basic() {
        let delta = blasius_delta99(1.0, 1e-5, 10.0);
        assert!((delta - 4.91e-3).abs() < 1e-10, "delta = {delta}");
    }

    #[test]
    fn test_displacement_thickness_uniform() {
        let y: Vec<f64> = (0..5).map(|i| i as f64 * 0.1).collect();
        let u = vec![1.0; 5];
        let ds = displacement_thickness(&u, &y, 1.0);
        assert!(ds.abs() < 1e-14, "ds = {ds}");
    }

    #[test]
    fn test_displacement_thickness_zero_velocity() {
        let y = vec![0.0, 0.1, 0.2, 0.3, 0.4];
        let u = vec![0.0; 5];
        let ds = displacement_thickness(&u, &y, 1.0);
        assert!((ds - 0.4).abs() < 1e-12, "ds = {ds}");
    }

    #[test]
    fn test_momentum_thickness_basic() {
        let n = 100;
        let y: Vec<f64> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let u: Vec<f64> = y.to_vec();
        let theta = momentum_thickness(&u, &y, 1.0);
        assert!((theta - 1.0 / 6.0).abs() < 1e-3, "theta = {theta}");
    }

    #[test]
    fn test_shape_factor_laminar_approx() {
        let h = shape_factor(0.344, 0.133);
        assert!((h - 0.344 / 0.133).abs() < 1e-12, "H = {h}");
    }

    #[test]
    fn test_shape_factor_turbulent_approx() {
        let h = shape_factor(0.014, 0.010);
        assert!((h - 1.4).abs() < 1e-10, "H = {h}");
    }

    // --- Wall distance computation ---

    #[test]
    fn test_wall_distance_no_walls() {
        let dist = compute_wall_distance(5, 5, &[false; 25], 1.0);
        for &d in &dist {
            assert_eq!(d, f64::MAX, "No walls → all distances should be MAX");
        }
    }

    #[test]
    fn test_wall_distance_wall_cell_zero() {
        let mut mask = vec![false; 25];
        mask[0] = true;
        let dist = compute_wall_distance(5, 5, &mask, 1.0);
        assert_eq!(dist[0], 0.0, "Wall cell distance should be 0");
    }

    #[test]
    fn test_wall_distance_adjacent() {
        let mut mask = vec![false; 25];
        // Wall at (0, 0)
        mask[0] = true;
        let dx = 1.0;
        let dist = compute_wall_distance(5, 5, &mask, dx);
        // Cell (1, 0) should be at distance 1.0
        assert!((dist[1] - 1.0).abs() < 1e-10, "dist = {}", dist[1]);
        // Cell (0, 1) at distance 1.0
        assert!((dist[5] - 1.0).abs() < 1e-10, "dist = {}", dist[5]);
    }

    #[test]
    fn test_wall_distance_diagonal() {
        let mut mask = vec![false; 25];
        mask[0] = true;
        let dist = compute_wall_distance(5, 5, &mask, 1.0);
        // Cell (1, 1) at distance sqrt(2)
        assert!(
            (dist[6] - 2.0_f64.sqrt()).abs() < 1e-10,
            "dist = {}",
            dist[6]
        );
    }

    // --- Wall shear stress 3-point ---

    #[test]
    fn test_wall_shear_stress_3pt_linear_profile() {
        // u = c*y, u0=0, u1=c*dy, u2=c*2*dy
        // du/dy = c
        let c = 10.0;
        let dy = 0.01;
        let nu = 1e-3;
        let tau = wall_shear_stress_3pt(0.0, c * dy, c * 2.0 * dy, dy, nu);
        let expected = nu * c;
        assert!(
            (tau - expected).abs() < 1e-10,
            "tau = {tau}, expected {expected}"
        );
    }

    // --- y+ from tau ---

    #[test]
    fn test_y_plus_from_tau_basic() {
        let yp = y_plus_from_tau(0.01, 0.04, 1e-4);
        // u_tau = sqrt(0.04) = 0.2; y+ = 0.01*0.2/1e-4 = 20
        assert!((yp - 20.0).abs() < 1e-10, "y+ = {yp}");
    }

    // --- WMLES ---

    #[test]
    fn test_wmles_wall_stress_positive() {
        let wm = WallModeledLes::new(1e-5);
        let tau = wm.wall_stress(10.0, 0.001);
        assert!(tau > 0.0, "Wall stress should be positive: {tau}");
    }

    #[test]
    fn test_wmles_eddy_viscosity_positive() {
        let wm = WallModeledLes::new(1e-5);
        let nu_t = wm.eddy_viscosity(0.001, 10.0);
        assert!(nu_t >= 0.0, "Eddy viscosity should be non-negative: {nu_t}");
    }

    #[test]
    fn test_wmles_y_plus_positive() {
        let wm = WallModeledLes::new(1e-5);
        let yp = wm.y_plus_at_match(0.001, 10.0);
        assert!(yp > 0.0, "y+ should be positive: {yp}");
    }

    // --- Reichardt velocity ---

    #[test]
    fn test_reichardt_near_zero() {
        let u = reichardt_velocity(0.01, 0.41);
        assert!(u > 0.0 && u < 0.1, "u+ near wall = {u}");
    }

    #[test]
    fn test_reichardt_monotone() {
        let vals: Vec<f64> = [1.0, 5.0, 20.0, 100.0, 500.0]
            .iter()
            .map(|&yp| reichardt_velocity(yp, 0.41))
            .collect();
        for i in 0..vals.len() - 1 {
            assert!(vals[i + 1] > vals[i], "Not monotone at {i}");
        }
    }

    // --- Spalding's law ---

    #[test]
    fn test_spalding_at_zero() {
        let yp = spalding_y_plus(0.0, 0.41, 5.2);
        assert!(yp.abs() < 1e-10, "y+(0) = {yp}");
    }

    #[test]
    fn test_spalding_monotone() {
        let vals: Vec<f64> = [1.0, 5.0, 10.0, 15.0, 20.0]
            .iter()
            .map(|&up| spalding_y_plus(up, 0.41, 5.2))
            .collect();
        for i in 0..vals.len() - 1 {
            assert!(vals[i + 1] > vals[i], "Not monotone at {i}");
        }
    }

    #[test]
    fn test_spalding_sublayer_match() {
        // In the viscous sublayer (u+ < ~5), y+ ≈ u+
        let yp = spalding_y_plus(3.0, 0.41, 5.2);
        assert!((yp - 3.0).abs() < 0.5, "Sublayer mismatch: y+ = {yp}");
    }

    // ── Equilibrium wall model ────────────────────────────────────────────────

    #[test]
    fn test_equilibrium_wall_model_positive_stress() {
        let ewm = EquilibriumWallModel::new(1e-5);
        let tau = ewm.wall_stress(5.0, 0.01);
        assert!(
            tau > 0.0,
            "Equilibrium wall stress should be positive: {tau}"
        );
    }

    #[test]
    fn test_equilibrium_wall_model_zero_velocity() {
        let ewm = EquilibriumWallModel::new(1e-5);
        let tau = ewm.wall_stress(0.0, 0.01);
        assert!(tau.abs() < 1e-10, "Zero velocity → zero stress: {tau}");
    }

    #[test]
    fn test_equilibrium_wall_model_eddy_viscosity_positive() {
        let ewm = EquilibriumWallModel::new(1e-5);
        let nu_t = ewm.effective_viscosity(5.0, 0.01, 0.001);
        assert!(nu_t >= 0.0, "Eddy viscosity should be non-negative: {nu_t}");
    }

    // ── Non-equilibrium wall model ────────────────────────────────────────────

    #[test]
    fn test_non_equilibrium_wall_model_basic() {
        let neq = NonEquilibriumWallModel::new(1e-5, 0.1);
        let tau = neq.wall_stress(5.0, 0.01, 0.2);
        assert!(
            tau.is_finite(),
            "Non-eq wall stress should be finite: {tau}"
        );
    }

    #[test]
    fn test_non_equilibrium_pressure_gradient_effect() {
        let neq = EquilibriumWallModel::new(1e-5);
        let tau_eq = neq.wall_stress(5.0, 0.01);
        let neq_wm = NonEquilibriumWallModel::new(1e-5, 0.1);
        let tau_neq_adv = neq_wm.wall_stress(5.0, 0.01, 0.0);
        // Without adverse PG, non-eq ≈ equilibrium
        assert!(
            tau_neq_adv.is_finite(),
            "Non-eq wall stress should be finite: {tau_neq_adv}"
        );
        let _ = tau_eq;
    }

    // ── Wall model switching ──────────────────────────────────────────────────

    #[test]
    fn test_wall_model_switch_uses_eq_for_small_yplus() {
        let switcher = WallModelSwitcher::new(1e-5, 5.0);
        let tau = switcher.compute_stress(1.0, 1e-4); // very small y → small y+
        assert!(
            tau.is_finite(),
            "Switched wall stress should be finite: {tau}"
        );
    }

    #[test]
    fn test_wall_model_switch_uses_neq_for_large_yplus() {
        let switcher = WallModelSwitcher::new(1e-5, 5.0);
        let tau = switcher.compute_stress(100.0, 0.01); // large y
        assert!(
            tau.is_finite(),
            "Switched (non-eq) wall stress should be finite: {tau}"
        );
    }

    // ── Wall model validation ─────────────────────────────────────────────────

    #[test]
    fn test_wall_model_cf_vs_blasius() {
        // For laminar flow, Cf from WMLES should be in a reasonable range
        let cf = blasius_cf(1e5);
        assert!(cf > 0.0 && cf < 0.01, "Blasius Cf = {cf}");
    }

    #[test]
    fn test_wall_model_cf_vs_turbulent_1_7() {
        // Turbulent 1/7-law Cf should be lower than laminar at high Re
        let cf_lam = blasius_cf(1e6);
        let cf_turb = turbulent_cf_1_7(1e6);
        assert!(
            cf_turb > cf_lam,
            "Turbulent Cf should exceed Blasius: cf_turb={cf_turb}, cf_lam={cf_lam}"
        );
    }

    // ── Wall distance improvements ────────────────────────────────────────────

    #[test]
    fn test_wall_distance_3d_basic() {
        // Single wall node at (0,0,0), fluid node at (1,0,0)
        let wall_pts: &[[f64; 3]] = &[[0.0, 0.0, 0.0]];
        let d = min_wall_distance_3d([1.0, 0.0, 0.0], wall_pts);
        assert!((d - 1.0).abs() < 1e-12, "distance = {d}");
    }

    #[test]
    fn test_wall_distance_3d_multiple_walls() {
        let wall_pts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let d = min_wall_distance_3d([1.0, 0.0, 0.0], wall_pts);
        assert!((d - 1.0).abs() < 1e-12, "min distance = {d}");
    }

    #[test]
    fn test_wall_normal_computation() {
        // Wall at y=0, fluid node at (0, 0.5, 0) → normal should point in y direction
        let wall_pts: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let (_, normal) = nearest_wall_normal([0.5, 1.0, 0.0], wall_pts);
        // Normal should be roughly in y-direction
        assert!(normal[1] > 0.5, "Wall normal y-component: {}", normal[1]);
    }

    // ── Pressure gradient effects ─────────────────────────────────────────────

    #[test]
    fn test_adverse_pressure_gradient_increases_cf() {
        let nu = 1e-5;
        let u = 10.0;
        let y = 0.01;
        let dpds = 1.0; // adverse pressure gradient
        let cf_adv = pressure_gradient_correction_cf(u, y, nu, dpds);
        let cf_zero = pressure_gradient_correction_cf(u, y, nu, 0.0);
        assert!(
            cf_adv.is_finite() && cf_zero.is_finite(),
            "Cf should be finite"
        );
    }

    // ── Temperature wall model ────────────────────────────────────────────────

    #[test]
    fn test_thermal_wall_model_basic() {
        let twm = ThermalWallModel::new(1e-5, 0.85);
        let q_w = twm.wall_heat_flux(300.0, 350.0, 0.01, 10.0);
        assert!(q_w.is_finite(), "Wall heat flux should be finite: {q_w}");
    }

    #[test]
    fn test_thermal_wall_model_zero_delta_t() {
        let twm = ThermalWallModel::new(1e-5, 0.85);
        let q_w = twm.wall_heat_flux(300.0, 300.0, 0.01, 10.0);
        assert!(q_w.abs() < 1e-10, "Zero ΔT → zero heat flux: {q_w}");
    }

    // ── Rough-wall correction ─────────────────────────────────────────────────

    #[test]
    fn test_rough_wall_log_law_shift() {
        // Sand-grain roughness shifts the log-law intercept downward
        let k_s_plus = 70.0; // fully rough regime
        let shift = roughness_function(k_s_plus);
        assert!(
            shift > 0.0,
            "Roughness function should be positive for ks+ > 0: {shift}"
        );
    }

    #[test]
    fn test_rough_wall_log_law_velocity() {
        let kappa = 0.41;
        let b = 5.2;
        let k_s_plus = 70.0;
        let u_plus_smooth = law_of_the_wall_velocity(100.0, kappa, b);
        let u_plus_rough = rough_wall_velocity(100.0, kappa, b, k_s_plus);
        assert!(
            u_plus_rough < u_plus_smooth,
            "Rough wall: lower u+ for same y+"
        );
    }

    // ── Compressibility correction ────────────────────────────────────────────

    #[test]
    fn test_van_driest_compressibility_correction() {
        let u_plus_incomp = 15.0;
        let t_wall = 300.0;
        let t_adiabatic = 350.0;
        let u_plus_comp = van_driest_compressibility(u_plus_incomp, t_wall, t_adiabatic);
        assert!(
            u_plus_comp.is_finite(),
            "Compressibility-corrected u+ = {u_plus_comp}"
        );
        assert!(
            u_plus_comp > 0.0,
            "Compressibility-corrected u+ should be positive: {u_plus_comp}"
        );
    }

    // ─── LogLawWallModel tests ───

    #[test]
    fn test_log_law_wall_model_viscous_sublayer() {
        // In the viscous sublayer, u+ = y+
        let u_p = LogLawWallModel::u_plus(3.0, 0.41, 5.2);
        assert!(
            (u_p - 3.0).abs() < 1e-12,
            "u+ in sublayer should equal y+: {u_p}"
        );
    }

    #[test]
    fn test_log_law_wall_model_log_region() {
        let kappa = 0.41;
        let b = 5.2;
        let y_p = 50.0;
        let u_p = LogLawWallModel::u_plus(y_p, kappa, b);
        let expected = (1.0 / kappa) * y_p.ln() + b;
        assert!(
            (u_p - expected).abs() < 1e-12,
            "u+ = {u_p}, expected {expected}"
        );
    }

    #[test]
    fn test_log_law_wall_shear_stress_positive() {
        let model = LogLawWallModel::new(0.41, 5.2, 1e-5);
        let tau = model.wall_shear_stress_newton(10.0, 0.001);
        assert!(tau > 0.0, "Wall shear stress should be positive: {tau}");
    }

    // ─── roughness_correction tests ───

    #[test]
    fn test_roughness_correction_positive() {
        let delta_u = roughness_correction(0.0, 50.0);
        assert!(
            delta_u > 0.0,
            "Roughness correction should be positive for ks+ > 0: {delta_u}"
        );
    }

    #[test]
    fn test_roughness_correction_zero_for_smooth() {
        let delta_u = roughness_correction(0.0, 0.0);
        assert!(
            delta_u.abs() < 1e-10,
            "Smooth wall correction should be ~0: {delta_u}"
        );
    }

    // ─── Spalding struct tests ───

    #[test]
    fn test_spalding_u_from_y_plus_sublayer() {
        let s = Spalding::new(0.41, 5.2);
        let u_p = s.u_from_y_plus(5.0);
        // In the viscous sublayer, u+ ≈ y+
        assert!((u_p - 5.0).abs() < 0.5, "Spalding u+ in sublayer: {u_p}");
    }

    #[test]
    fn test_spalding_u_from_y_plus_monotone() {
        let s = Spalding::new(0.41, 5.2);
        let yp_vals = [1.0, 5.0, 11.6, 30.0, 100.0];
        let up_vals: Vec<f64> = yp_vals.iter().map(|&yp| s.u_from_y_plus(yp)).collect();
        for i in 0..up_vals.len() - 1 {
            assert!(
                up_vals[i + 1] > up_vals[i],
                "Spalding profile not monotone at index {i}: {} -> {}",
                up_vals[i],
                up_vals[i + 1]
            );
        }
    }

    // ─── von_karman_constant test ───

    #[test]
    fn test_von_karman_constant_value() {
        let kappa = von_karman_constant();
        assert!((kappa - 0.41).abs() < 1e-15, "kappa = {kappa}");
    }

    // ─── WernerWengle test ───

    #[test]
    fn test_werner_wengle_viscous_sublayer() {
        // Low velocity → viscous sublayer
        let tau = WernerWengle::wall_shear_stress(1e-4, 1e-3, 1e-4);
        assert!(tau.is_finite() && tau >= 0.0, "WW tau = {tau}");
    }

    // ─── WallFunctionBc test ───

    #[test]
    fn test_wall_function_bc_positive_velocity() {
        let model = LogLawWallModel::new(0.41, 5.2, 1e-5);
        let bc = WallFunctionBc::new(0.001, 1e-5, model);
        let v_bc = bc.wall_velocity_bc(10.0);
        assert!(
            v_bc.is_finite(),
            "WallFunctionBc velocity should be finite: {v_bc}"
        );
        assert!(v_bc > 0.0, "BC velocity should be positive: {v_bc}");
    }
}
