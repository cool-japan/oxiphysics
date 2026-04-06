#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermal Lattice Boltzmann method for natural convection.
//!
//! Implements a Double Distribution Function (DDF) approach where a second
//! set of distributions `g_i` solves the advection-diffusion equation for
//! temperature, while `f_i` solves the mass/momentum equations.
//!
//! Buoyancy is incorporated via the Boussinesq approximation.
//!
//! References:
//! - He, X., Chen, S., & Doolen, G. D. (1998). *JCP* 146, 282–300.
//! - Shan, X. (1997). *Phys. Rev. E* 55, 2780.

use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};

/// Speed of sound squared in lattice units (cs² = 1/3).
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// ThermalLbmCell
// ---------------------------------------------------------------------------

/// A single node in the thermal LBM grid.
///
/// Stores both the hydrodynamic distributions `f` and the temperature
/// distributions `g` for the D2Q9 velocity set.
#[allow(dead_code)]
pub struct ThermalLbmCell {
    /// Hydrodynamic distribution functions f\[q\] for q = 0..9.
    pub f: [f64; 9],
    /// Temperature distribution functions g\[q\] for q = 0..9.
    pub g: [f64; 9],
    /// Macroscopic density ρ.
    pub rho: f64,
    /// Macroscopic x-velocity.
    pub ux: f64,
    /// Macroscopic y-velocity.
    pub uy: f64,
    /// Macroscopic temperature T.
    pub temperature: f64,
}

impl ThermalLbmCell {
    /// Create a new cell initialised to equilibrium at given density, velocity, and temperature.
    pub fn new(rho: f64, ux: f64, uy: f64, temperature: f64) -> Self {
        let mut cell = Self {
            f: [0.0; 9],
            g: [0.0; 9],
            rho,
            ux,
            uy,
            temperature,
        };
        for i in 0..9 {
            cell.f[i] = feq(i, rho, ux, uy);
            cell.g[i] = geq(i, temperature, ux, uy);
        }
        cell
    }

    /// Compute macroscopic density from distributions.
    pub fn compute_rho(&self) -> f64 {
        self.f.iter().sum()
    }

    /// Compute macroscopic temperature from g distributions.
    pub fn compute_temperature(&self) -> f64 {
        self.g.iter().sum()
    }

    /// Compute macroscopic velocity components (ux, uy).
    pub fn compute_velocity(&self) -> (f64, f64) {
        let rho = self.compute_rho();
        if rho.abs() < f64::EPSILON {
            return (0.0, 0.0);
        }
        let mut ux = 0.0;
        let mut uy = 0.0;
        for i in 0..9 {
            let cx = D2Q9_VELOCITIES[i][0] as f64;
            let cy = D2Q9_VELOCITIES[i][1] as f64;
            ux += cx * self.f[i];
            uy += cy * self.f[i];
        }
        (ux / rho, uy / rho)
    }
}

// ---------------------------------------------------------------------------
// ThermalD2Q9 — main solver
// ---------------------------------------------------------------------------

/// D2Q9 thermal LBM solver using the Double Distribution Function approach.
///
/// Solves coupled mass/momentum (`f`) and advection-diffusion (`g`) equations.
/// Natural convection is driven by Boussinesq buoyancy.
pub struct ThermalD2Q9 {
    /// Hydrodynamic distribution functions, indexed \[node\]\[direction\].
    pub f: Vec<[f64; 9]>,
    /// Temperature distribution functions, indexed \[node\]\[direction\].
    pub g: Vec<[f64; 9]>,
    /// Macroscopic density at each node.
    pub rho: Vec<f64>,
    /// Macroscopic x-velocity at each node.
    pub ux: Vec<f64>,
    /// Macroscopic y-velocity at each node.
    pub uy: Vec<f64>,
    /// Macroscopic temperature at each node.
    pub temperature: Vec<f64>,
    /// Number of lattice nodes in x.
    pub nx: usize,
    /// Number of lattice nodes in y.
    pub ny: usize,
    /// Hydrodynamic relaxation time τ.
    pub tau: f64,
    /// Thermal relaxation time τ_t.
    pub tau_t: f64,
    /// Thermal diffusivity α = (τ_t − 0.5) / 3.
    pub alpha: f64,
    /// Kinematic viscosity ν = (τ − 0.5) / 3.
    pub nu: f64,
    /// Thermal expansion coefficient β.
    pub beta: f64,
    /// Reference temperature T_ref for Boussinesq approximation.
    pub t_ref: f64,
    /// Gravitational acceleration (y-component, typically negative).
    pub gravity: f64,
}

impl ThermalD2Q9 {
    /// Create a new solver with uniform initial conditions.
    ///
    /// # Arguments
    /// * `nx`, `ny` — grid dimensions
    /// * `tau` — hydrodynamic relaxation time (> 0.5)
    /// * `tau_t` — thermal relaxation time (> 0.5)
    /// * `beta` — thermal expansion coefficient
    /// * `t_ref` — reference temperature
    /// * `gravity` — gravitational acceleration (y-direction)
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        tau: f64,
        tau_t: f64,
        beta: f64,
        t_ref: f64,
        gravity: f64,
    ) -> Self {
        let n = nx * ny;
        let nu = (tau - 0.5) / 3.0;
        let alpha = (tau_t - 0.5) / 3.0;
        let f_eq = feq_uniform(1.0, 0.0, 0.0);
        let g_eq = geq_uniform(t_ref, 0.0, 0.0);
        Self {
            f: vec![f_eq; n],
            g: vec![g_eq; n],
            rho: vec![1.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            temperature: vec![t_ref; n],
            nx,
            ny,
            tau,
            tau_t,
            alpha,
            nu,
            beta,
            t_ref,
            gravity,
        }
    }

    /// Return the flat node index for grid position (ix, iy).
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Perform a single BGK collision step for both f and g distributions.
    pub fn collision_step(&mut self) {
        let n = self.nx * self.ny;
        for node in 0..n {
            let rho = self.rho[node];
            let ux = self.ux[node];
            let uy = self.uy[node];
            let temp = self.temperature[node];

            // Boussinesq body force (y-direction)
            let fy = boussinesq_force(rho, self.beta, self.gravity, temp, self.t_ref);

            for i in 0..9 {
                let cy = D2Q9_VELOCITIES[i][1] as f64;
                let fi_eq = feq(i, rho, ux, uy);
                // Force correction term
                let force_term = D2Q9_WEIGHTS[i] * cy * fy / CS2;
                self.f[node][i] += -(self.f[node][i] - fi_eq) / self.tau + force_term;

                let gi_eq = geq(i, temp, ux, uy);
                self.g[node][i] += -(self.g[node][i] - gi_eq) / self.tau_t;
            }
        }
    }

    /// Perform BGK collision for temperature distribution only.
    pub fn collision_thermal(&mut self) {
        let n = self.nx * self.ny;
        for node in 0..n {
            let temp = self.temperature[node];
            let ux = self.ux[node];
            let uy = self.uy[node];
            for i in 0..9 {
                let gi_eq = geq(i, temp, ux, uy);
                self.g[node][i] += -(self.g[node][i] - gi_eq) / self.tau_t;
            }
        }
    }

    /// Perform streaming (propagation) step with periodic boundaries.
    pub fn streaming_step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;

        let f_old = self.f.clone();
        let g_old = self.g.clone();

        for iy in 0..ny {
            for ix in 0..nx {
                let node = iy * nx + ix;
                for i in 0..9 {
                    let cx = D2Q9_VELOCITIES[i][0];
                    let cy = D2Q9_VELOCITIES[i][1];
                    let sx = ((ix as i64 - cx as i64).rem_euclid(nx as i64)) as usize;
                    let sy = ((iy as i64 - cy as i64).rem_euclid(ny as i64)) as usize;
                    let src = sy * nx + sx;
                    self.f[node][i] = f_old[src][i];
                    self.g[node][i] = g_old[src][i];
                }
            }
        }
        // suppress unused warning
        let _ = n;
    }

    /// Update macroscopic variables (ρ, u, T) from distributions.
    pub fn update_macroscopic(&mut self) {
        let n = self.nx * self.ny;
        for node in 0..n {
            let mut rho = 0.0;
            let mut ux = 0.0;
            let mut uy = 0.0;
            let mut temp = 0.0;
            for i in 0..9 {
                rho += self.f[node][i];
                let cx = D2Q9_VELOCITIES[i][0] as f64;
                let cy = D2Q9_VELOCITIES[i][1] as f64;
                ux += cx * self.f[node][i];
                uy += cy * self.f[node][i];
                temp += self.g[node][i];
            }
            if rho > f64::EPSILON {
                ux /= rho;
                uy /= rho;
            }
            self.rho[node] = rho;
            self.ux[node] = ux;
            self.uy[node] = uy;
            self.temperature[node] = temp;
        }
    }

    /// Compute Nusselt number from wall heat flux on the left wall (ix = 0).
    ///
    /// Nu = q_wall * L / (k * ΔT)
    /// where q_wall is estimated from the temperature gradient at the wall.
    pub fn nusselt_number(&self, t_hot: f64, t_cold: f64) -> f64 {
        nusselt_number(&self.temperature, self.nx, self.ny, t_hot, t_cold)
    }

    /// Compute the Rayleigh number Ra = g β ΔT L³ / (ν α).
    pub fn rayleigh_number(&self, delta_t: f64, length: f64) -> f64 {
        rayleigh_number(
            self.gravity.abs(),
            self.beta,
            delta_t,
            length,
            self.nu,
            self.alpha,
        )
    }

    /// Compute the Prandtl number Pr = ν / α.
    pub fn prandtl_number(&self) -> f64 {
        prandtl_number(self.nu, self.alpha)
    }
}

// ---------------------------------------------------------------------------
// Standalone physics functions
// ---------------------------------------------------------------------------

/// Compute the D2Q9 hydrodynamic equilibrium distribution f_i^eq.
///
/// f_eq_i = w_i ρ \[1 + (e_i·u)/cs² + (e_i·u)²/(2cs⁴) − u²/(2cs²)\]
pub fn feq(i: usize, rho: f64, ux: f64, uy: f64) -> f64 {
    let cx = D2Q9_VELOCITIES[i][0] as f64;
    let cy = D2Q9_VELOCITIES[i][1] as f64;
    let eu = cx * ux + cy * uy;
    let u2 = ux * ux + uy * uy;
    D2Q9_WEIGHTS[i] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
}

/// Compute uniform f_eq array at given (rho, ux, uy).
fn feq_uniform(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let mut out = [0.0; 9];
    for i in 0..9 {
        out[i] = feq(i, rho, ux, uy);
    }
    out
}

/// Compute the D2Q9 thermal equilibrium distribution g_i^eq.
///
/// g_eq_i = w_i T \[1 + (e_i·u)/cs²\]
pub fn geq(i: usize, temperature: f64, ux: f64, uy: f64) -> f64 {
    let cx = D2Q9_VELOCITIES[i][0] as f64;
    let cy = D2Q9_VELOCITIES[i][1] as f64;
    let eu = cx * ux + cy * uy;
    D2Q9_WEIGHTS[i] * temperature * (1.0 + eu / CS2)
}

/// Compute uniform g_eq array at given (temperature, ux, uy).
fn geq_uniform(temperature: f64, ux: f64, uy: f64) -> [f64; 9] {
    let mut out = [0.0; 9];
    for i in 0..9 {
        out[i] = geq(i, temperature, ux, uy);
    }
    out
}

/// Boussinesq buoyancy force in the y-direction.
///
/// F_y = ρ β g (T − T_ref)
///
/// # Arguments
/// * `rho` — local density
/// * `beta` — thermal expansion coefficient
/// * `gravity` — gravitational acceleration (positive upward convention)
/// * `temperature` — local temperature T
/// * `t_ref` — reference temperature T_ref
pub fn boussinesq_force(rho: f64, beta: f64, gravity: f64, temperature: f64, t_ref: f64) -> f64 {
    rho * beta * gravity * (temperature - t_ref)
}

/// Nusselt number from temperature field wall gradient.
///
/// Estimates the heat transfer coefficient using a first-order finite difference
/// of the temperature in the x-direction at the hot wall (ix = 0).
///
/// Nu = (−dT/dx|_wall × L) / ΔT
pub fn nusselt_number(temperature: &[f64], nx: usize, ny: usize, t_hot: f64, t_cold: f64) -> f64 {
    let delta_t = t_hot - t_cold;
    if delta_t.abs() < f64::EPSILON {
        return 0.0;
    }
    let mut flux_sum = 0.0;
    for iy in 0..ny {
        let t0 = temperature[iy * nx]; // wall temperature
        let t1 = temperature[iy * nx + 1]; // one node inside
        let grad = t1 - t0; // dT/dx (negative if hot wall on left)
        flux_sum += grad.abs();
    }
    let avg_grad = flux_sum / ny as f64;
    // Normalise by delta_t; L = nx - 1 (lattice units)
    avg_grad * (nx as f64 - 1.0) / delta_t
}

/// Rayleigh number: Ra = g β ΔT L³ / (ν α).
///
/// # Arguments
/// * `g` — gravitational acceleration magnitude
/// * `beta` — thermal expansion coefficient
/// * `delta_t` — temperature difference ΔT
/// * `length` — characteristic length L
/// * `nu` — kinematic viscosity
/// * `alpha` — thermal diffusivity
pub fn rayleigh_number(g: f64, beta: f64, delta_t: f64, length: f64, nu: f64, alpha: f64) -> f64 {
    g * beta * delta_t * length.powi(3) / (nu * alpha)
}

/// Prandtl number: Pr = ν / α.
///
/// # Arguments
/// * `nu` — kinematic viscosity
/// * `alpha` — thermal diffusivity
pub fn prandtl_number(nu: f64, alpha: f64) -> f64 {
    nu / alpha
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // ── equilibrium weights sum ───────────────────────────────────────────

    #[test]
    fn test_feq_sum_equals_rho() {
        let rho = 1.2;
        let sum: f64 = (0..9).map(|i| feq(i, rho, 0.1, -0.05)).sum();
        assert!((sum - rho).abs() < 1e-12, "Σf_eq should equal rho");
    }

    #[test]
    fn test_geq_sum_equals_temperature() {
        let temp = 2.5;
        let sum: f64 = (0..9).map(|i| geq(i, temp, 0.0, 0.0)).sum();
        assert!((sum - temp).abs() < 1e-12, "Σg_eq should equal T");
    }

    #[test]
    fn test_feq_zero_velocity() {
        // At rest, f_eq_0 = 4/9 * rho
        let rho = 1.0;
        let f0 = feq(0, rho, 0.0, 0.0);
        assert!((f0 - 4.0 / 9.0).abs() < EPS);
    }

    #[test]
    fn test_geq_zero_velocity() {
        // At rest, g_eq_0 = 4/9 * T
        let temp = 3.0;
        let g0 = geq(0, temp, 0.0, 0.0);
        assert!((g0 - 4.0 / 9.0 * temp).abs() < EPS);
    }

    #[test]
    fn test_feq_momentum() {
        // Σ e_ix * f_eq_i = rho * ux
        let rho = 1.0;
        let ux = 0.1;
        let uy = 0.0;
        let mom_x: f64 = (0..9)
            .map(|i| D2Q9_VELOCITIES[i][0] as f64 * feq(i, rho, ux, uy))
            .sum();
        assert!((mom_x - rho * ux).abs() < 1e-12);
    }

    #[test]
    fn test_feq_momentum_y() {
        let rho = 1.0;
        let ux = 0.0;
        let uy = 0.05;
        let mom_y: f64 = (0..9)
            .map(|i| D2Q9_VELOCITIES[i][1] as f64 * feq(i, rho, ux, uy))
            .sum();
        assert!((mom_y - rho * uy).abs() < 1e-12);
    }

    // ── Boussinesq force ──────────────────────────────────────────────────

    #[test]
    fn test_boussinesq_force_zero_at_ref() {
        let f = boussinesq_force(1.0, 0.01, -9.8, 300.0, 300.0);
        assert!(f.abs() < EPS);
    }

    #[test]
    fn test_boussinesq_force_positive_above_ref() {
        // Hot fluid rises: T > T_ref, gravity < 0 → force should be in direction of gravity sign
        let f = boussinesq_force(1.0, 0.01, 9.8, 310.0, 300.0);
        assert!(
            f > 0.0,
            "buoyancy force should be positive for T > T_ref with g > 0"
        );
    }

    #[test]
    fn test_boussinesq_force_negative_below_ref() {
        let f = boussinesq_force(1.0, 0.01, 9.8, 290.0, 300.0);
        assert!(f < 0.0);
    }

    #[test]
    fn test_boussinesq_force_scales_linearly_with_beta() {
        let f1 = boussinesq_force(1.0, 0.01, 9.8, 310.0, 300.0);
        let f2 = boussinesq_force(1.0, 0.02, 9.8, 310.0, 300.0);
        assert!((f2 / f1 - 2.0).abs() < EPS);
    }

    #[test]
    fn test_boussinesq_force_scales_with_rho() {
        let f1 = boussinesq_force(1.0, 0.01, 9.8, 310.0, 300.0);
        let f2 = boussinesq_force(2.0, 0.01, 9.8, 310.0, 300.0);
        assert!((f2 / f1 - 2.0).abs() < EPS);
    }

    // ── Rayleigh number ───────────────────────────────────────────────────

    #[test]
    fn test_rayleigh_number_positive() {
        let ra = rayleigh_number(9.8, 3.4e-3, 10.0, 0.1, 1e-6, 1.43e-7);
        assert!(ra > 0.0);
    }

    #[test]
    fn test_rayleigh_number_zero_delta_t() {
        let ra = rayleigh_number(9.8, 3.4e-3, 0.0, 0.1, 1e-6, 1.43e-7);
        assert!(ra.abs() < EPS);
    }

    #[test]
    fn test_rayleigh_number_cubic_length() {
        let ra1 = rayleigh_number(9.8, 1e-3, 1.0, 1.0, 1e-6, 1e-6);
        let ra2 = rayleigh_number(9.8, 1e-3, 1.0, 2.0, 1e-6, 1e-6);
        assert!((ra2 / ra1 - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_rayleigh_lattice_units() {
        // In lattice units: g=0.001, beta=0.01, dT=1.0, L=32, nu=0.1667, alpha=0.1667
        let nu = 1.0 / 6.0;
        let alpha = 1.0 / 6.0;
        let ra = rayleigh_number(0.001, 0.01, 1.0, 32.0, nu, alpha);
        assert!(ra > 0.0);
    }

    // ── Prandtl number ────────────────────────────────────────────────────

    #[test]
    fn test_prandtl_number_air() {
        // Air: Pr ≈ 0.7
        let pr = prandtl_number(1.5e-5, 2.1e-5);
        assert!((pr - 1.5e-5 / 2.1e-5).abs() < EPS);
    }

    #[test]
    fn test_prandtl_number_unity() {
        let pr = prandtl_number(0.1, 0.1);
        assert!((pr - 1.0).abs() < EPS);
    }

    #[test]
    fn test_prandtl_number_water() {
        // Water Pr ≈ 7
        let pr = prandtl_number(7e-7, 1e-7);
        assert!((pr - 7.0).abs() < 1e-12);
    }

    // ── ThermalLbmCell ───────────────────────────────────────────────────

    #[test]
    fn test_cell_new_rho() {
        let cell = ThermalLbmCell::new(1.0, 0.0, 0.0, 300.0);
        assert!((cell.compute_rho() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_cell_new_temperature() {
        let cell = ThermalLbmCell::new(1.0, 0.0, 0.0, 300.0);
        assert!((cell.compute_temperature() - 300.0).abs() < 1e-12);
    }

    #[test]
    fn test_cell_velocity_at_rest() {
        let cell = ThermalLbmCell::new(1.0, 0.0, 0.0, 300.0);
        let (ux, uy) = cell.compute_velocity();
        assert!(ux.abs() < 1e-12);
        assert!(uy.abs() < 1e-12);
    }

    #[test]
    fn test_cell_velocity_nonzero() {
        let cell = ThermalLbmCell::new(1.0, 0.05, -0.03, 300.0);
        let (ux, uy) = cell.compute_velocity();
        assert!((ux - 0.05).abs() < 1e-12);
        assert!((uy - (-0.03)).abs() < 1e-12);
    }

    // ── ThermalD2Q9 solver ───────────────────────────────────────────────

    #[test]
    fn test_solver_new() {
        let solver = ThermalD2Q9::new(8, 8, 0.6, 0.6, 0.01, 300.0, -9.8);
        assert_eq!(solver.nx, 8);
        assert_eq!(solver.ny, 8);
        assert_eq!(solver.f.len(), 64);
        assert_eq!(solver.g.len(), 64);
    }

    #[test]
    fn test_solver_prandtl() {
        let solver = ThermalD2Q9::new(8, 8, 0.8, 0.6, 0.01, 300.0, -9.8);
        let pr = solver.prandtl_number();
        let nu = (0.8 - 0.5) / 3.0;
        let alpha = (0.6 - 0.5) / 3.0;
        assert!((pr - nu / alpha).abs() < 1e-12);
    }

    #[test]
    fn test_solver_rayleigh() {
        let solver = ThermalD2Q9::new(32, 32, 0.6, 0.6, 0.01, 300.0, 9.8);
        let ra = solver.rayleigh_number(10.0, 32.0);
        assert!(ra > 0.0);
    }

    #[test]
    fn test_solver_update_macroscopic_conserves_mass() {
        let mut solver = ThermalD2Q9::new(4, 4, 0.6, 0.6, 0.01, 300.0, -9.8);
        solver.update_macroscopic();
        let total_rho: f64 = solver.rho.iter().sum();
        // 16 nodes × rho=1.0
        assert!((total_rho - 16.0).abs() < 1e-10);
    }

    #[test]
    fn test_solver_collision_thermal_conserves_temperature_sum() {
        let mut solver = ThermalD2Q9::new(4, 4, 0.6, 0.8, 0.01, 300.0, -9.8);
        solver.update_macroscopic();
        let t_before: f64 = solver.temperature.iter().sum();
        solver.collision_thermal();
        solver.update_macroscopic();
        let t_after: f64 = solver.temperature.iter().sum();
        // BGK should conserve total temperature (Σg_i is conserved)
        assert!((t_after - t_before).abs() < 1e-9);
    }

    #[test]
    fn test_solver_idx() {
        let solver = ThermalD2Q9::new(10, 8, 0.6, 0.6, 0.01, 0.0, -9.8);
        assert_eq!(solver.idx(0, 0), 0);
        assert_eq!(solver.idx(5, 3), 35);
        assert_eq!(solver.idx(9, 7), 79);
    }

    #[test]
    fn test_nusselt_number_isothermal_zero() {
        // Uniform temperature → zero gradient → zero Nusselt
        let temp = vec![1.0_f64; 16];
        let nu = nusselt_number(&temp, 4, 4, 2.0, 0.0);
        assert!(nu.abs() < EPS);
    }

    #[test]
    fn test_nusselt_number_positive() {
        // Linear temperature profile: T(x) = t_hot - x * (t_hot - t_cold) / (nx - 1)
        let nx = 4;
        let ny = 4;
        let t_hot = 1.0;
        let t_cold = 0.0;
        let mut temp = vec![0.0_f64; nx * ny];
        for iy in 0..ny {
            for ix in 0..nx {
                temp[iy * nx + ix] = t_hot - ix as f64 * (t_hot - t_cold) / (nx as f64 - 1.0);
            }
        }
        let nu_val = nusselt_number(&temp, nx, ny, t_hot, t_cold);
        assert!(nu_val > 0.0);
    }

    #[test]
    fn test_geq_sum_velocity() {
        // Σ e_ix * g_eq_i = T * ux
        let temp = 2.0;
        let ux = 0.1;
        let uy = 0.0;
        let flux_x: f64 = (0..9)
            .map(|i| D2Q9_VELOCITIES[i][0] as f64 * geq(i, temp, ux, uy))
            .sum();
        assert!((flux_x - temp * ux).abs() < 1e-12);
    }

    #[test]
    fn test_streaming_periodic() {
        // After streaming, total f-mass is conserved
        let mut solver = ThermalD2Q9::new(4, 4, 0.6, 0.6, 0.01, 1.0, -9.8);
        let mass_before: f64 = solver.f.iter().flat_map(|row| row.iter()).sum();
        solver.streaming_step();
        let mass_after: f64 = solver.f.iter().flat_map(|row| row.iter()).sum();
        assert!((mass_after - mass_before).abs() < 1e-10);
    }

    #[test]
    fn test_nu_alpha_from_tau() {
        let tau = 0.8;
        let tau_t = 0.7;
        let solver = ThermalD2Q9::new(4, 4, tau, tau_t, 0.0, 0.0, 0.0);
        let expected_nu = (tau - 0.5) / 3.0;
        let expected_alpha = (tau_t - 0.5) / 3.0;
        assert!((solver.nu - expected_nu).abs() < EPS);
        assert!((solver.alpha - expected_alpha).abs() < EPS);
    }
}
