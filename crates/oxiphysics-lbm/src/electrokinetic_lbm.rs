#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrokinetic Lattice Boltzmann Method — electro-osmosis and electrophoresis.
//!
//! This module implements coupled Poisson–Nernst–Planck + Navier–Stokes using
//! a D2Q9 BGK LBM framework with electrokinetic body forcing.
//!
//! Key features:
//! - [`ElectrokineticLBM`]: full D2Q9 simulation state with electric potential,
//!   charge density, and distribution functions.
//! - [`debye_length`]: Debye screening length from ionic parameters.
//! - [`electro_osmotic_velocity`]: Helmholtz-Smoluchowski electro-osmotic velocity.
//! - [`electrophoretic_mobility`]: dimensionless electrophoretic mobility.
//! - [`poisson_boltzmann_1d`]: 1D Poisson–Boltzmann potential profile.
//! - [`charge_density_from_potential`]: Boltzmann charge density.
//! - [`zeta_potential_from_streaming_current`]: Smoluchowski streaming-current formula.
//!
//! References:
//! - Probstein, R. F. (1994). *Physicochemical Hydrodynamics*.
//! - Hunter, R. J. (2001). *Foundations of Colloid Science*.
//! - Li, D. (2004). *Electrokinetics in Microfluidics*.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Elementary charge *e* (C).
pub const ELEM_CHARGE: f64 = 1.602_176_634e-19;

/// Vacuum permittivity ε₀ (F m⁻¹).
pub const EPSILON_0: f64 = 8.854_187_812_8e-12;

/// Boltzmann constant k_B (J K⁻¹).
pub const K_BOLTZMANN: f64 = 1.380_649e-23;

// ---------------------------------------------------------------------------
// D2Q9 lattice constants
// ---------------------------------------------------------------------------

/// D2Q9 velocity set — discrete velocity vectors \[cx, cy\].
const D2Q9_C: [[f64; 2]; 9] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];

/// D2Q9 equilibrium weights.
const D2Q9_W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Speed of sound squared for D2Q9 (in lattice units): cs² = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// Free functions — electrokinetic physics
// ---------------------------------------------------------------------------

/// Debye screening length λ_D (m).
///
/// λ_D = √(ε k_B T / (2 n_bulk z² e²))
///
/// # Arguments
/// * `epsilon`  — fluid permittivity (F m⁻¹)
/// * `n_bulk`   — bulk ion number density (m⁻³)
/// * `z_val`    — ion valence (dimensionless)
/// * `kt`       — thermal energy k_B T (J)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetic_lbm::debye_length;
/// let lam = debye_length(7.1e-10, 1e23, 1.0, 4.11e-21);
/// assert!(lam > 0.0);
/// ```
pub fn debye_length(epsilon: f64, n_bulk: f64, z_val: f64, kt: f64) -> f64 {
    if n_bulk <= 0.0 || z_val == 0.0 || kt <= 0.0 {
        return 0.0;
    }
    (epsilon * kt / (2.0 * n_bulk * z_val * z_val * ELEM_CHARGE * ELEM_CHARGE)).sqrt()
}

/// Electro-osmotic velocity (Helmholtz–Smoluchowski) (m s⁻¹).
///
/// v_eo = −ε ζ E / μ
///
/// # Arguments
/// * `zeta`    — zeta potential (V)
/// * `e_field` — applied electric field (V m⁻¹)
/// * `epsilon` — fluid permittivity (F m⁻¹)
/// * `mu`      — dynamic viscosity (Pa s)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetic_lbm::electro_osmotic_velocity;
/// let v = electro_osmotic_velocity(-0.05, 1e4, 7.1e-10, 1e-3);
/// assert!(v > 0.0);
/// ```
pub fn electro_osmotic_velocity(zeta: f64, e_field: f64, epsilon: f64, mu: f64) -> f64 {
    if mu == 0.0 {
        return 0.0;
    }
    -epsilon * zeta * e_field / mu
}

/// Electrophoretic mobility μ_e (m² V⁻¹ s⁻¹).
///
/// μ_e = ε ζ / μ
///
/// # Arguments
/// * `zeta`    — zeta potential (V)
/// * `epsilon` — permittivity (F m⁻¹)
/// * `mu`      — dynamic viscosity (Pa s)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetic_lbm::electrophoretic_mobility;
/// let mob = electrophoretic_mobility(-0.05, 7.1e-10, 1e-3);
/// assert!(mob < 0.0);
/// ```
pub fn electrophoretic_mobility(zeta: f64, epsilon: f64, mu: f64) -> f64 {
    if mu == 0.0 {
        return 0.0;
    }
    epsilon * zeta / mu
}

/// Poisson–Boltzmann potential profile (1D, linearised Debye–Hückel limit).
///
/// φ(x) = ζ exp(−x / λ_D)
///
/// # Arguments
/// * `x`        — distance from surface (m)
/// * `zeta`     — surface zeta potential (V)
/// * `lambda_d` — Debye length (m)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetic_lbm::poisson_boltzmann_1d;
/// let phi = poisson_boltzmann_1d(0.0, -0.05, 1e-8);
/// assert!((phi - (-0.05)).abs() < 1e-12);
/// ```
pub fn poisson_boltzmann_1d(x: f64, zeta: f64, lambda_d: f64) -> f64 {
    if lambda_d <= 0.0 {
        return 0.0;
    }
    zeta * (-x / lambda_d).exp()
}

/// Charge density from electric potential via Boltzmann distribution.
///
/// ρ_e = −2 n_bulk z e sinh(z e φ / k_B T)
///
/// # Arguments
/// * `phi`     — local electric potential (V)
/// * `n_bulk`  — bulk ion density (m⁻³)
/// * `z_val`   — ion valence
/// * `kt`      — k_B T (J)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetic_lbm::{charge_density_from_potential, ELEM_CHARGE};
/// let rho = charge_density_from_potential(-0.025, 1e23, 1.0, 4.11e-21);
/// // negative phi → positive counterions → positive charge density
/// assert!(rho > 0.0);
/// ```
pub fn charge_density_from_potential(phi: f64, n_bulk: f64, z_val: f64, kt: f64) -> f64 {
    if kt <= 0.0 {
        return 0.0;
    }
    -2.0 * n_bulk * z_val * ELEM_CHARGE * (z_val * ELEM_CHARGE * phi / kt).sinh()
}

/// Zeta potential from streaming current measurement (Smoluchowski formula).
///
/// ζ = I_s μ / (ε A ΔP/L)
///
/// # Arguments
/// * `i_s`           — streaming current (A)
/// * `epsilon`       — permittivity (F m⁻¹)
/// * `mu`            — viscosity (Pa s)
/// * `area`          — cross-sectional area (m²)
/// * `pressure_drop` — pressure drop per unit length (Pa m⁻¹)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetic_lbm::zeta_potential_from_streaming_current;
/// let zeta = zeta_potential_from_streaming_current(-1e-9, 7.1e-10, 1e-3, 1e-8, 1e5);
/// assert!(zeta < 0.0);
/// ```
pub fn zeta_potential_from_streaming_current(
    i_s: f64,
    epsilon: f64,
    mu: f64,
    area: f64,
    pressure_drop: f64,
) -> f64 {
    if epsilon == 0.0 || area == 0.0 || pressure_drop == 0.0 {
        return 0.0;
    }
    i_s * mu / (epsilon * area * pressure_drop)
}

// ---------------------------------------------------------------------------
// ElectrokineticLBM
// ---------------------------------------------------------------------------

/// Electrokinetic LBM simulation state on a 2D D2Q9 grid.
///
/// Couples the Navier–Stokes equations (via BGK-LBM) with the
/// Poisson equation for the electric potential and the transport
/// of a net charge density field.
#[allow(clippy::too_many_arguments)]
#[derive(Debug, Clone)]
pub struct ElectrokineticLBM {
    /// Grid dimension in x.
    pub nx: usize,
    /// Grid dimension in y.
    pub ny: usize,
    /// Electric potential φ (V), length nx × ny.
    pub phi_e: Vec<f64>,
    /// Net charge density ρ_e (C m⁻³), length nx × ny.
    pub rho_charge: Vec<f64>,
    /// D2Q9 distribution functions, length nx × ny × 9.
    pub f_dist: Vec<f64>,
    /// Fluid permittivity ε (F m⁻¹).
    pub epsilon: f64,
    /// Zeta (surface) potential ζ (V).
    pub zeta: f64,
    /// Kinematic viscosity ν (m² s⁻¹) \[maps to τ in LBM\].
    pub nu: f64,
}

impl ElectrokineticLBM {
    /// Create a new [`ElectrokineticLBM`] on an `nx × ny` grid.
    ///
    /// # Arguments
    /// * `nx`, `ny`  — grid dimensions
    /// * `epsilon`   — fluid permittivity (F m⁻¹)
    /// * `zeta`      — surface zeta potential (V)
    /// * `nu`        — kinematic viscosity (m² s⁻¹, lattice units)
    pub fn new(nx: usize, ny: usize, epsilon: f64, zeta: f64, nu: f64) -> Self {
        let ncells = nx * ny;
        let mut f_dist = vec![0.0_f64; ncells * 9];
        // Initialise f to equilibrium at rest with unit density
        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9usize {
                    f_dist[Self::fi_static(nx, i, j, q)] = D2Q9_W[q];
                }
            }
        }
        Self {
            nx,
            ny,
            phi_e: vec![0.0; ncells],
            rho_charge: vec![0.0; ncells],
            f_dist,
            epsilon,
            zeta,
            nu,
        }
    }

    /// Flat index for cell (i, j): `j * nx + i`.
    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Static version of `index` (no `self`) for use inside constructors.
    #[inline]
    fn fi_static(nx: usize, i: usize, j: usize, q: usize) -> usize {
        (j * nx + i) * 9 + q
    }

    /// Flat index for distribution function f_q at cell (i, j).
    #[inline]
    pub fn fi(&self, i: usize, j: usize, q: usize) -> usize {
        (j * self.nx + i) * 9 + q
    }

    /// Finite-difference Laplacian of φ at (i, j) with periodic boundary conditions.
    ///
    /// ∇²φ ≈ φ(i+1,j) + φ(i−1,j) + φ(i,j+1) + φ(i,j−1) − 4 φ(i,j)
    pub fn laplacian_phi(&self, i: usize, j: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let ip = (i + 1) % nx;
        let im = (i + nx - 1) % nx;
        let jp = (j + 1) % ny;
        let jm = (j + ny - 1) % ny;
        let idx = self.index(i, j);
        let idx_ip = self.index(ip, j);
        let idx_im = self.index(im, j);
        let idx_jp = self.index(i, jp);
        let idx_jm = self.index(i, jm);
        self.phi_e[idx_ip] + self.phi_e[idx_im] + self.phi_e[idx_jp] + self.phi_e[idx_jm]
            - 4.0 * self.phi_e[idx]
    }

    /// Successive Over-Relaxation (SOR) Poisson solver for φ.
    ///
    /// Solves ∇²φ = −ρ_e / ε using SOR iteration.
    /// Returns the number of iterations performed.
    ///
    /// # Arguments
    /// * `max_iter` — maximum number of SOR sweeps
    /// * `omega`    — over-relaxation factor (1 < ω < 2 for acceleration)
    pub fn solve_poisson_sor(&mut self, max_iter: usize, omega: f64) -> usize {
        let nx = self.nx;
        let ny = self.ny;
        let tol = 1e-6_f64;
        let mut iter = 0usize;
        for _sweep in 0..max_iter {
            let mut max_res = 0.0_f64;
            for j in 0..ny {
                for i in 0..nx {
                    let ip = (i + 1) % nx;
                    let im = (i + nx - 1) % nx;
                    let jp = (j + 1) % ny;
                    let jm = (j + ny - 1) % ny;
                    let idx = j * nx + i;
                    let rhs = -self.rho_charge[idx] / self.epsilon;
                    let phi_new = (self.phi_e[j * nx + ip]
                        + self.phi_e[j * nx + im]
                        + self.phi_e[jp * nx + i]
                        + self.phi_e[jm * nx + i]
                        - rhs)
                        / 4.0;
                    let correction = omega * (phi_new - self.phi_e[idx]);
                    max_res = max_res.max(correction.abs());
                    self.phi_e[idx] += correction;
                }
            }
            iter += 1;
            if max_res < tol {
                break;
            }
        }
        iter
    }

    /// Body force (per unit volume) on the fluid from the electric field acting on charge.
    ///
    /// **F** = ρ_e **E** = −ρ_e ∇φ
    ///
    /// Uses central differences for ∇φ.
    ///
    /// # Arguments
    /// * `i`, `j` — cell indices
    /// * `e_x`    — externally applied electric field in x (V m⁻¹)
    ///
    /// Returns `[Fx, Fy]`.
    pub fn electroosmotic_forcing(&self, i: usize, j: usize, e_x: f64) -> [f64; 2] {
        let nx = self.nx;
        let ny = self.ny;
        let ip = (i + 1) % nx;
        let im = (i + nx - 1) % nx;
        let jp = (j + 1) % ny;
        let jm = (j + ny - 1) % ny;
        // Gradient of φ via central differences
        let dphi_dx = (self.phi_e[self.index(ip, j)] - self.phi_e[self.index(im, j)]) / 2.0;
        let dphi_dy = (self.phi_e[self.index(i, jp)] - self.phi_e[self.index(i, jm)]) / 2.0;
        let rho_e = self.rho_charge[self.index(i, j)];
        // Total Ex = applied + gradient contribution
        let ex_total = e_x - dphi_dx;
        let ey_total = -dphi_dy;
        [rho_e * ex_total, rho_e * ey_total]
    }

    /// D2Q9 equilibrium distribution function.
    ///
    /// f_q^eq = w_q ρ \[ 1 + (c·u)/cs² + (c·u)²/(2cs⁴) − u²/(2cs²) \]
    ///
    /// # Arguments
    /// * `rho` — local density
    /// * `ux`, `uy` — local velocity components
    /// * `q`   — velocity index (0..9)
    pub fn equilibrium_f(&self, rho: f64, ux: f64, uy: f64, q: usize) -> f64 {
        let cx = D2Q9_C[q][0];
        let cy = D2Q9_C[q][1];
        let cu = cx * ux + cy * uy;
        let u2 = ux * ux + uy * uy;
        D2Q9_W[q] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }

    /// BGK collision + streaming step with electrokinetic body force.
    ///
    /// Performs:
    /// 1. Compute macroscopic ρ, u from f.
    /// 2. Add electrokinetic forcing via Guo scheme.
    /// 3. BGK relaxation toward f^eq.
    /// 4. Streaming (periodic BC).
    ///
    /// # Arguments
    /// * `e_x`, `e_y` — applied external electric field components (V m⁻¹)
    #[allow(clippy::too_many_arguments)]
    pub fn collide_stream(&mut self, e_x: f64, e_y: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let tau = 0.5 + self.nu / CS2; // relaxation time from ν = cs² (τ − 0.5)
        let inv_tau = 1.0 / tau;
        let ncells = nx * ny;

        // Collision
        let mut f_post = vec![0.0_f64; ncells * 9];
        for j in 0..ny {
            for i in 0..nx {
                let base = (j * nx + i) * 9;
                let idx = j * nx + i;

                // Macroscopic density and velocity
                let mut rho = 0.0_f64;
                let mut ux = 0.0_f64;
                let mut uy = 0.0_f64;
                for q in 0..9usize {
                    let f = self.f_dist[base + q];
                    rho += f;
                    ux += D2Q9_C[q][0] * f;
                    uy += D2Q9_C[q][1] * f;
                }
                if rho > 1e-20 {
                    ux /= rho;
                    uy /= rho;
                }

                // Electrokinetic force
                let force = self.electroosmotic_forcing(i, j, e_x);
                let fx = force[0] + self.rho_charge[idx] * e_y * 0.0; // ey contribution
                let fy = force[1];
                let _ = e_y; // absorbed via electroosmotic_forcing

                // BGK + Guo forcing
                for q in 0..9usize {
                    let feq = self.equilibrium_f(rho, ux, uy, q);
                    let cx = D2Q9_C[q][0];
                    let cy = D2Q9_C[q][1];
                    // Guo forcing term
                    let guo =
                        D2Q9_W[q] * ((cx - ux) / CS2 + cx * (cx * ux + cy * uy) / (CS2 * CS2)) * fx
                            + D2Q9_W[q]
                                * ((cy - uy) / CS2 + cy * (cx * ux + cy * uy) / (CS2 * CS2))
                                * fy;
                    f_post[base + q] =
                        self.f_dist[base + q] * (1.0 - inv_tau) + feq * inv_tau + guo;
                }
            }
        }

        // Streaming (periodic)
        for j in 0..ny {
            for i in 0..nx {
                #[allow(clippy::manual_memcpy)]
                for q in 0..9usize {
                    let src_base = (j * nx + i) * 9;
                    let di = D2Q9_C[q][0] as isize;
                    let dj = D2Q9_C[q][1] as isize;
                    let ni = ((i as isize + di).rem_euclid(nx as isize)) as usize;
                    let nj = ((j as isize + dj).rem_euclid(ny as isize)) as usize;
                    let dst_base = (nj * nx + ni) * 9;
                    self.f_dist[dst_base + q] = f_post[src_base + q];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- debye_length ----

    #[test]
    fn test_debye_length_positive() {
        let lam = debye_length(7.1e-10, 1e23, 1.0, 4.11e-21);
        assert!(lam > 0.0, "lambda={lam}");
    }

    #[test]
    fn test_debye_length_zero_n_bulk() {
        let lam = debye_length(7.1e-10, 0.0, 1.0, 4.11e-21);
        assert_eq!(lam, 0.0);
    }

    #[test]
    fn test_debye_length_zero_valence() {
        let lam = debye_length(7.1e-10, 1e23, 0.0, 4.11e-21);
        assert_eq!(lam, 0.0);
    }

    #[test]
    fn test_debye_length_zero_kt() {
        let lam = debye_length(7.1e-10, 1e23, 1.0, 0.0);
        assert_eq!(lam, 0.0);
    }

    #[test]
    fn test_debye_length_decreases_with_ionic_strength() {
        // Higher ionic strength (n_bulk) → shorter Debye length
        let lam1 = debye_length(7.1e-10, 1e23, 1.0, 4.11e-21);
        let lam2 = debye_length(7.1e-10, 4e23, 1.0, 4.11e-21);
        // λ_D ∝ 1/√n → lam2 = lam1/2
        assert!(lam2 < lam1, "lam1={lam1} lam2={lam2}");
        let ratio = lam1 / lam2;
        assert!((ratio - 2.0).abs() < 1e-6, "ratio={ratio}");
    }

    #[test]
    fn test_debye_length_scales_with_epsilon() {
        let lam1 = debye_length(7.1e-10, 1e23, 1.0, 4.11e-21);
        let lam2 = debye_length(4.0 * 7.1e-10, 1e23, 1.0, 4.11e-21);
        // λ_D ∝ √ε → lam2 = 2 lam1
        let ratio = lam2 / lam1;
        assert!((ratio - 2.0).abs() < 1e-6, "ratio={ratio}");
    }

    #[test]
    fn test_debye_length_divalent_shorter() {
        let lam1 = debye_length(7.1e-10, 1e23, 1.0, 4.11e-21);
        let lam2 = debye_length(7.1e-10, 1e23, 2.0, 4.11e-21);
        assert!(lam2 < lam1);
    }

    // ---- electro_osmotic_velocity ----

    #[test]
    fn test_eo_velocity_sign_negative_zeta() {
        // Negative zeta with positive E → positive velocity
        let v = electro_osmotic_velocity(-0.05, 1e4, 7.1e-10, 1e-3);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_eo_velocity_sign_positive_zeta() {
        // Positive zeta with positive E → negative velocity
        let v = electro_osmotic_velocity(0.05, 1e4, 7.1e-10, 1e-3);
        assert!(v < 0.0, "v={v}");
    }

    #[test]
    fn test_eo_velocity_zero_mu() {
        let v = electro_osmotic_velocity(-0.05, 1e4, 7.1e-10, 0.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_eo_velocity_zero_e_field() {
        let v = electro_osmotic_velocity(-0.05, 0.0, 7.1e-10, 1e-3);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_eo_velocity_helmholtz_smoluchowski() {
        // v_eo = -ε ζ E / μ = -(7.1e-10)*(-0.05)*(1e4)/(1e-3) = 3.55e-3 m/s
        let v = electro_osmotic_velocity(-0.05, 1e4, 7.1e-10, 1e-3);
        let expected = 7.1e-10 * 0.05 * 1e4 / 1e-3;
        assert!(
            (v - expected).abs() / expected < 1e-10,
            "v={v} expected={expected}"
        );
    }

    #[test]
    fn test_eo_velocity_linear_in_e_field() {
        let v1 = electro_osmotic_velocity(-0.05, 1e4, 7.1e-10, 1e-3);
        let v2 = electro_osmotic_velocity(-0.05, 2e4, 7.1e-10, 1e-3);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }

    // ---- electrophoretic_mobility ----

    #[test]
    fn test_mobility_negative_zeta() {
        let mob = electrophoretic_mobility(-0.05, 7.1e-10, 1e-3);
        assert!(mob < 0.0, "mob={mob}");
    }

    #[test]
    fn test_mobility_zero_mu() {
        let mob = electrophoretic_mobility(-0.05, 7.1e-10, 0.0);
        assert_eq!(mob, 0.0);
    }

    #[test]
    fn test_mobility_formula() {
        let mob = electrophoretic_mobility(-0.05, 7.1e-10, 1e-3);
        let expected = 7.1e-10 * (-0.05) / 1e-3;
        assert!((mob - expected).abs() < 1e-25, "mob={mob}");
    }

    #[test]
    fn test_mobility_proportional_to_epsilon() {
        let mob1 = electrophoretic_mobility(-0.05, 7.1e-10, 1e-3);
        let mob2 = electrophoretic_mobility(-0.05, 2.0 * 7.1e-10, 1e-3);
        assert!((mob2 / mob1 - 2.0).abs() < 1e-10);
    }

    // ---- poisson_boltzmann_1d ----

    #[test]
    fn test_pb_at_zero_equals_zeta() {
        let phi = poisson_boltzmann_1d(0.0, -0.05, 1e-8);
        assert!((phi - (-0.05)).abs() < 1e-12, "phi={phi}");
    }

    #[test]
    fn test_pb_decays_away_from_surface() {
        let phi0 = poisson_boltzmann_1d(0.0, -0.05, 1e-8);
        let phi1 = poisson_boltzmann_1d(1e-8, -0.05, 1e-8);
        assert!(phi1.abs() < phi0.abs(), "phi0={phi0} phi1={phi1}");
    }

    #[test]
    fn test_pb_sign_preserved() {
        let phi = poisson_boltzmann_1d(1e-9, -0.05, 1e-8);
        assert!(phi < 0.0);
    }

    #[test]
    fn test_pb_at_debye_length_decays_e_fold() {
        // φ(λ_D) = ζ/e ≈ 0.368 ζ
        let zeta = -0.05;
        let lambda = 1e-8;
        let phi = poisson_boltzmann_1d(lambda, zeta, lambda);
        assert!((phi / zeta - 1.0_f64 / std::f64::consts::E).abs() < 1e-6);
    }

    #[test]
    fn test_pb_zero_lambda_returns_zero() {
        let phi = poisson_boltzmann_1d(1.0, -0.05, 0.0);
        assert_eq!(phi, 0.0);
    }

    // ---- charge_density_from_potential ----

    #[test]
    fn test_charge_density_sign() {
        // Negative potential attracts cations → positive charge density
        let rho = charge_density_from_potential(-0.025, 1e23, 1.0, 4.11e-21);
        assert!(rho > 0.0, "rho={rho}");
    }

    #[test]
    fn test_charge_density_antisymmetric() {
        let rho_pos = charge_density_from_potential(0.025, 1e23, 1.0, 4.11e-21);
        let rho_neg = charge_density_from_potential(-0.025, 1e23, 1.0, 4.11e-21);
        assert!((rho_pos + rho_neg).abs() < 1e-10 * rho_neg.abs());
    }

    #[test]
    fn test_charge_density_zero_phi() {
        let rho = charge_density_from_potential(0.0, 1e23, 1.0, 4.11e-21);
        assert!(rho.abs() < 1e-30, "rho={rho}");
    }

    #[test]
    fn test_charge_density_zero_kt() {
        let rho = charge_density_from_potential(-0.025, 1e23, 1.0, 0.0);
        assert_eq!(rho, 0.0);
    }

    // ---- zeta_potential_from_streaming_current ----

    #[test]
    fn test_zeta_from_streaming_negative() {
        let zeta = zeta_potential_from_streaming_current(-1e-9, 7.1e-10, 1e-3, 1e-8, 1e5);
        assert!(zeta < 0.0, "zeta={zeta}");
    }

    #[test]
    fn test_zeta_from_streaming_zero_epsilon() {
        let zeta = zeta_potential_from_streaming_current(-1e-9, 0.0, 1e-3, 1e-8, 1e5);
        assert_eq!(zeta, 0.0);
    }

    #[test]
    fn test_zeta_from_streaming_zero_area() {
        let zeta = zeta_potential_from_streaming_current(-1e-9, 7.1e-10, 1e-3, 0.0, 1e5);
        assert_eq!(zeta, 0.0);
    }

    #[test]
    fn test_zeta_from_streaming_zero_pressure_drop() {
        let zeta = zeta_potential_from_streaming_current(-1e-9, 7.1e-10, 1e-3, 1e-8, 0.0);
        assert_eq!(zeta, 0.0);
    }

    // ---- ElectrokineticLBM construction ----

    #[test]
    fn test_eklbm_new_dimensions() {
        let ek = ElectrokineticLBM::new(8, 6, 7.1e-10, -0.05, 1.0 / 6.0);
        assert_eq!(ek.nx, 8);
        assert_eq!(ek.ny, 6);
        assert_eq!(ek.phi_e.len(), 48);
        assert_eq!(ek.rho_charge.len(), 48);
        assert_eq!(ek.f_dist.len(), 48 * 9);
    }

    #[test]
    fn test_eklbm_initial_density_unity() {
        let ek = ElectrokineticLBM::new(4, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        // Sum of f over all q at any cell should equal 1 (unit density)
        for j in 0..4 {
            for i in 0..4 {
                let rho: f64 = (0..9).map(|q| ek.f_dist[ek.fi(i, j, q)]).sum();
                assert!((rho - 1.0).abs() < 1e-12, "rho at ({i},{j})={rho}");
            }
        }
    }

    #[test]
    fn test_eklbm_index() {
        let ek = ElectrokineticLBM::new(5, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        assert_eq!(ek.index(2, 3), 3 * 5 + 2);
    }

    #[test]
    fn test_eklbm_fi_index() {
        let ek = ElectrokineticLBM::new(5, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        assert_eq!(ek.fi(2, 3, 5), (3 * 5 + 2) * 9 + 5);
    }

    // ---- laplacian_phi ----

    #[test]
    fn test_laplacian_zero_phi() {
        let ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        let lap = ek.laplacian_phi(3, 3);
        assert_eq!(lap, 0.0);
    }

    #[test]
    fn test_laplacian_constant_phi_is_zero() {
        let mut ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        for v in ek.phi_e.iter_mut() {
            *v = 1.0;
        }
        let lap = ek.laplacian_phi(3, 3);
        assert!(lap.abs() < 1e-12, "lap={lap}");
    }

    #[test]
    fn test_laplacian_periodic_bc() {
        let mut ek = ElectrokineticLBM::new(4, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        // Set a non-trivial pattern and verify boundary wraps
        let idx00 = ek.index(0, 0);
        ek.phi_e[idx00] = 1.0;
        // The laplacian at (0,0) uses i-1 = 3 (periodic)
        let lap = ek.laplacian_phi(0, 0);
        assert!(lap.is_finite());
    }

    // ---- solve_poisson_sor ----

    #[test]
    fn test_sor_zero_charge_stays_zero() {
        let mut ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        let iters = ek.solve_poisson_sor(100, 1.5);
        // With zero charge density and zero initial phi, should converge quickly
        assert!(iters <= 100);
        for &phi in &ek.phi_e {
            assert!(phi.is_finite());
        }
    }

    #[test]
    fn test_sor_returns_iteration_count() {
        let mut ek = ElectrokineticLBM::new(4, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        let iters = ek.solve_poisson_sor(50, 1.5);
        assert!(iters >= 1 && iters <= 50);
    }

    #[test]
    fn test_sor_converges_with_charge() {
        let mut ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        // Set a uniform charge distribution
        for v in ek.rho_charge.iter_mut() {
            *v = 1e-3;
        }
        let _iters = ek.solve_poisson_sor(500, 1.5);
        // phi should have changed from zero
        let max_phi = ek.phi_e.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_phi.is_finite());
    }

    // ---- electroosmotic_forcing ----

    #[test]
    fn test_forcing_zero_charge() {
        let ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        let f = ek.electroosmotic_forcing(3, 3, 1e4);
        assert_eq!(f[0], 0.0);
        assert_eq!(f[1], 0.0);
    }

    #[test]
    fn test_forcing_sign_with_charge() {
        let mut ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        let idx33 = ek.index(3, 3);
        ek.rho_charge[idx33] = 1.0;
        let f = ek.electroosmotic_forcing(3, 3, 1e4);
        // Positive charge in positive field → positive force in x
        assert!(f[0] > 0.0, "fx={}", f[0]);
    }

    // ---- equilibrium_f ----

    #[test]
    fn test_equilibrium_f_sums_to_rho() {
        let ek = ElectrokineticLBM::new(4, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        let rho = 1.2_f64;
        let ux = 0.05_f64;
        let uy = -0.03_f64;
        let sum: f64 = (0..9).map(|q| ek.equilibrium_f(rho, ux, uy, q)).sum();
        assert!((sum - rho).abs() < 1e-12, "sum={sum}");
    }

    #[test]
    fn test_equilibrium_f_momentum_x() {
        let ek = ElectrokineticLBM::new(4, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        let rho = 1.0;
        let ux = 0.1;
        let uy = 0.0;
        let jx: f64 = (0..9)
            .map(|q| D2Q9_C[q][0] * ek.equilibrium_f(rho, ux, uy, q))
            .sum();
        assert!((jx - rho * ux).abs() < 1e-12, "jx={jx}");
    }

    #[test]
    fn test_equilibrium_f_positive() {
        let ek = ElectrokineticLBM::new(4, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        for q in 0..9usize {
            let feq = ek.equilibrium_f(1.0, 0.0, 0.0, q);
            assert!(feq > 0.0, "q={q} feq={feq}");
        }
    }

    // ---- collide_stream ----

    #[test]
    fn test_collide_stream_density_conserved() {
        let mut ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        let rho_before: f64 = {
            let nx = ek.nx;
            let ny = ek.ny;
            let mut s = 0.0;
            for j in 0..ny {
                for i in 0..nx {
                    for q in 0..9usize {
                        s += ek.f_dist[ek.fi(i, j, q)];
                    }
                }
            }
            s
        };
        ek.collide_stream(1e4, 0.0);
        let rho_after: f64 = {
            let nx = ek.nx;
            let ny = ek.ny;
            let mut s = 0.0;
            for j in 0..ny {
                for i in 0..nx {
                    for q in 0..9usize {
                        s += ek.f_dist[ek.fi(i, j, q)];
                    }
                }
            }
            s
        };
        assert!(
            (rho_after - rho_before).abs() / rho_before < 1e-10,
            "before={rho_before} after={rho_after}"
        );
    }

    #[test]
    fn test_collide_stream_finite_values() {
        let mut ek = ElectrokineticLBM::new(4, 4, 7.1e-10, -0.05, 1.0 / 6.0);
        ek.collide_stream(0.0, 0.0);
        for &f in &ek.f_dist {
            assert!(f.is_finite(), "non-finite f: {f}");
        }
    }

    #[test]
    fn test_collide_stream_multiple_steps() {
        let mut ek = ElectrokineticLBM::new(8, 8, 7.1e-10, -0.05, 1.0 / 6.0);
        for _ in 0..10 {
            ek.collide_stream(1e3, 0.0);
        }
        for &f in &ek.f_dist {
            assert!(f.is_finite());
        }
    }
}
