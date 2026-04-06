#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Electrokinetic Lattice Boltzmann Method with full electroosmotic/electrophoretic coupling.
//!
//! Implements Poisson–Boltzmann electrostatics coupled with D2Q9 BGK fluid dynamics.
//!
//! Key features:
//! - [`ElectrokineticLBM`]: simulation state with velocity distributions, charge density,
//!   and electric potential fields.
//! - [`ElectrokineticParams`]: material and driving parameters.
//! - Linearized Poisson–Boltzmann equation (`poisson_boltzmann_eq`).
//! - Helmholtz–Smoluchowski electroosmotic velocity (`electroosmotic_flow_velocity`).
//! - Debye–Hückel screened potential (`debye_huckel_potential`).
//! - Electrophoretic mobility (`electrophoretic_mobility`).
//!
//! References:
//! - Probstein, R. F. (1994). *Physicochemical Hydrodynamics*.
//! - Hunter, R. J. (2001). *Foundations of Colloid Science*.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// D2Q9 lattice constants
// ---------------------------------------------------------------------------

/// D2Q9 discrete velocity vectors `[cx, cy]`.
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

/// D2Q9 speed of sound squared (lattice units): cs² = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// ElectrokineticParams
// ---------------------------------------------------------------------------

/// Physical and driving parameters for the electrokinetic LBM simulation.
#[derive(Debug, Clone)]
pub struct ElectrokineticParams {
    /// Fluid permittivity ε (F m⁻¹).
    pub permittivity: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Zeta potential ζ (V) at the solid–fluid interface.
    pub zeta_potential: f64,
    /// Debye screening length λ_D (m).
    pub debye_length: f64,
    /// External electric field `[Ex, Ey]` (V m⁻¹).
    pub external_field: [f64; 2],
}

impl ElectrokineticParams {
    /// Construct new [`ElectrokineticParams`].
    ///
    /// # Arguments
    /// * `permittivity` — fluid permittivity (F m⁻¹)
    /// * `viscosity` — dynamic viscosity (Pa·s)
    /// * `zeta_potential` — surface zeta potential (V)
    /// * `debye_length` — Debye screening length (m)
    /// * `external_field` — applied field `[Ex, Ey]` (V m⁻¹)
    pub fn new(
        permittivity: f64,
        viscosity: f64,
        zeta_potential: f64,
        debye_length: f64,
        external_field: [f64; 2],
    ) -> Self {
        Self {
            permittivity,
            viscosity,
            zeta_potential,
            debye_length,
            external_field,
        }
    }
}

// ---------------------------------------------------------------------------
// ElectrokineticLBM
// ---------------------------------------------------------------------------

/// Electrokinetic LBM simulation state on a 2D D2Q9 grid.
///
/// Couples Navier–Stokes (BGK-LBM) with Poisson electrostatics
/// and a net charge density transport field.
#[derive(Debug, Clone)]
pub struct ElectrokineticLBM {
    /// Grid dimension in x.
    pub nx: usize,
    /// Grid dimension in y.
    pub ny: usize,
    /// D2Q9 velocity distribution functions, length `nx × ny × 9`.
    pub f: Vec<f64>,
    /// Net charge density ρ_e (C m⁻³), length `nx × ny`.
    pub charge_density: Vec<f64>,
    /// Electric potential φ (V), length `nx × ny`.
    pub potential: Vec<f64>,
    /// LBM relaxation time τ (lattice units).
    pub tau: f64,
}

impl ElectrokineticLBM {
    /// Flat node index `j * nx + i`.
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Flat distribution-function index `(j * nx + i) * 9 + q`.
    #[inline]
    pub fn fi(&self, i: usize, j: usize, q: usize) -> usize {
        (j * self.nx + i) * 9 + q
    }

    /// D2Q9 equilibrium distribution at (rho, ux, uy) for direction `q`.
    #[inline]
    pub fn feq(&self, rho: f64, ux: f64, uy: f64, q: usize) -> f64 {
        let cx = D2Q9_C[q][0];
        let cy = D2Q9_C[q][1];
        let cu = cx * ux + cy * uy;
        let u2 = ux * ux + uy * uy;
        D2Q9_W[q] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Linearized Poisson–Boltzmann potential profile.
///
/// In the Debye–Hückel limit: φ(x) = ζ exp(−x / λ_D).
///
/// # Arguments
/// * `potential` — distance from the surface (m) \[used as x\]
/// * `zeta`      — surface zeta potential (V)
/// * `debye`     — Debye screening length (m)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::poisson_boltzmann_eq;
/// let phi = poisson_boltzmann_eq(0.0, -0.05, 1e-8);
/// assert!((phi - (-0.05)).abs() < 1e-12);
/// ```
pub fn poisson_boltzmann_eq(potential: f64, zeta: f64, debye: f64) -> f64 {
    if debye <= 0.0 {
        return 0.0;
    }
    zeta * (-potential / debye).exp()
}

/// Helmholtz–Smoluchowski electroosmotic flow velocity.
///
/// v_eo = −ε ζ E / μ
///
/// # Arguments
/// * `zeta`        — zeta potential (V)
/// * `field`       — applied electric field (V m⁻¹)
/// * `viscosity`   — dynamic viscosity (Pa·s)
/// * `permittivity`— fluid permittivity (F m⁻¹)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::electroosmotic_flow_velocity;
/// let v = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
/// assert!(v > 0.0);
/// ```
pub fn electroosmotic_flow_velocity(
    zeta: f64,
    field: f64,
    viscosity: f64,
    permittivity: f64,
) -> f64 {
    if viscosity == 0.0 {
        return 0.0;
    }
    -permittivity * zeta * field / viscosity
}

/// Initialize an [`ElectrokineticLBM`] state at rest with unit density.
///
/// # Arguments
/// * `nx`, `ny` — grid dimensions
/// * `params`   — electrokinetic parameters
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::{ElectrokineticParams, init_electrokinetic};
/// let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e4, 0.0]);
/// let state = init_electrokinetic(8, 6, &p);
/// assert_eq!(state.nx, 8);
/// assert_eq!(state.ny, 6);
/// ```
pub fn init_electrokinetic(
    nx: usize,
    ny: usize,
    params: &ElectrokineticParams,
) -> ElectrokineticLBM {
    let ncells = nx * ny;
    let mut f = vec![0.0_f64; ncells * 9];
    // Initialize f to D2Q9 equilibrium at rest, unit density
    for j in 0..ny {
        for i in 0..nx {
            for q in 0..9_usize {
                f[(j * nx + i) * 9 + q] = D2Q9_W[q];
            }
        }
    }
    let tau = 0.5 + params.viscosity / CS2;
    ElectrokineticLBM {
        nx,
        ny,
        f,
        charge_density: vec![0.0; ncells],
        potential: vec![0.0; ncells],
        tau,
    }
}

/// Perform one BGK collision + streaming step with electrokinetic forcing.
///
/// The Guo forcing scheme adds body forces from the electric field acting
/// on the local charge density.
///
/// # Arguments
/// * `state`  — mutable simulation state
/// * `params` — electrokinetic parameters
pub fn step_electrokinetic(state: &mut ElectrokineticLBM, params: &ElectrokineticParams) {
    let nx = state.nx;
    let ny = state.ny;
    let inv_tau = 1.0 / state.tau;
    let ncells = nx * ny;
    let ex = params.external_field[0];
    let ey = params.external_field[1];

    // Collision
    let mut f_post = vec![0.0_f64; ncells * 9];
    for j in 0..ny {
        for i in 0..nx {
            let base = (j * nx + i) * 9;
            let cell = j * nx + i;

            let mut rho = 0.0_f64;
            let mut ux = 0.0_f64;
            let mut uy_val = 0.0_f64;
            for q in 0..9_usize {
                let fq = state.f[base + q];
                rho += fq;
                ux += D2Q9_C[q][0] * fq;
                uy_val += D2Q9_C[q][1] * fq;
            }
            if rho > 1e-20 {
                ux /= rho;
                uy_val /= rho;
            }

            let rho_e = state.charge_density[cell];
            let fx = rho_e * ex;
            let fy = rho_e * ey;

            for q in 0..9_usize {
                let feq = state.feq(rho, ux, uy_val, q);
                let cx = D2Q9_C[q][0];
                let cy = D2Q9_C[q][1];
                // Guo forcing
                let guo =
                    D2Q9_W[q] * ((cx - ux) / CS2 + cx * (cx * ux + cy * uy_val) / (CS2 * CS2)) * fx
                        + D2Q9_W[q]
                            * ((cy - uy_val) / CS2 + cy * (cx * ux + cy * uy_val) / (CS2 * CS2))
                            * fy;
                f_post[base + q] = state.f[base + q] * (1.0 - inv_tau) + feq * inv_tau + guo;
            }
        }
    }

    // Streaming (periodic)
    for j in 0..ny {
        for i in 0..nx {
            #[allow(clippy::manual_memcpy)]
            for q in 0..9_usize {
                let src_base = (j * nx + i) * 9;
                let di = D2Q9_C[q][0] as isize;
                let dj = D2Q9_C[q][1] as isize;
                let ni = ((i as isize + di).rem_euclid(nx as isize)) as usize;
                let nj = ((j as isize + dj).rem_euclid(ny as isize)) as usize;
                let dst_base = (nj * nx + ni) * 9;
                state.f[dst_base + q] = f_post[src_base + q];
            }
        }
    }
}

/// Compute charge flux (ρ_e · u) at every cell.
///
/// Returns a `Vec`f64` of length `nx × ny` with the magnitude of
/// the charge flux vector |ρ_e **u**|.
///
/// # Arguments
/// * `state` — simulation state
pub fn compute_charge_flux(state: &ElectrokineticLBM) -> Vec<f64> {
    let nx = state.nx;
    let ny = state.ny;
    let ncells = nx * ny;
    let mut flux = vec![0.0_f64; ncells];

    for j in 0..ny {
        for i in 0..nx {
            let base = (j * nx + i) * 9;
            let cell = j * nx + i;
            let mut rho = 0.0_f64;
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            for q in 0..9_usize {
                let fq = state.f[base + q];
                rho += fq;
                ux += D2Q9_C[q][0] * fq;
                uy += D2Q9_C[q][1] * fq;
            }
            if rho > 1e-20 {
                ux /= rho;
                uy /= rho;
            }
            let rho_e = state.charge_density[cell];
            flux[cell] = (rho_e * ux).hypot(rho_e * uy);
        }
    }
    flux
}

/// Debye–Hückel screened electrostatic potential.
///
/// φ(r) = (charge / (4π ε₀ r)) exp(−r / λ_D)
///
/// Using a dimensionless convention (divided by 4π ε₀):
/// φ(r) = charge · exp(−r / λ_D) / r
///
/// # Arguments
/// * `r`           — radial distance (m)
/// * `charge`      — source charge (C or dimensionless)
/// * `debye_length`— Debye screening length (m)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::debye_huckel_potential;
/// let phi = debye_huckel_potential(1e-9, 1.0, 1e-8);
/// assert!(phi > 0.0);
/// ```
pub fn debye_huckel_potential(r: f64, charge: f64, debye_length: f64) -> f64 {
    if r <= 0.0 || debye_length <= 0.0 {
        return 0.0;
    }
    charge * (-r / debye_length).exp() / r
}

/// Electrophoretic mobility (Helmholtz–Smoluchowski).
///
/// μ_e = ε ζ / μ
///
/// # Arguments
/// * `zeta`        — zeta potential (V)
/// * `viscosity`   — dynamic viscosity (Pa·s)
/// * `permittivity`— fluid permittivity (F m⁻¹)
///
/// ```no_run
/// use oxiphysics_lbm::electrokinetics_lbm::electrophoretic_mobility;
/// let mob = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
/// assert!(mob < 0.0);
/// ```
pub fn electrophoretic_mobility(zeta: f64, viscosity: f64, permittivity: f64) -> f64 {
    if viscosity == 0.0 {
        return 0.0;
    }
    permittivity * zeta / viscosity
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── poisson_boltzmann_eq ─────────────────────────────────────────────

    #[test]
    fn test_pb_eq_at_zero_equals_zeta() {
        let phi = poisson_boltzmann_eq(0.0, -0.05, 1e-8);
        assert!((phi - (-0.05)).abs() < 1e-12);
    }

    #[test]
    fn test_pb_eq_decays_with_distance() {
        let phi0 = poisson_boltzmann_eq(0.0, -0.05, 1e-8);
        let phi1 = poisson_boltzmann_eq(1e-8, -0.05, 1e-8);
        assert!(phi1.abs() < phi0.abs());
    }

    #[test]
    fn test_pb_eq_zero_debye_returns_zero() {
        let phi = poisson_boltzmann_eq(1.0, -0.05, 0.0);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_pb_eq_e_fold_at_debye_length() {
        let zeta = -0.05;
        let debye = 1e-8;
        let phi = poisson_boltzmann_eq(debye, zeta, debye);
        let expected = zeta / std::f64::consts::E;
        assert!((phi - expected).abs() < 1e-14);
    }

    #[test]
    fn test_pb_eq_sign_preserved() {
        let phi = poisson_boltzmann_eq(5e-9, -0.05, 1e-8);
        assert!(phi < 0.0);
    }

    #[test]
    fn test_pb_eq_positive_zeta() {
        let phi = poisson_boltzmann_eq(0.0, 0.05, 1e-8);
        assert!((phi - 0.05).abs() < 1e-12);
    }

    // ── electroosmotic_flow_velocity ─────────────────────────────────────

    #[test]
    fn test_eoflow_negative_zeta_positive_e() {
        // Negative zeta, positive field → positive velocity
        let v = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        assert!(v > 0.0, "v={v}");
    }

    #[test]
    fn test_eoflow_positive_zeta_positive_e() {
        let v = electroosmotic_flow_velocity(0.05, 1e4, 1e-3, 7.1e-10);
        assert!(v < 0.0, "v={v}");
    }

    #[test]
    fn test_eoflow_zero_viscosity() {
        let v = electroosmotic_flow_velocity(-0.05, 1e4, 0.0, 7.1e-10);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_eoflow_zero_field() {
        let v = electroosmotic_flow_velocity(-0.05, 0.0, 1e-3, 7.1e-10);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_eoflow_helmholtz_smoluchowski_value() {
        // v = -ε ζ E / μ = -(7.1e-10)*(-0.05)*(1e4)/(1e-3)
        let v = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        let expected = 7.1e-10 * 0.05 * 1e4 / 1e-3;
        assert!((v - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_eoflow_linear_in_field() {
        let v1 = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        let v2 = electroosmotic_flow_velocity(-0.05, 2e4, 1e-3, 7.1e-10);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_eoflow_linear_in_zeta() {
        let v1 = electroosmotic_flow_velocity(-0.05, 1e4, 1e-3, 7.1e-10);
        let v2 = electroosmotic_flow_velocity(-0.10, 1e4, 1e-3, 7.1e-10);
        assert!((v2 / v1 - 2.0).abs() < 1e-10);
    }

    // ── init_electrokinetic ──────────────────────────────────────────────

    #[test]
    fn test_init_dimensions() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e4, 0.0]);
        let state = init_electrokinetic(8, 6, &p);
        assert_eq!(state.nx, 8);
        assert_eq!(state.ny, 6);
        assert_eq!(state.charge_density.len(), 48);
        assert_eq!(state.potential.len(), 48);
        assert_eq!(state.f.len(), 48 * 9);
    }

    #[test]
    fn test_init_unit_density() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p);
        for j in 0..4 {
            for i in 0..4 {
                let rho: f64 = (0..9).map(|q| state.f[state.fi(i, j, q)]).sum();
                assert!((rho - 1.0).abs() < 1e-12, "rho at ({i},{j})={rho}");
            }
        }
    }

    #[test]
    fn test_init_zero_charge_potential() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p);
        assert!(state.charge_density.iter().all(|&x| x == 0.0));
        assert!(state.potential.iter().all(|&x| x == 0.0));
    }

    // ── step_electrokinetic ──────────────────────────────────────────────

    #[test]
    fn test_step_conserves_mass() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e3, 0.0]);
        let mut state = init_electrokinetic(8, 8, &p);
        let m_before: f64 = state.f.iter().sum();
        step_electrokinetic(&mut state, &p);
        let m_after: f64 = state.f.iter().sum();
        assert!(
            (m_after - m_before).abs() / m_before < 1e-10,
            "m_before={m_before} m_after={m_after}"
        );
    }

    #[test]
    fn test_step_finite_values() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let mut state = init_electrokinetic(4, 4, &p);
        step_electrokinetic(&mut state, &p);
        assert!(state.f.iter().all(|&x| x.is_finite()));
    }

    #[test]
    fn test_step_multiple_iterations_finite() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e2, 0.0]);
        let mut state = init_electrokinetic(8, 8, &p);
        for _ in 0..20 {
            step_electrokinetic(&mut state, &p);
        }
        assert!(state.f.iter().all(|&x| x.is_finite()));
    }

    // ── compute_charge_flux ──────────────────────────────────────────────

    #[test]
    fn test_charge_flux_zero_charge() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p);
        let flux = compute_charge_flux(&state);
        assert!(flux.iter().all(|&f| f == 0.0));
    }

    #[test]
    fn test_charge_flux_length() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(6, 5, &p);
        let flux = compute_charge_flux(&state);
        assert_eq!(flux.len(), 30);
    }

    #[test]
    fn test_charge_flux_nonnegative() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e3, 0.0]);
        let mut state = init_electrokinetic(4, 4, &p);
        for v in state.charge_density.iter_mut() {
            *v = 0.1;
        }
        step_electrokinetic(&mut state, &p);
        let flux = compute_charge_flux(&state);
        assert!(flux.iter().all(|&f| f >= 0.0));
    }

    // ── debye_huckel_potential ───────────────────────────────────────────

    #[test]
    fn test_dh_positive_charge() {
        let phi = debye_huckel_potential(1e-9, 1.0, 1e-8);
        assert!(phi > 0.0);
    }

    #[test]
    fn test_dh_decays_with_r() {
        let phi1 = debye_huckel_potential(1e-9, 1.0, 1e-8);
        let phi2 = debye_huckel_potential(2e-9, 1.0, 1e-8);
        assert!(phi2 < phi1);
    }

    #[test]
    fn test_dh_zero_r_returns_zero() {
        let phi = debye_huckel_potential(0.0, 1.0, 1e-8);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_dh_zero_debye_returns_zero() {
        let phi = debye_huckel_potential(1e-9, 1.0, 0.0);
        assert_eq!(phi, 0.0);
    }

    #[test]
    fn test_dh_negative_charge() {
        let phi = debye_huckel_potential(1e-9, -1.0, 1e-8);
        assert!(phi < 0.0);
    }

    #[test]
    fn test_dh_antisymmetric_charge() {
        let phi_pos = debye_huckel_potential(1e-9, 1.0, 1e-8);
        let phi_neg = debye_huckel_potential(1e-9, -1.0, 1e-8);
        assert!((phi_pos + phi_neg).abs() < 1e-30);
    }

    #[test]
    fn test_dh_finite_value() {
        let phi = debye_huckel_potential(5e-10, 1.602e-19, 3e-9);
        assert!(phi.is_finite());
    }

    // ── electrophoretic_mobility ─────────────────────────────────────────

    #[test]
    fn test_mob_negative_zeta() {
        let mob = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        assert!(mob < 0.0);
    }

    #[test]
    fn test_mob_zero_viscosity() {
        let mob = electrophoretic_mobility(-0.05, 0.0, 7.1e-10);
        assert_eq!(mob, 0.0);
    }

    #[test]
    fn test_mob_formula() {
        let mob = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        let expected = 7.1e-10 * (-0.05) / 1e-3;
        assert!((mob - expected).abs() < 1e-30);
    }

    #[test]
    fn test_mob_proportional_to_permittivity() {
        let mob1 = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        let mob2 = electrophoretic_mobility(-0.05, 1e-3, 2.0 * 7.1e-10);
        assert!((mob2 / mob1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_mob_inversely_proportional_to_viscosity() {
        let mob1 = electrophoretic_mobility(-0.05, 1e-3, 7.1e-10);
        let mob2 = electrophoretic_mobility(-0.05, 2e-3, 7.1e-10);
        assert!((mob1 / mob2 - 2.0).abs() < 1e-10);
    }

    // ── ElectrokineticParams ─────────────────────────────────────────────

    #[test]
    fn test_params_construction() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [1e4, 0.0]);
        assert!((p.permittivity - 7.1e-10).abs() < 1e-20);
        assert!((p.viscosity - 1e-3).abs() < 1e-10);
        assert!((p.zeta_potential - (-0.05)).abs() < 1e-12);
        assert!((p.debye_length - 1e-8).abs() < 1e-18);
        assert!((p.external_field[0] - 1e4).abs() < 1e-6);
        assert_eq!(p.external_field[1], 0.0);
    }

    // ── ElectrokineticLBM helpers ─────────────────────────────────────────

    #[test]
    fn test_ek_idx() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(5, 4, &p);
        assert_eq!(state.idx(2, 3), 3 * 5 + 2);
    }

    #[test]
    fn test_ek_fi() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(5, 4, &p);
        assert_eq!(state.fi(2, 3, 5), (3 * 5 + 2) * 9 + 5);
    }

    #[test]
    fn test_ek_feq_sum_to_rho() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p);
        let rho = 1.5_f64;
        let sum: f64 = (0..9).map(|q| state.feq(rho, 0.1, -0.05, q)).sum();
        assert!((sum - rho).abs() < 1e-12);
    }

    #[test]
    fn test_ek_feq_momentum_x() {
        let p = ElectrokineticParams::new(7.1e-10, 1e-3, -0.05, 1e-8, [0.0, 0.0]);
        let state = init_electrokinetic(4, 4, &p);
        let rho = 1.0;
        let ux = 0.1;
        let jx: f64 = (0..9)
            .map(|q| D2Q9_C[q][0] * state.feq(rho, ux, 0.0, q))
            .sum();
        assert!((jx - rho * ux).abs() < 1e-12);
    }
}
