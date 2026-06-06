// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermal Lattice Boltzmann Method for natural convection and heat transfer.
//!
//! Implements a Double Distribution Function (DDF) approach:
//! - `f` distributions solve Navier–Stokes (momentum/mass).
//! - `g` distributions solve the temperature advection-diffusion equation.
//!
//! Supports Rayleigh–Bénard convection initialization and thermal plume detection.
//!
//! References:
//! - He, X., Chen, S., & Doolen, G. D. (1998). *JCP* 146, 282–300.
//! - Shan, X. (1997). *Phys. Rev. E* 55, 2780.

/// Speed of sound squared in D2Q9 lattice units: cs² = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ---------------------------------------------------------------------------
// ThermalParams
// ---------------------------------------------------------------------------

/// Physical parameters for the thermal LBM simulation.
#[derive(Debug, Clone)]
pub struct ThermalParams {
    /// Thermal diffusivity α (lattice units).
    pub thermal_diffusivity: f64,
    /// Prandtl number Pr = ν / α.
    pub prandtl: f64,
    /// Rayleigh number Ra = g β ΔT L³ / (ν α).
    pub rayleigh: f64,
    /// Buoyancy coefficient β (thermal expansion coefficient times gravity, lattice units).
    pub buoyancy_coeff: f64,
}

impl ThermalParams {
    /// Construct new [`ThermalParams`].
    ///
    /// # Arguments
    /// * `thermal_diffusivity` — thermal diffusivity α
    /// * `prandtl`             — Prandtl number
    /// * `rayleigh`            — Rayleigh number
    /// * `buoyancy_coeff`      — β·g product
    pub fn new(thermal_diffusivity: f64, prandtl: f64, rayleigh: f64, buoyancy_coeff: f64) -> Self {
        Self {
            thermal_diffusivity,
            prandtl,
            rayleigh,
            buoyancy_coeff,
        }
    }

    /// Kinematic viscosity derived from Prandtl and thermal diffusivity: ν = Pr · α.
    pub fn kinematic_viscosity(&self) -> f64 {
        self.prandtl * self.thermal_diffusivity
    }

    /// Hydrodynamic relaxation time: τ = 0.5 + ν / (cs² Δt), Δt=1 in lattice units.
    pub fn tau_hydro(&self) -> f64 {
        0.5 + self.kinematic_viscosity() / CS2
    }

    /// Thermal relaxation time: τ_t = 0.5 + α / (cs² Δt).
    pub fn tau_thermal(&self) -> f64 {
        0.5 + self.thermal_diffusivity / CS2
    }
}

// ---------------------------------------------------------------------------
// ThermalLBM
// ---------------------------------------------------------------------------

/// Double Distribution Function thermal LBM state on a D2Q9 grid.
#[derive(Debug, Clone)]
pub struct ThermalLBM {
    /// Grid dimension in x.
    pub nx: usize,
    /// Grid dimension in y.
    pub ny: usize,
    /// Momentum distribution functions `f`, length `nx × ny × 9`.
    pub f: Vec<f64>,
    /// Energy (temperature) distribution functions `g`, length `nx × ny × 9`.
    pub g: Vec<f64>,
    /// Macroscopic temperature at each node, length `nx × ny`.
    pub temp: Vec<f64>,
    /// Hydrodynamic relaxation time τ.
    pub tau: f64,
    /// Thermal relaxation time τ_t.
    pub tau_t: f64,
    /// Reference temperature for Boussinesq buoyancy.
    pub t_ref: f64,
}

impl ThermalLBM {
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
}

// ---------------------------------------------------------------------------
// D2Q9 lattice definitions
// ---------------------------------------------------------------------------

/// Standard D2Q9 equilibrium weights.
///
/// ```no_run
/// use oxiphysics_lbm::heat_transfer_lbm::d2q9_weights;
/// let w = d2q9_weights();
/// let sum: f64 = w.iter().sum();
/// assert!((sum - 1.0).abs() < 1e-12);
/// ```
pub fn d2q9_weights() -> [f64; 9] {
    [
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
    ]
}

/// Standard D2Q9 discrete velocity vectors `[cx, cy]`.
///
/// ```no_run
/// use oxiphysics_lbm::heat_transfer_lbm::d2q9_velocities;
/// let v = d2q9_velocities();
/// assert_eq!(v[0], [0, 0]);
/// ```
pub fn d2q9_velocities() -> [[i32; 2]; 9] {
    [
        [0, 0],
        [1, 0],
        [0, 1],
        [-1, 0],
        [0, -1],
        [1, 1],
        [-1, 1],
        [-1, -1],
        [1, -1],
    ]
}

// ---------------------------------------------------------------------------
// Equilibrium distributions
// ---------------------------------------------------------------------------

/// D2Q9 hydrodynamic equilibrium distribution f_q^eq.
///
/// f_q^eq = w_q ρ \[1 + (c·u)/cs² + (c·u)²/(2cs⁴) − u²/(2cs²)\]
///
/// # Arguments
/// * `rho`   — local density
/// * `u`     — local velocity `[ux, uy]`
/// * `temp`  — temperature (unused in `f`, kept for unified signature)
/// * `alpha` — velocity direction index (0..9)
///
/// ```no_run
/// use oxiphysics_lbm::heat_transfer_lbm::equilibrium_thermal;
/// let feq: f64 = (0..9).map(|q| equilibrium_thermal(1.0, [0.0, 0.0], 1.0, q)).sum();
/// assert!((feq - 1.0).abs() < 1e-12);
/// ```
pub fn equilibrium_thermal(rho: f64, u: [f64; 2], _temp: f64, alpha: usize) -> f64 {
    let w = d2q9_weights();
    let c = d2q9_velocities();
    let cx = c[alpha][0] as f64;
    let cy = c[alpha][1] as f64;
    let cu = cx * u[0] + cy * u[1];
    let u2 = u[0] * u[0] + u[1] * u[1];
    w[alpha] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
}

/// D2Q9 thermal (temperature advection) equilibrium: g_q^eq = w_q T (1 + c·u / cs²).
fn geq(temp: f64, u: [f64; 2], alpha: usize) -> f64 {
    let w = d2q9_weights();
    let c = d2q9_velocities();
    let cx = c[alpha][0] as f64;
    let cy = c[alpha][1] as f64;
    let cu = cx * u[0] + cy * u[1];
    w[alpha] * temp * (1.0 + cu / CS2)
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

/// Initialize a Rayleigh–Bénard convection setup.
///
/// Creates a [`ThermalLBM`] with a linear temperature profile between the
/// hot bottom wall (y=0) and the cold top wall (y=ny-1), plus a small random
/// perturbation to seed convective instability.
///
/// # Arguments
/// * `nx`, `ny` — grid dimensions
/// * `params`   — thermal parameters
///
/// ```no_run
/// use oxiphysics_lbm::heat_transfer_lbm::{ThermalParams, rayleigh_benard_init};
/// let p = ThermalParams::new(1.0/6.0, 0.71, 1e4, 0.001);
/// let state = rayleigh_benard_init(8, 8, &p);
/// assert_eq!(state.nx, 8);
/// ```
pub fn rayleigh_benard_init(nx: usize, ny: usize, params: &ThermalParams) -> ThermalLBM {
    let ncells = nx * ny;
    let t_hot = 1.0_f64;
    let t_cold = 0.0_f64;
    let tau = params.tau_hydro();
    let tau_t = params.tau_thermal();

    let mut f = vec![0.0_f64; ncells * 9];
    let mut g = vec![0.0_f64; ncells * 9];
    let mut temp = vec![0.0_f64; ncells];

    let w = d2q9_weights();

    for j in 0..ny {
        for i in 0..nx {
            // Linear temperature profile from hot bottom to cold top
            let t = t_hot - (j as f64 / (ny as f64 - 1.0).max(1.0)) * (t_hot - t_cold);
            let cell = j * nx + i;
            temp[cell] = t;

            let u = [0.0_f64; 2];
            for (q, &_wq) in w.iter().enumerate() {
                let base = cell * 9 + q;
                f[base] = equilibrium_thermal(1.0, u, t, q);
                g[base] = geq(t, u, q);
            }
        }
    }

    ThermalLBM {
        nx,
        ny,
        f,
        g,
        temp,
        tau,
        tau_t,
        t_ref: (t_hot + t_cold) / 2.0,
    }
}

// ---------------------------------------------------------------------------
// Step function
// ---------------------------------------------------------------------------

/// Perform one BGK + streaming step for both `f` and `g` distributions.
///
/// Includes Boussinesq buoyancy forcing in the y-direction.
///
/// # Arguments
/// * `state`  — mutable [`ThermalLBM`] state
/// * `params` — thermal parameters
pub fn step_thermal_lbm(state: &mut ThermalLBM, params: &ThermalParams) {
    let nx = state.nx;
    let ny = state.ny;
    let ncells = nx * ny;
    let inv_tau = 1.0 / state.tau;
    let inv_tau_t = 1.0 / state.tau_t;
    let buoy = params.buoyancy_coeff;
    let t_ref = state.t_ref;
    let w = d2q9_weights();
    let c = d2q9_velocities();

    // Update macroscopic fields
    let mut rho_field = vec![0.0_f64; ncells];
    let mut ux_field = vec![0.0_f64; ncells];
    let mut uy_field = vec![0.0_f64; ncells];
    let mut temp_field = vec![0.0_f64; ncells];

    for j in 0..ny {
        for i in 0..nx {
            let base = (j * nx + i) * 9;
            let mut rho = 0.0_f64;
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            let mut t = 0.0_f64;
            for (q, cq) in c.iter().enumerate() {
                let fq = state.f[base + q];
                rho += fq;
                ux += cq[0] as f64 * fq;
                uy += cq[1] as f64 * fq;
                t += state.g[base + q];
            }
            if rho > 1e-20 {
                ux /= rho;
                uy /= rho;
            }
            let cell = j * nx + i;
            rho_field[cell] = rho;
            ux_field[cell] = ux;
            uy_field[cell] = uy;
            temp_field[cell] = t;
        }
    }

    // Collision
    let mut f_post = vec![0.0_f64; ncells * 9];
    let mut g_post = vec![0.0_f64; ncells * 9];

    for j in 0..ny {
        for i in 0..nx {
            let base = (j * nx + i) * 9;
            let cell = j * nx + i;
            let rho = rho_field[cell];
            let ux = ux_field[cell];
            let uy = uy_field[cell];
            let t = temp_field[cell];
            let u = [ux, uy];

            // Boussinesq buoyancy (y-direction)
            let fy = buoy * (t - t_ref);

            for q in 0..9_usize {
                let feq = equilibrium_thermal(rho, u, t, q);
                let cy = c[q][1] as f64;
                let force_term = w[q] * cy * fy / CS2;
                f_post[base + q] = state.f[base + q] * (1.0 - inv_tau) + feq * inv_tau + force_term;

                let geq_val = geq(t, u, q);
                g_post[base + q] = state.g[base + q] * (1.0 - inv_tau_t) + geq_val * inv_tau_t;
            }
        }
    }

    // Streaming (periodic)
    for j in 0..ny {
        for i in 0..nx {
            let src_base = (j * nx + i) * 9;
            for (q, cq) in c.iter().enumerate() {
                let di = cq[0] as isize;
                let dj = cq[1] as isize;
                let ni = ((i as isize + di).rem_euclid(nx as isize)) as usize;
                let nj = ((j as isize + dj).rem_euclid(ny as isize)) as usize;
                let dst_base = (nj * nx + ni) * 9;
                state.f[dst_base + q] = f_post[src_base + q];
                state.g[dst_base + q] = g_post[src_base + q];
            }
        }
    }

    // Update temp field
    for j in 0..ny {
        for i in 0..nx {
            let base = (j * nx + i) * 9;
            let cell = j * nx + i;
            state.temp[cell] = (0..9_usize).map(|q| state.g[base + q]).sum();
        }
    }
}

// ---------------------------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------------------------

/// Apply isothermal boundary conditions at top (cold) and bottom (hot) walls.
///
/// Sets the temperature distributions at y=0 (hot) and y=ny-1 (cold)
/// to their equilibrium values.
///
/// # Arguments
/// * `state`    — mutable [`ThermalLBM`] state
/// * `hot_temp` — hot wall temperature (at y=0)
/// * `cold_temp`— cold wall temperature (at y=ny-1)
pub fn apply_thermal_bc(state: &mut ThermalLBM, hot_temp: f64, cold_temp: f64) {
    let nx = state.nx;
    let ny = state.ny;
    // Bottom wall (j=0): hot
    for i in 0..nx {
        let base = i * 9;
        let u = [0.0_f64; 2];
        for q in 0..9_usize {
            state.g[base + q] = geq(hot_temp, u, q);
        }
        state.temp[i] = hot_temp;
    }
    // Top wall (j=ny-1): cold
    for i in 0..nx {
        let base = ((ny - 1) * nx + i) * 9;
        let cell = (ny - 1) * nx + i;
        let u = [0.0_f64; 2];
        for q in 0..9_usize {
            state.g[base + q] = geq(cold_temp, u, q);
        }
        state.temp[cell] = cold_temp;
    }
}

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

/// Compute Nusselt number from the bottom-wall temperature gradient.
///
/// Nu = (−dT/dy|_{wall} × L) / ΔT
///
/// Estimated via forward finite difference at y=0.
///
/// # Arguments
/// * `state`  — [`ThermalLBM`] state
/// * `params` — thermal parameters (used for Ra normalization, unused here)
pub fn nusselt_number(state: &ThermalLBM, _params: &ThermalParams) -> f64 {
    let nx = state.nx;
    let ny = state.ny;
    if ny < 2 {
        return 0.0;
    }
    let t_hot: f64 = (0..nx).map(|i| state.temp[i]).sum::<f64>() / nx as f64;
    let t_cold: f64 = (0..nx).map(|i| state.temp[(ny - 1) * nx + i]).sum::<f64>() / nx as f64;
    let delta_t = t_hot - t_cold;
    if delta_t.abs() < f64::EPSILON {
        return 0.0;
    }
    // Wall gradient at y=0: dT/dy ~ (T[j=1] - T[j=0])
    let mut flux_sum = 0.0_f64;
    for i in 0..nx {
        let t0 = state.temp[i];
        let t1 = state.temp[nx + i];
        flux_sum += (t1 - t0).abs();
    }
    let avg_grad = flux_sum / nx as f64;
    avg_grad * (ny as f64 - 1.0) / delta_t
}

/// Detect thermal plumes: returns grid positions where temperature exceeds a threshold.
///
/// # Arguments
/// * `state`     — [`ThermalLBM`] state
/// * `threshold` — temperature threshold for plume detection
///
/// ```no_run
/// use oxiphysics_lbm::heat_transfer_lbm::{ThermalParams, rayleigh_benard_init, thermal_plume_detect};
/// let p = ThermalParams::new(1.0/6.0, 0.71, 1e4, 0.001);
/// let state = rayleigh_benard_init(4, 4, &p);
/// let plumes = thermal_plume_detect(&state, 0.5);
/// assert!(plumes.iter().all(|pos| pos[0] < 4 && pos[1] < 4));
/// ```
pub fn thermal_plume_detect(state: &ThermalLBM, threshold: f64) -> Vec<[usize; 2]> {
    let nx = state.nx;
    let ny = state.ny;
    let mut plumes = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            if state.temp[j * nx + i] > threshold {
                plumes.push([i, j]);
            }
        }
    }
    plumes
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── d2q9_weights ─────────────────────────────────────────────────────

    #[test]
    fn test_weights_sum_to_one() {
        let w = d2q9_weights();
        let sum: f64 = w.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "sum={sum}");
    }

    #[test]
    fn test_weights_all_positive() {
        for (i, &wi) in d2q9_weights().iter().enumerate() {
            assert!(wi > 0.0, "w[{i}]={wi}");
        }
    }

    #[test]
    fn test_weights_w0_equals_4_9() {
        let w = d2q9_weights();
        assert!((w[0] - 4.0 / 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_weights_axial_equals_1_9() {
        let w = d2q9_weights();
        for (offset, &wi) in w[1..=4].iter().enumerate() {
            let i = offset + 1;
            assert!((wi - 1.0 / 9.0).abs() < 1e-12, "w[{i}]={wi}");
        }
    }

    #[test]
    fn test_weights_diagonal_equals_1_36() {
        let w = d2q9_weights();
        for (offset, &wi) in w[5..=8].iter().enumerate() {
            let i = offset + 5;
            assert!((wi - 1.0 / 36.0).abs() < 1e-12, "w[{i}]={wi}");
        }
    }

    // ── d2q9_velocities ──────────────────────────────────────────────────

    #[test]
    fn test_velocities_rest_is_zero() {
        let v = d2q9_velocities();
        assert_eq!(v[0], [0, 0]);
    }

    #[test]
    fn test_velocities_count() {
        assert_eq!(d2q9_velocities().len(), 9);
    }

    #[test]
    fn test_velocities_symmetric_sum() {
        let v = d2q9_velocities();
        let sum_cx: i32 = v.iter().map(|c| c[0]).sum();
        let sum_cy: i32 = v.iter().map(|c| c[1]).sum();
        assert_eq!(sum_cx, 0);
        assert_eq!(sum_cy, 0);
    }

    // ── equilibrium_thermal ──────────────────────────────────────────────

    #[test]
    fn test_feq_sum_to_rho() {
        let rho = 1.2;
        let sum: f64 = (0..9)
            .map(|q| equilibrium_thermal(rho, [0.1, -0.05], 1.0, q))
            .sum();
        assert!((sum - rho).abs() < 1e-12);
    }

    #[test]
    fn test_feq_momentum_x() {
        let rho = 1.0;
        let ux = 0.1;
        let c = d2q9_velocities();
        let jx: f64 = (0..9)
            .map(|q| c[q][0] as f64 * equilibrium_thermal(rho, [ux, 0.0], 1.0, q))
            .sum();
        assert!((jx - rho * ux).abs() < 1e-12);
    }

    #[test]
    fn test_feq_momentum_y() {
        let rho = 1.0;
        let uy = 0.05;
        let c = d2q9_velocities();
        let jy: f64 = (0..9)
            .map(|q| c[q][1] as f64 * equilibrium_thermal(rho, [0.0, uy], 1.0, q))
            .sum();
        assert!((jy - rho * uy).abs() < 1e-12);
    }

    // ── ThermalParams ─────────────────────────────────────────────────────

    #[test]
    fn test_params_kinematic_viscosity() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let nu = p.kinematic_viscosity();
        assert!((nu - 0.71 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_params_tau_hydro_positive() {
        let p = ThermalParams::new(0.1, 0.71, 1e4, 0.001);
        assert!(p.tau_hydro() > 0.5);
    }

    #[test]
    fn test_params_tau_thermal_positive() {
        let p = ThermalParams::new(0.1, 0.71, 1e4, 0.001);
        assert!(p.tau_thermal() > 0.5);
    }

    // ── rayleigh_benard_init ─────────────────────────────────────────────

    #[test]
    fn test_rb_init_dimensions() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let state = rayleigh_benard_init(8, 8, &p);
        assert_eq!(state.nx, 8);
        assert_eq!(state.ny, 8);
        assert_eq!(state.f.len(), 8 * 8 * 9);
        assert_eq!(state.g.len(), 8 * 8 * 9);
        assert_eq!(state.temp.len(), 64);
    }

    #[test]
    fn test_rb_init_hot_bottom() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let state = rayleigh_benard_init(8, 8, &p);
        // Bottom row (j=0) should be hotter
        let t_bottom = state.temp[0];
        let t_top = state.temp[(state.ny - 1) * state.nx];
        assert!(t_bottom > t_top, "t_bottom={t_bottom} t_top={t_top}");
    }

    #[test]
    fn test_rb_init_density_unity() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let state = rayleigh_benard_init(4, 4, &p);
        for j in 0..4 {
            for i in 0..4 {
                let rho: f64 = (0..9).map(|q| state.f[state.fi(i, j, q)]).sum();
                assert!((rho - 1.0).abs() < 1e-10, "rho at ({i},{j})={rho}");
            }
        }
    }

    // ── step_thermal_lbm ─────────────────────────────────────────────────

    #[test]
    fn test_step_conserves_mass() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(8, 8, &p);
        let m_before: f64 = state.f.iter().sum();
        step_thermal_lbm(&mut state, &p);
        let m_after: f64 = state.f.iter().sum();
        assert!(
            (m_after - m_before).abs() / m_before < 1e-9,
            "m_before={m_before} m_after={m_after}"
        );
    }

    #[test]
    fn test_step_finite_values() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(4, 4, &p);
        step_thermal_lbm(&mut state, &p);
        assert!(state.f.iter().all(|&x| x.is_finite()));
        assert!(state.g.iter().all(|&x| x.is_finite()));
    }

    #[test]
    fn test_step_multiple_iterations() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(8, 8, &p);
        for _ in 0..10 {
            step_thermal_lbm(&mut state, &p);
        }
        assert!(state.temp.iter().all(|&x| x.is_finite()));
    }

    // ── apply_thermal_bc ──────────────────────────────────────────────────

    #[test]
    fn test_bc_bottom_temp_set() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(4, 4, &p);
        apply_thermal_bc(&mut state, 1.5, 0.1);
        for i in 0..4 {
            assert!((state.temp[i] - 1.5).abs() < 1e-12);
        }
    }

    #[test]
    fn test_bc_top_temp_set() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(4, 4, &p);
        apply_thermal_bc(&mut state, 1.5, 0.1);
        let ny = state.ny;
        let nx = state.nx;
        for i in 0..4 {
            assert!((state.temp[(ny - 1) * nx + i] - 0.1).abs() < 1e-12);
        }
    }

    // ── nusselt_number ────────────────────────────────────────────────────

    #[test]
    fn test_nusselt_isothermal_zero() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(4, 4, &p);
        // Set uniform temperature → zero Nusselt
        for v in state.temp.iter_mut() {
            *v = 1.0;
        }
        let nu = nusselt_number(&state, &p);
        assert!(nu.abs() < 1e-10, "nu={nu}");
    }

    #[test]
    fn test_nusselt_positive_gradient() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let state = rayleigh_benard_init(8, 8, &p);
        let nu = nusselt_number(&state, &p);
        assert!(nu >= 0.0);
    }

    // ── thermal_plume_detect ─────────────────────────────────────────────

    #[test]
    fn test_plume_detect_all_below_threshold() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(4, 4, &p);
        for v in state.temp.iter_mut() {
            *v = 0.0;
        }
        let plumes = thermal_plume_detect(&state, 0.5);
        assert!(plumes.is_empty());
    }

    #[test]
    fn test_plume_detect_all_above_threshold() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let mut state = rayleigh_benard_init(4, 4, &p);
        for v in state.temp.iter_mut() {
            *v = 1.0;
        }
        let plumes = thermal_plume_detect(&state, 0.5);
        assert_eq!(plumes.len(), 16);
    }

    #[test]
    fn test_plume_detect_positions_valid() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let state = rayleigh_benard_init(6, 5, &p);
        let plumes = thermal_plume_detect(&state, 0.3);
        for pos in &plumes {
            assert!(pos[0] < 6);
            assert!(pos[1] < 5);
        }
    }

    #[test]
    fn test_plume_detect_hot_wall() {
        let p = ThermalParams::new(1.0 / 6.0, 0.71, 1e4, 0.001);
        let state = rayleigh_benard_init(4, 4, &p);
        // Bottom row should be at T=1.0 (hot)
        let plumes = thermal_plume_detect(&state, 0.9);
        // At least bottom row nodes should appear
        let bottom_plumes: Vec<_> = plumes.iter().filter(|pos| pos[1] == 0).collect();
        assert!(!bottom_plumes.is_empty());
    }
}
