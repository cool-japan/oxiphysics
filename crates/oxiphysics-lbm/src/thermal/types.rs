//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};

use super::functions::*;

/// Configuration for the de Vahl Davis (1983) natural convection benchmark.
///
/// Square cavity: height = width = `n`.  Left wall hot, right wall cold,
/// top and bottom insulated.
///
/// Reference: de Vahl Davis, G. (1983). Natural convection of air in a square
/// cavity: A benchmark numerical solution. *Int. J. Numer. Methods Fluids* 3,
/// 249–264.
#[allow(dead_code)]
pub struct DeVahlDavisSetup {
    /// Number of lattice nodes per side (square cavity).
    pub n: usize,
    /// Hot wall temperature (left wall).
    pub t_hot: f64,
    /// Cold wall temperature (right wall).
    pub t_cold: f64,
    /// Rayleigh number.
    pub ra: f64,
    /// Prandtl number (air ≈ 0.71).
    pub pr: f64,
}
#[allow(dead_code)]
impl DeVahlDavisSetup {
    /// Create a new de Vahl Davis cavity setup.
    pub fn new(n: usize, t_hot: f64, t_cold: f64, ra: f64, pr: f64) -> Self {
        Self {
            n,
            t_hot,
            t_cold,
            ra,
            pr,
        }
    }
    /// LBM transport parameters derived from Ra and Pr.
    ///
    /// Returns `(tau_f, tau_t, beta, gravity)` where:
    /// - `tau_f` = fluid relaxation time
    /// - `tau_t` = thermal relaxation time
    /// - `beta`  = thermal expansion coefficient
    /// - `gravity` = gravitational body force vector `[0, -g_lbm]`
    pub fn lbm_parameters(&self) -> (f64, f64, f64, [f64; 2]) {
        let l = self.n as f64;
        let u_ref = 1e-2_f64;
        let g_lbm = 1e-4_f64;
        let nu = (self.pr / self.ra).sqrt() * u_ref * l;
        let alpha = nu / self.pr;
        let tau_f = 3.0 * nu + 0.5;
        let tau_t = 3.0 * alpha + 0.5;
        let delta_t = (self.t_hot - self.t_cold).abs().max(1e-30);
        let beta = self.ra * nu * alpha / (g_lbm * delta_t * l * l * l);
        (tau_f, tau_t, beta, [0.0, -g_lbm])
    }
    /// Initialize a linear temperature profile from hot (left) to cold (right).
    ///
    /// `temperature[y * n + x] = T_hot + (T_cold - T_hot) * x / (n-1)`
    pub fn initial_temperature(&self) -> Vec<f64> {
        let n = self.n;
        let mut temperature = vec![0.0; n * n];
        for y in 0..n {
            for x in 0..n {
                let frac = x as f64 / (n - 1).max(1) as f64;
                temperature[y * n + x] = self.t_hot + (self.t_cold - self.t_hot) * frac;
            }
        }
        temperature
    }
    /// Apply hot-wall BC on left column (`x = 0`).
    pub fn apply_hot_wall(&self, thermal: &mut ThermalD2Q9, velocities: &[[f64; 2]]) {
        thermal.apply_temperature_bc_left(self.t_hot, velocities);
    }
    /// Apply cold-wall BC on right column (`x = n-1`).
    pub fn apply_cold_wall(&self, thermal: &mut ThermalD2Q9, velocities: &[[f64; 2]]) {
        thermal.apply_temperature_bc_right(self.t_cold, velocities);
    }
    /// Apply insulating BCs on top and bottom walls using bounce-back.
    pub fn apply_insulating_walls(&self, thermal: &mut ThermalD2Q9) {
        thermal.apply_insulated_bc_bottom();
    }
    /// Nusselt number from the average temperature gradient on the hot wall.
    ///
    /// Uses a one-sided finite difference in x at `x = 0`.
    pub fn nusselt_hot_wall(&self, temperature: &[f64]) -> f64 {
        de_vahl_davis_nusselt(temperature, self.n, self.n, self.t_hot, self.t_cold)
    }
    /// Check whether the LBM parameters satisfy basic stability requirements.
    ///
    /// Returns `true` if `0.5 < τ_f < 2.0` and `0.5 < τ_t < 2.0`.
    pub fn check_stability(&self) -> bool {
        let (tau_f, tau_t, _beta, _g) = self.lbm_parameters();
        tau_f > 0.5 && tau_f < 2.0 && tau_t > 0.5 && tau_t < 2.0
    }
    /// Dimensionless temperature: `θ = (T - T_cold) / (T_hot - T_cold)`.
    ///
    /// Returns a field of θ ∈ \[0, 1\].
    pub fn dimensionless_temperature(&self, temperature: &[f64]) -> Vec<f64> {
        let delta_t = (self.t_hot - self.t_cold).abs().max(1e-30);
        temperature
            .iter()
            .map(|&t| (t - self.t_cold) / delta_t)
            .collect()
    }
    /// Horizontal mid-plane temperature profile at `y = n/2`.
    ///
    /// Returns a vector of `n` temperature values along the horizontal centre line.
    pub fn mid_plane_x_temperature(&self, temperature: &[f64]) -> Vec<f64> {
        let y_mid = self.n / 2;
        (0..self.n)
            .map(|x| temperature[y_mid * self.n + x])
            .collect()
    }
    /// Vertical mid-plane temperature profile at `x = n/2`.
    pub fn mid_plane_y_temperature(&self, temperature: &[f64]) -> Vec<f64> {
        let x_mid = self.n / 2;
        (0..self.n)
            .map(|y| temperature[y * self.n + x_mid])
            .collect()
    }
}
/// Utilities for computing the Nusselt number from thermal LBM data.
#[allow(dead_code)]
pub struct NusseltComputation;
#[allow(dead_code)]
impl NusseltComputation {
    /// Compute local Nusselt number at the bottom wall from a 2D temperature
    /// field stored as `temperature[y * nx + x]`.
    ///
    /// Uses a one-sided finite difference at y=0:
    ///
    /// ```text
    /// Nu_local(x) = -dT/dy|_{y=0} * H / ΔT
    ///             = -(T(x,1) - T(x,0)) * H / ΔT
    /// ```
    pub fn local_nusselt_bottom(
        temperature: &[f64],
        nx: usize,
        ny: usize,
        t_hot: f64,
        t_cold: f64,
    ) -> Vec<f64> {
        let h = ny as f64;
        let delta_t = (t_hot - t_cold).abs().max(1e-30);
        let mut nu_local = Vec::with_capacity(nx);
        for x in 0..nx {
            let t0 = temperature[x];
            let t1 = temperature[nx + x];
            let grad = t1 - t0;
            nu_local.push((-grad).abs() * h / delta_t);
        }
        nu_local
    }
    /// Compute local Nusselt number at the top wall.
    pub fn local_nusselt_top(
        temperature: &[f64],
        nx: usize,
        ny: usize,
        t_hot: f64,
        t_cold: f64,
    ) -> Vec<f64> {
        let h = ny as f64;
        let delta_t = (t_hot - t_cold).abs().max(1e-30);
        let mut nu_local = Vec::with_capacity(nx);
        let y_top = ny - 1;
        for x in 0..nx {
            let t_top = temperature[y_top * nx + x];
            let t_below = temperature[(y_top - 1) * nx + x];
            let grad = t_top - t_below;
            nu_local.push((-grad).abs() * h / delta_t);
        }
        nu_local
    }
    /// Compute the average Nusselt number from local values.
    pub fn average_nusselt(local_nu: &[f64]) -> f64 {
        if local_nu.is_empty() {
            return 0.0;
        }
        local_nu.iter().sum::<f64>() / local_nu.len() as f64
    }
    /// Compute Nusselt number via volume-averaged temperature gradient.
    ///
    /// ```text
    /// Nu = 1 + <u_y * T> * H / (α * ΔT)
    /// ```
    ///
    /// where the average is over the entire domain.
    pub fn nusselt_volume_averaged(
        temperature: &[f64],
        uy: &[f64],
        nx: usize,
        ny: usize,
        alpha: f64,
        t_hot: f64,
        t_cold: f64,
    ) -> f64 {
        let h = ny as f64;
        let delta_t = (t_hot - t_cold).abs().max(1e-30);
        let n = nx * ny;
        let mut sum_uy_t = 0.0;
        for k in 0..n {
            sum_uy_t += uy[k] * temperature[k];
        }
        let avg_uy_t = sum_uy_t / n as f64;
        1.0 + avg_uy_t * h / (alpha * delta_t)
    }
}
/// Configuration for natural convection in a differentially heated cavity.
///
/// Left wall hot, right wall cold, top and bottom insulated.
#[allow(dead_code)]
pub struct NaturalConvectionSetup {
    /// Number of lattice nodes in x.
    pub nx: usize,
    /// Number of lattice nodes in y.
    pub ny: usize,
    /// Hot wall (left) temperature.
    pub t_hot: f64,
    /// Cold wall (right) temperature.
    pub t_cold: f64,
    /// Rayleigh number.
    pub ra: f64,
    /// Prandtl number.
    pub pr: f64,
}
#[allow(dead_code)]
impl NaturalConvectionSetup {
    /// Create a new natural convection setup.
    pub fn new(nx: usize, ny: usize, t_hot: f64, t_cold: f64, ra: f64, pr: f64) -> Self {
        Self {
            nx,
            ny,
            t_hot,
            t_cold,
            ra,
            pr,
        }
    }
    /// Derive LBM parameters from Ra, Pr.
    ///
    /// Returns `(tau_f, tau_t, beta, gravity)`.
    pub fn lbm_parameters(&self) -> (f64, f64, f64, [f64; 2]) {
        let l = self.nx as f64;
        let u_ref = 1e-2_f64;
        let g_ref = 1e-4_f64;
        let nu = (self.pr / self.ra).sqrt() * u_ref * l;
        let alpha = nu / self.pr;
        let tau_f = 3.0 * nu + 0.5;
        let tau_t = 3.0 * alpha + 0.5;
        let delta_t = (self.t_hot - self.t_cold).abs().max(1e-30);
        let beta = self.ra * nu * alpha / (g_ref * delta_t * l * l * l);
        let gravity = [0.0, -g_ref];
        (tau_f, tau_t, beta, gravity)
    }
    /// Initialize a linear temperature profile from hot (left) to cold (right).
    pub fn initial_temperature(&self) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut temperature = vec![0.0; n];
        for y in 0..ny {
            for x in 0..nx {
                let k = y * nx + x;
                let frac = x as f64 / (nx - 1).max(1) as f64;
                temperature[k] = self.t_hot + (self.t_cold - self.t_hot) * frac;
            }
        }
        temperature
    }
}
/// D2Q9 thermal distribution function using the Double Distribution Function
/// (DDF) approach.
///
/// The thermal field is governed by the advection–diffusion equation.
/// The equilibrium distribution is:
///
/// ```text
/// g_eq_i = w_i * T * (1 + (e_i · u) / cs²)
/// ```
///
/// and the BGK collision reads:
///
/// ```text
/// g_post_i = g_i - (g_i - g_eq_i) / tau_t
/// ```
pub struct ThermalD2Q9 {
    /// Temperature distribution, `g[node_index][direction]`.
    pub g: Vec<[f64; 9]>,
    /// Macroscopic temperature at each node.
    pub temperature: Vec<f64>,
    /// Number of lattice nodes in x.
    pub nx: usize,
    /// Number of lattice nodes in y.
    pub ny: usize,
    /// Thermal relaxation time τ_t.
    pub tau_t: f64,
    /// Thermal diffusivity α = (τ_t − 0.5) / 3.
    pub alpha: f64,
}
impl ThermalD2Q9 {
    /// Create a new thermal field initialised to equilibrium at temperature
    /// `t_init` with zero velocity everywhere.
    pub fn new(nx: usize, ny: usize, tau_t: f64, t_init: f64) -> Self {
        let n = nx * ny;
        let alpha = (tau_t - 0.5) / 3.0;
        let g_eq: [f64; 9] = {
            let mut arr = [0.0; 9];
            for (i, w) in D2Q9_WEIGHTS.iter().enumerate() {
                arr[i] = w * t_init;
            }
            arr
        };
        Self {
            g: vec![g_eq; n],
            temperature: vec![t_init; n],
            nx,
            ny,
            tau_t,
            alpha,
        }
    }
    /// Flat node index for lattice position `(x, y)`.
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Equilibrium distribution for direction `i`:
    ///
    /// ```text
    /// g_eq_i = w_i * T * (1 + (e_i · u) / cs²)
    /// ```
    pub fn equilibrium(&self, temp: f64, u: [f64; 2], i: usize) -> f64 {
        let w = D2Q9_WEIGHTS[i];
        let c = D2Q9_VELOCITIES[i];
        let e_dot_u = c[0] as f64 * u[0] + c[1] as f64 * u[1];
        w * temp * (1.0 + e_dot_u / CS2)
    }
    /// Second-order equilibrium distribution with velocity squared terms.
    ///
    /// ```text
    /// g_eq_i = w_i * T * (1 + (e_i · u)/cs² + (e_i · u)²/(2cs⁴) − u²/(2cs²))
    /// ```
    #[allow(dead_code)]
    pub fn equilibrium_second_order(&self, temp: f64, u: [f64; 2], i: usize) -> f64 {
        let w = D2Q9_WEIGHTS[i];
        let c = D2Q9_VELOCITIES[i];
        let e_dot_u = c[0] as f64 * u[0] + c[1] as f64 * u[1];
        let u2 = u[0] * u[0] + u[1] * u[1];
        w * temp * (1.0 + e_dot_u / CS2 + e_dot_u * e_dot_u / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    }
    /// BGK collision: `g_post = g - (g - g_eq) / tau_t`.
    ///
    /// `velocities` is a flat slice of `[ux, uy]` per node (length = nx*ny).
    pub fn collide(&mut self, velocities: &[[f64; 2]]) {
        let n = self.nx * self.ny;
        for k in 0..n {
            let t = self.temperature[k];
            let u = velocities[k];
            for i in 0..9 {
                let geq = self.equilibrium(t, u, i);
                self.g[k][i] -= (self.g[k][i] - geq) / self.tau_t;
            }
        }
    }
    /// Streaming step with periodic boundaries in both directions.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut g_new = vec![[0.0_f64; 9]; n];
        for y in 0..ny {
            for x in 0..nx {
                let k_src = self.idx(x, y);
                for i in 0..9 {
                    let c = D2Q9_VELOCITIES[i];
                    let xd = ((x as isize + c[0] as isize).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + c[1] as isize).rem_euclid(ny as isize)) as usize;
                    let k_dst = self.idx(xd, yd);
                    g_new[k_dst][i] = self.g[k_src][i];
                }
            }
        }
        self.g = g_new;
    }
    /// Compute macroscopic temperature: `T = Σ_i g_i` (zeroth moment).
    pub fn compute_temperature(&mut self) {
        let n = self.nx * self.ny;
        for k in 0..n {
            self.temperature[k] = self.g[k].iter().sum();
        }
    }
    /// Apply fixed-temperature boundary condition at the bottom wall (`y = 0`)
    /// using anti-bounce-back (equilibrium re-initialisation).
    ///
    /// All nodes on the bottom row are set to equilibrium distributions
    /// corresponding to temperature `t_bot` and the local fluid velocity.
    pub fn apply_temperature_bc_bottom(&mut self, t_bot: f64, velocities: &[[f64; 2]]) {
        for x in 0..self.nx {
            let k = self.idx(x, 0);
            let u = velocities[k];
            for i in 0..9 {
                self.g[k][i] = self.equilibrium(t_bot, u, i);
            }
            self.temperature[k] = t_bot;
        }
    }
    /// Apply fixed-temperature boundary condition at the top wall (`y = ny-1`)
    /// using equilibrium re-initialisation.
    pub fn apply_temperature_bc_top(&mut self, t_top: f64, velocities: &[[f64; 2]]) {
        let y_top = self.ny - 1;
        for x in 0..self.nx {
            let k = self.idx(x, y_top);
            let u = velocities[k];
            for i in 0..9 {
                self.g[k][i] = self.equilibrium(t_top, u, i);
            }
            self.temperature[k] = t_top;
        }
    }
    /// Apply fixed-temperature BC on the left wall (`x = 0`).
    #[allow(dead_code)]
    pub fn apply_temperature_bc_left(&mut self, t_left: f64, velocities: &[[f64; 2]]) {
        for y in 0..self.ny {
            let k = self.idx(0, y);
            let u = velocities[k];
            for i in 0..9 {
                self.g[k][i] = self.equilibrium(t_left, u, i);
            }
            self.temperature[k] = t_left;
        }
    }
    /// Apply fixed-temperature BC on the right wall (`x = nx-1`).
    #[allow(dead_code)]
    pub fn apply_temperature_bc_right(&mut self, t_right: f64, velocities: &[[f64; 2]]) {
        let x_right = self.nx - 1;
        for y in 0..self.ny {
            let k = self.idx(x_right, y);
            let u = velocities[k];
            for i in 0..9 {
                self.g[k][i] = self.equilibrium(t_right, u, i);
            }
            self.temperature[k] = t_right;
        }
    }
    /// Apply insulated (zero-flux) boundary at the bottom wall using
    /// bounce-back on the thermal distributions.
    #[allow(dead_code)]
    pub fn apply_insulated_bc_bottom(&mut self) {
        const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
        for x in 0..self.nx {
            let k = self.idx(x, 0);
            let mut tmp = self.g[k];
            for i in 0..9 {
                self.g[k][i] = tmp[OPP[i]];
            }
            let _ = &mut tmp;
        }
    }
    /// One full thermal step: collide → stream → compute_temperature.
    pub fn step(&mut self, velocities: &[[f64; 2]]) {
        self.collide(velocities);
        self.stream();
        self.compute_temperature();
    }
    /// Compute the average temperature over the entire domain.
    #[allow(dead_code)]
    pub fn average_temperature(&self) -> f64 {
        let n = self.nx * self.ny;
        self.temperature.iter().sum::<f64>() / n as f64
    }
    /// Compute the temperature variance over the domain.
    #[allow(dead_code)]
    pub fn temperature_variance(&self) -> f64 {
        let avg = self.average_temperature();
        let n = self.nx * self.ny;
        self.temperature
            .iter()
            .map(|&t| (t - avg) * (t - avg))
            .sum::<f64>()
            / n as f64
    }
    /// Compute heat flux in x-direction at a node: q_x = Σ_i e_ix * g_i.
    #[allow(dead_code)]
    pub fn heat_flux_x(&self, k: usize) -> f64 {
        let mut qx = 0.0;
        for i in 0..9 {
            qx += D2Q9_VELOCITIES[i][0] as f64 * self.g[k][i];
        }
        qx
    }
    /// Compute heat flux in y-direction at a node: q_y = Σ_i e_iy * g_i.
    #[allow(dead_code)]
    pub fn heat_flux_y(&self, k: usize) -> f64 {
        let mut qy = 0.0;
        for i in 0..9 {
            qy += D2Q9_VELOCITIES[i][1] as f64 * self.g[k][i];
        }
        qy
    }
}
/// High-level combined thermal LBM simulation object using the
/// Double-Distribution-Function (DDF) approach.
///
/// Couples a D2Q9 flow solver with a D2Q9 thermal advection-diffusion solver.
/// The temperature field is advected by the flow velocity and diffused with
/// thermal diffusivity `α = (τ_t − 0.5) / 3`.
///
/// The Boussinesq approximation couples the temperature back to the flow via
/// buoyancy forces when `boussinesq` is provided.
///
/// # Example setup (Rayleigh-Bénard):
///
/// 1. Create `ThermalLbm::new(nx, ny, tau_t, t_init)`.
/// 2. Each time step: call `step_thermal(&velocities)`.
/// 3. Apply BCs: `apply_bc_bottom(T_hot)`, `apply_bc_top(T_cold)`.
/// 4. Read `thermal.temperature[k]` for coupling to the flow solver.
#[allow(dead_code)]
pub struct ThermalLbm {
    /// Underlying D2Q9 thermal distribution field.
    pub thermal: ThermalD2Q9,
    /// Optional Boussinesq coupling (buoyancy body force).
    pub boussinesq: Option<BoussinesqCoupling>,
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub ny: usize,
    /// Time step counter.
    pub step_count: usize,
}
#[allow(dead_code)]
impl ThermalLbm {
    /// Create a new `ThermalLbm` with uniform initial temperature `t_init`.
    pub fn new(nx: usize, ny: usize, tau_t: f64, t_init: f64) -> Self {
        let thermal = ThermalD2Q9::new(nx, ny, tau_t, t_init);
        Self {
            thermal,
            boussinesq: None,
            nx,
            ny,
            step_count: 0,
        }
    }
    /// Enable Boussinesq buoyancy coupling.
    pub fn with_boussinesq(mut self, gravity: [f64; 2], beta: f64, t_ref: f64) -> Self {
        self.boussinesq = Some(BoussinesqCoupling::new(gravity, beta, t_ref));
        self
    }
    /// Perform one thermal time step: collide → stream → compute_temperature.
    pub fn step_thermal(&mut self, velocities: &[[f64; 2]]) {
        self.thermal.step(velocities);
        self.step_count += 1;
    }
    /// Apply fixed-temperature Dirichlet BC on bottom wall.
    pub fn apply_bc_bottom(&mut self, t_bot: f64, velocities: &[[f64; 2]]) {
        self.thermal.apply_temperature_bc_bottom(t_bot, velocities);
    }
    /// Apply fixed-temperature Dirichlet BC on top wall.
    pub fn apply_bc_top(&mut self, t_top: f64, velocities: &[[f64; 2]]) {
        self.thermal.apply_temperature_bc_top(t_top, velocities);
    }
    /// Apply fixed-temperature BC on left wall.
    pub fn apply_bc_left(&mut self, t_left: f64, velocities: &[[f64; 2]]) {
        self.thermal.apply_temperature_bc_left(t_left, velocities);
    }
    /// Apply fixed-temperature BC on right wall.
    pub fn apply_bc_right(&mut self, t_right: f64, velocities: &[[f64; 2]]) {
        self.thermal.apply_temperature_bc_right(t_right, velocities);
    }
    /// Get a reference to the temperature field.
    pub fn temperature(&self) -> &[f64] {
        &self.thermal.temperature
    }
    /// Average temperature over the entire domain.
    pub fn average_temperature(&self) -> f64 {
        self.thermal.average_temperature()
    }
    /// Temperature variance over the domain.
    pub fn temperature_variance(&self) -> f64 {
        self.thermal.temperature_variance()
    }
    /// Compute local Nusselt number at the bottom wall.
    ///
    /// Wraps `NusseltComputation::local_nusselt_bottom`.
    pub fn local_nusselt_bottom(&self, t_hot: f64, t_cold: f64) -> Vec<f64> {
        NusseltComputation::local_nusselt_bottom(
            &self.thermal.temperature,
            self.nx,
            self.ny,
            t_hot,
            t_cold,
        )
    }
    /// Compute average Nusselt number at the bottom wall.
    pub fn average_nusselt_bottom(&self, t_hot: f64, t_cold: f64) -> f64 {
        let local = self.local_nusselt_bottom(t_hot, t_cold);
        NusseltComputation::average_nusselt(&local)
    }
    /// Compute buoyancy forces for the entire domain.
    ///
    /// Returns `None` if no Boussinesq coupling is configured.
    pub fn compute_buoyancy_forces(&self, densities: &[f64]) -> Option<Vec<[f64; 2]>> {
        self.boussinesq
            .as_ref()
            .map(|b| b.compute_forces(densities, &self.thermal.temperature))
    }
    /// Apply Guo buoyancy forcing to D2Q9 distribution at node `k`.
    ///
    /// No-op if no Boussinesq coupling is configured.
    pub fn apply_buoyancy_at(&self, f: &mut [f64; 9], k: usize, rho: f64, u: [f64; 2], tau: f64) {
        if let Some(b) = &self.boussinesq {
            let temp = self.thermal.temperature[k];
            b.apply_to_distribution(f, rho, u, temp, tau);
        }
    }
    /// Energy equilibrium distribution at a node.
    ///
    /// Wraps `ThermalD2Q9::equilibrium`.
    pub fn equilibrium(&self, temp: f64, u: [f64; 2], i: usize) -> f64 {
        self.thermal.equilibrium(temp, u, i)
    }
    /// Reset temperature field to a uniform value.
    pub fn reset_temperature(&mut self, t_uniform: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        for k in 0..n {
            self.thermal.temperature[k] = t_uniform;
            for i in 0..9 {
                self.thermal.g[k][i] = self.thermal.equilibrium(t_uniform, [0.0, 0.0], i);
            }
        }
    }
    /// Set temperature at a specific node (and re-initialize g to equilibrium).
    pub fn set_temperature_at(&mut self, x: usize, y: usize, t: f64, u: [f64; 2]) {
        let k = self.thermal.idx(x, y);
        self.thermal.temperature[k] = t;
        for i in 0..9 {
            self.thermal.g[k][i] = self.thermal.equilibrium(t, u, i);
        }
    }
}
/// Boussinesq buoyancy coupling for thermal LBM.
///
/// The buoyancy body force is:
///
/// ```text
/// F = ρ · g_vec · β · (T − T_ref)
/// ```
///
/// where `g_vec` is the gravitational acceleration vector (in LBM units),
/// `β` is the thermal expansion coefficient, and `T_ref` is the reference
/// temperature.
pub struct BoussinesqCoupling {
    /// Gravitational acceleration vector (e.g. `[0.0, -9.81e-4]` in LBM units).
    pub gravity: [f64; 2],
    /// Thermal expansion coefficient β.
    pub thermal_expansion: f64,
    /// Reference temperature T_ref.
    pub t_reference: f64,
}
impl BoussinesqCoupling {
    /// Create a new `BoussinesqCoupling`.
    pub fn new(gravity: [f64; 2], beta: f64, t_ref: f64) -> Self {
        Self {
            gravity,
            thermal_expansion: beta,
            t_reference: t_ref,
        }
    }
    /// Compute the buoyancy body force `[Fx, Fy]` for the given density and
    /// temperature.
    ///
    /// `F = ρ · g_vec · β · (T − T_ref)`
    pub fn buoyancy_force(&self, rho: f64, temperature: f64) -> [f64; 2] {
        let factor = rho * self.thermal_expansion * (temperature - self.t_reference);
        [self.gravity[0] * factor, self.gravity[1] * factor]
    }
    /// Apply buoyancy to a D2Q9 velocity distribution using the Guo forcing
    /// scheme.
    ///
    /// The Guo forcing term for direction `i` is:
    ///
    /// ```text
    /// F_i = w_i * (1 − 1/(2τ)) * [(e_i − u)/cs² + (e_i·u)/cs⁴ · e_i] · F
    /// ```
    ///
    /// This term is added in-place to `f[i]` for all 9 directions.
    pub fn apply_to_distribution(
        &self,
        f: &mut [f64; 9],
        rho: f64,
        u: [f64; 2],
        temp: f64,
        tau: f64,
    ) {
        let force = self.buoyancy_force(rho, temp);
        let cs4 = CS2 * CS2;
        let prefactor = 1.0 - 1.0 / (2.0 * tau);
        for i in 0..9 {
            let w = D2Q9_WEIGHTS[i];
            let c = D2Q9_VELOCITIES[i];
            let ex = c[0] as f64;
            let ey = c[1] as f64;
            let e_dot_u = ex * u[0] + ey * u[1];
            let bracket_x = (ex - u[0]) / CS2 + e_dot_u / cs4 * ex;
            let bracket_y = (ey - u[1]) / CS2 + e_dot_u / cs4 * ey;
            let bracket_dot_f = bracket_x * force[0] + bracket_y * force[1];
            f[i] += w * prefactor * bracket_dot_f;
        }
    }
    /// Compute buoyancy forces for the entire domain.
    ///
    /// Returns a `Vec<[f64; 2]>` of forces, one per cell.
    #[allow(dead_code)]
    pub fn compute_forces(&self, densities: &[f64], temperatures: &[f64]) -> Vec<[f64; 2]> {
        densities
            .iter()
            .zip(temperatures.iter())
            .map(|(&rho, &t)| self.buoyancy_force(rho, t))
            .collect()
    }
}
/// Configuration and helper utilities for the Rayleigh–Bénard convection
/// scenario (heated bottom plate, cooled top plate).
pub struct RayleighBenardSetup {
    /// Number of lattice nodes in x.
    pub nx: usize,
    /// Number of lattice nodes in y.
    pub ny: usize,
    /// Bottom wall (hot) temperature.
    pub t_hot: f64,
    /// Top wall (cold) temperature.
    pub t_cold: f64,
    /// Rayleigh number Ra = g β ΔT H³ / (ν α).
    pub ra: f64,
    /// Prandtl number Pr = ν / α.
    pub pr: f64,
}
impl RayleighBenardSetup {
    /// Create a new `RayleighBenardSetup`.
    pub fn new(nx: usize, ny: usize, t_hot: f64, t_cold: f64, ra: f64, pr: f64) -> Self {
        Self {
            nx,
            ny,
            t_hot,
            t_cold,
            ra,
            pr,
        }
    }
    /// Derive LBM transport parameters from Ra and Pr.
    ///
    /// Using a reference velocity `U_ref = 1e-2` and the channel height
    /// `H = ny` as the characteristic length:
    ///
    /// ```text
    /// ν  = sqrt(Pr / Ra) · U_ref · H
    /// α  = ν / Pr
    /// τ_f = 3ν + 0.5
    /// τ_t = 3α + 0.5
    /// β  = Ra · ν · α / (g_ref · ΔT · H³)   [with g_ref = 1e-4 in LBM units]
    /// ```
    ///
    /// Returns `(tau_f, tau_t, beta)`.
    pub fn lbm_parameters(&self) -> (f64, f64, f64) {
        let h = self.ny as f64;
        let u_ref = 1e-2_f64;
        let g_ref = 1e-4_f64;
        let nu = (self.pr / self.ra).sqrt() * u_ref * h;
        let alpha = nu / self.pr;
        let tau_f = 3.0 * nu + 0.5;
        let tau_t = 3.0 * alpha + 0.5;
        let delta_t = (self.t_hot - self.t_cold).abs().max(1e-30);
        let beta = self.ra * nu * alpha / (g_ref * delta_t * h * h * h);
        (tau_f, tau_t, beta)
    }
    /// Estimate the Nusselt number from the vertical temperature gradient at
    /// the boundaries.
    ///
    /// ```text
    /// Nu = -dT/dy|_wall · H / ΔT
    /// ```
    ///
    /// The gradient at the bottom wall is approximated by a one-sided finite
    /// difference between rows `y = 0` and `y = 1`.
    ///
    /// `temperature` is a 2-D array indexed `temperature[y][x]`.
    pub fn nusselt_number(&self, temperature: &[Vec<f64>]) -> f64 {
        let h = self.ny as f64;
        let delta_t = (self.t_hot - self.t_cold).abs().max(1e-30);
        let mut grad_sum = 0.0_f64;
        for x in 0..self.nx {
            grad_sum += temperature[1][x] - temperature[0][x];
        }
        let avg_grad = grad_sum / self.nx as f64;
        (-avg_grad).abs() * h / delta_t
    }
    /// Initialize a linear temperature profile from hot (bottom) to cold (top).
    #[allow(dead_code)]
    pub fn initial_temperature_profile(&self) -> Vec<Vec<f64>> {
        let mut temperature = vec![vec![0.0; self.nx]; self.ny];
        for y in 0..self.ny {
            let frac = y as f64 / (self.ny - 1).max(1) as f64;
            let t = self.t_hot + (self.t_cold - self.t_hot) * frac;
            for x in 0..self.nx {
                temperature[y][x] = t;
            }
        }
        temperature
    }
}
/// Conjugate heat transfer for thermal LBM.
///
/// Models heat conduction in a solid region that is thermally coupled
/// to the fluid region at the interface.  The solid has its own thermal
/// diffusivity `alpha_s` and a separate temperature field.
#[allow(dead_code)]
pub struct ConjugateHeatTransfer {
    /// Temperature field in the solid domain.
    pub solid_temperature: Vec<f64>,
    /// Number of nodes in x.
    pub nx: usize,
    /// Number of nodes in y.
    pub ny: usize,
    /// Solid thermal diffusivity.
    pub alpha_s: f64,
    /// Boolean mask: `true` if the cell is solid.
    pub is_solid: Vec<bool>,
}
#[allow(dead_code)]
impl ConjugateHeatTransfer {
    /// Create a new conjugate heat transfer field.
    ///
    /// All solid cells are initialized to `t_init`.
    pub fn new(nx: usize, ny: usize, alpha_s: f64, t_init: f64, is_solid: Vec<bool>) -> Self {
        let n = nx * ny;
        Self {
            solid_temperature: vec![t_init; n],
            nx,
            ny,
            alpha_s,
            is_solid,
        }
    }
    /// Flat index.
    #[inline]
    fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }
    /// Perform one diffusion step in the solid using a simple
    /// explicit finite-difference scheme:
    ///
    /// ```text
    /// T_new = T_old + α_s * (T_{x+1} + T_{x-1} + T_{y+1} + T_{y-1} - 4*T)
    /// ```
    ///
    /// At the solid-fluid interface, the fluid temperature is used as
    /// a Dirichlet boundary condition.
    pub fn diffuse_solid(&mut self, fluid_temperature: &[f64]) {
        let nx = self.nx;
        let ny = self.ny;
        let mut t_new = self.solid_temperature.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                if !self.is_solid[k] {
                    continue;
                }
                let t_center = self.solid_temperature[k];
                let get_neighbor = |nx2: usize, ny2: usize| -> f64 {
                    let kn = ny2 * nx + nx2;
                    if self.is_solid[kn] {
                        self.solid_temperature[kn]
                    } else {
                        fluid_temperature[kn]
                    }
                };
                let xp = if x + 1 < nx { x + 1 } else { 0 };
                let xm = if x > 0 { x - 1 } else { nx - 1 };
                let yp = if y + 1 < ny { y + 1 } else { 0 };
                let ym = if y > 0 { y - 1 } else { ny - 1 };
                let t_xp = get_neighbor(xp, y);
                let t_xm = get_neighbor(xm, y);
                let t_yp = get_neighbor(x, yp);
                let t_ym = get_neighbor(x, ym);
                let laplacian = t_xp + t_xm + t_yp + t_ym - 4.0 * t_center;
                t_new[k] = t_center + self.alpha_s * laplacian;
            }
        }
        self.solid_temperature = t_new;
    }
    /// Update fluid thermal BC at the interface: set the fluid thermal
    /// distribution to equilibrium at the solid temperature for interface cells.
    pub fn couple_to_fluid(&self, thermal: &mut ThermalD2Q9, velocities: &[[f64; 2]]) {
        let nx = self.nx;
        let ny = self.ny;
        for y in 0..ny {
            for x in 0..nx {
                let k = self.idx(x, y);
                if self.is_solid[k] {
                    continue;
                }
                let xp = if x + 1 < nx { x + 1 } else { 0 };
                let xm = if x > 0 { x - 1 } else { nx - 1 };
                let yp = if y + 1 < ny { y + 1 } else { 0 };
                let ym = if y > 0 { y - 1 } else { ny - 1 };
                let has_solid_neighbor = self.is_solid[yp * nx + x]
                    || self.is_solid[ym * nx + x]
                    || self.is_solid[y * nx + xp]
                    || self.is_solid[y * nx + xm];
                if has_solid_neighbor {
                    let mut t_sum = 0.0;
                    let mut count = 0.0;
                    for &nk in &[yp * nx + x, ym * nx + x, y * nx + xp, y * nx + xm] {
                        if self.is_solid[nk] {
                            t_sum += self.solid_temperature[nk];
                            count += 1.0;
                        }
                    }
                    if count > 0.0 {
                        let t_interface = t_sum / count;
                        let t_mixed = 0.5 * (thermal.temperature[k] + t_interface);
                        let u = velocities[k];
                        for i in 0..9 {
                            thermal.g[k][i] = thermal.equilibrium(t_mixed, u, i);
                        }
                        thermal.temperature[k] = t_mixed;
                    }
                }
            }
        }
    }
}
