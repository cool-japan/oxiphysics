// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compressible flow LBM.
//!
//! Provides a Mach-corrected double-distribution compressible LBM, a
//! nearly-inviscid Euler LBM, shock capturing, compressible Navier-Stokes,
//! a Rayleigh-Taylor two-population LBM and a 1-D reactive-Euler detonation
//! solver.
//!
//! All D2Q9 solvers in this module run at the lattice sound speed
//! `cs² = 1/3` ([`LATTICE_CS2`]); the D2Q9 weights only reproduce the
//! Hermite moments for that value.  Thermodynamics with an arbitrary
//! adiabatic exponent γ is carried by a separate total-energy population and
//! a pressure-correction force (see [`energy_equilibrium`]).

mod ddf;
#[cfg(test)]
mod ddf_tests;
mod detonation;
mod rayleigh_taylor;

pub use ddf::energy_equilibrium;
pub use detonation::DetonationLbm;
pub use rayleigh_taylor::{RayleighTaylorConfig, RayleighTaylorLbm};

use ddf::{DdfParams, DdfState};

/// Lattice (isothermal) speed of sound squared of the D2Q9 lattice, `1/3`.
pub const LATTICE_CS2: f64 = 1.0 / 3.0;

/// D2Q9 weights.
const W9: [f64; 9] = [
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
/// D2Q9 velocity ex.
const EX9: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 velocity ey.
const EY9: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

/// Third-order Hermite D2Q9 equilibrium.
///
/// `f_eq_i = w_i ρ [1 + c·u/cs² + ((c·u)² − cs²|u|²)/(2cs⁴)
///                + ((c·u)³ − 3cs²(c·u)|u|²)/(6cs⁶)]`.
///
/// The D2Q9 weights are only consistent with `cs2 = 1/3` ([`LATTICE_CS2`]);
/// for that value `Σ f_eq = ρ`, `Σ f_eq c = ρu`, `Σ f_eq cc = ρ(cs² I + uu)`
/// and the lattice-supported third moments `Σ f_eq c_x c_y² = ρ(cs² u_x +
/// u_x u_y²)` (and the x↔y mirror) hold exactly.  Other `cs2` values do not
/// conserve mass or momentum on D2Q9.
pub fn equilibrium_compressible(rho: f64, ux: f64, uy: f64, cs2: f64) -> [f64; 9] {
    let mut feq = [0.0f64; 9];
    let usq = ux * ux + uy * uy;
    let inv_cs2 = 1.0 / cs2;
    let inv_cs4 = inv_cs2 * inv_cs2;
    let inv_cs6 = inv_cs4 * inv_cs2;
    for q in 0..9 {
        let eu = f64::from(EX9[q]) * ux + f64::from(EY9[q]) * uy;
        // Third-order Hermite expansion (Ma³).
        feq[q] = W9[q]
            * rho
            * (1.0
                + eu * inv_cs2
                + 0.5 * (eu * eu * inv_cs4 - usq * inv_cs2)
                + eu * eu * eu * inv_cs6 / 6.0
                - eu * usq * inv_cs4 / 2.0);
    }
    feq
}

/// Compute local Mach number.
///
/// Returns Ma = |u| / cs.
pub fn local_mach(ux: f64, uy: f64, cs: f64) -> f64 {
    let u_mag = (ux * ux + uy * uy).sqrt();
    u_mag / cs
}

/// Local Mach number with respect to the adiabatic sound speed `sqrt(γT)`.
///
/// Returns 0 for a fluid at rest and `+∞` for a moving fluid at `T ≤ 0`.
fn adiabatic_mach(u: [f64; 2], gamma: f64, temp: f64) -> f64 {
    let u_mag = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let sound = (gamma * temp.max(0.0)).sqrt();
    if sound > 0.0 {
        u_mag / sound
    } else if u_mag > 0.0 {
        f64::INFINITY
    } else {
        0.0
    }
}

/// Compute entropy generation rate.
///
/// Returns sigma = rho * T * Ds/Dt estimate.
pub fn entropy_generation(rho: f64, temp: f64, du_dx: f64, dv_dy: f64, mu: f64) -> f64 {
    let phi_visc = 2.0 * mu * (du_dx * du_dx + dv_dy * dv_dy);
    phi_visc / (rho * temp + 1e-30)
}

/// High-Mach compressible LBM with a Mach-corrected equilibrium and an
/// energy double distribution.
///
/// * `f` (density/momentum) relaxes towards the D2Q9 Hermite equilibrium of
///   order `ma_order` (2 = standard, ≥3 adds the lattice-supported cubic
///   terms, see [`MaCorrection`]) at the lattice sound speed `cs² = 1/3`.
/// * `g` evolves the total energy `ρE = ρ(|u|²/2 + T/(γ−1))` with the
///   equilibrium [`energy_equilibrium`].
/// * The ideal-gas pressure `p = ρT` (`R = 1`, lattice units) is imposed via a
///   Guo pressure-correction force `F = −∇(ρT − ρcs²)`, so γ controls the
///   thermodynamics while the lattice stays at `cs² = 1/3`.
///
/// Both populations relax with `omega` (unit Prandtl number); use
/// [`NsCompressibleLbm`] for an independent thermal relaxation rate.  The
/// shear viscosity is `ν = cs²(1/ω − 1/2)`.  The domain is periodic.
#[derive(Debug, Clone)]
pub struct CompressibleLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Relaxation parameter omega (both populations).
    pub omega: f64,
    /// Density/momentum distribution functions f\[y\]\[x\]\[q\].
    pub f: Vec<Vec<Vec<f64>>>,
    /// Total-energy distribution functions g\[y\]\[x\]\[q\].
    pub g: Vec<Vec<Vec<f64>>>,
    /// Density field rho\[y\]\[x\].
    pub rho: Vec<Vec<f64>>,
    /// Velocity field u\[y\]\[x\]\[2\] (Guo-corrected).
    pub u: Vec<Vec<[f64; 2]>>,
    /// Temperature field T\[y\]\[x\] (`p = ρT`, lattice units).
    pub temp: Vec<Vec<f64>>,
    /// Pressure-correction force field F\[y\]\[x\] used in the last collision.
    pub force: Vec<Vec<[f64; 2]>>,
    /// Lattice (isothermal) speed of sound `sqrt(1/3)`.
    pub cs: f64,
    /// Adiabatic exponent gamma.
    pub gamma: f64,
    /// Hermite order of the `f` equilibrium (2 or 3).
    pub ma_order: usize,
}

impl CompressibleLbm {
    /// Create a new compressible LBM simulation at rest with `ρ = 1` and the
    /// reference temperature `T = cs² = 1/3`.
    pub fn new(nx: usize, ny: usize, omega: f64, gamma: f64, ma_order: usize) -> Self {
        let mut sim = Self {
            nx,
            ny,
            omega,
            f: vec![vec![vec![0.0f64; 9]; nx]; ny],
            g: vec![vec![vec![0.0f64; 9]; nx]; ny],
            rho: vec![vec![1.0f64; nx]; ny],
            u: vec![vec![[0.0f64; 2]; nx]; ny],
            temp: vec![vec![LATTICE_CS2; nx]; ny],
            force: vec![vec![[0.0f64; 2]; nx]; ny],
            cs: LATTICE_CS2.sqrt(),
            gamma,
            ma_order,
        };
        sim.initialize();
        sim
    }

    fn params(&self) -> DdfParams {
        DdfParams {
            omega_f: self.omega,
            omega_g: self.omega,
            gamma: self.gamma,
            ma_order: self.ma_order,
        }
    }

    fn state(&mut self) -> DdfState<'_> {
        DdfState {
            nx: self.nx,
            ny: self.ny,
            f: &mut self.f,
            g: &mut self.g,
            rho: &mut self.rho,
            u: &mut self.u,
            temp: &mut self.temp,
            force: &mut self.force,
        }
    }

    /// Initialize `f` and `g` at equilibrium from the current `rho`, `u` and
    /// `temp` fields.
    pub fn initialize(&mut self) {
        let params = self.params();
        ddf::initialize(&mut self.state(), params);
    }

    /// BGK collision of both populations (Guo-forced for `f`).
    pub fn collide(&mut self) {
        let params = self.params();
        ddf::collide(&mut self.state(), params);
    }

    /// Periodic streaming step of both populations.
    pub fn stream(&mut self) {
        ddf::stream_periodic(&mut self.f, self.nx, self.ny);
        ddf::stream_periodic(&mut self.g, self.nx, self.ny);
    }

    /// Update `ρ`, `u`, `T` and the pressure-correction force.
    pub fn update_macros(&mut self) {
        let gamma = self.gamma;
        ddf::update_macros(&mut self.state(), gamma);
    }

    /// Full time step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }

    /// Pressure `p = ρT` at (x, y).
    pub fn pressure_at(&self, x: usize, y: usize) -> f64 {
        self.rho[y][x] * self.temp[y][x]
    }

    /// Local Mach number at (x, y) with respect to the adiabatic sound speed
    /// `sqrt(γT)`.
    pub fn mach_at(&self, x: usize, y: usize) -> f64 {
        adiabatic_mach(self.u[y][x], self.gamma, self.temp[y][x])
    }

    /// Maximum Mach number in the domain.
    pub fn max_mach(&self) -> f64 {
        let mut ma_max = 0.0f64;
        for y in 0..self.ny {
            for x in 0..self.nx {
                ma_max = ma_max.max(self.mach_at(x, y));
            }
        }
        ma_max
    }

    /// Conserved totals `(mass, momentum, total energy)` over the domain.
    pub fn conserved_totals(&self) -> (f64, [f64; 2], f64) {
        ddf::totals(&self.f, &self.g, &self.rho, &self.u)
    }
}

/// Nearly-inviscid compressible Euler LBM (double distribution).
///
/// Same double-distribution model as [`CompressibleLbm`] with a third-order
/// equilibrium and a relaxation frequency close to 2: the shear viscosity is
/// `ν = cs²(1/ω − 1/2)` (≈ 8.4·10⁻⁴ at the default `ω = 1.99`).  `ω = 2`
/// exactly is the formally inviscid but neutrally unstable limit, so
/// `omega` must lie in `(0, 2)`.  The energy population `g` is collided and
/// streamed alongside `f`, and the pressure `p = ρT` enters through the
/// pressure-correction force.  The domain is periodic.
#[derive(Debug, Clone)]
pub struct EulerLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Density/momentum distribution functions.
    pub f: Vec<Vec<Vec<f64>>>,
    /// Total-energy distribution functions g.
    pub g: Vec<Vec<Vec<f64>>>,
    /// Density.
    pub rho: Vec<Vec<f64>>,
    /// Velocity (Guo-corrected).
    pub u: Vec<Vec<[f64; 2]>>,
    /// Total energy density `ρE = Σ g`.
    pub energy: Vec<Vec<f64>>,
    /// Temperature (`p = ρT`).
    pub temp: Vec<Vec<f64>>,
    /// Pressure-correction force used in the last collision.
    pub force: Vec<Vec<[f64; 2]>>,
    /// Lattice (isothermal) speed of sound `sqrt(1/3)`.
    pub cs: f64,
    /// Relaxation frequency (both populations), in `(0, 2)`.
    pub omega: f64,
    /// Adiabatic exponent γ.
    pub gamma: f64,
}

impl EulerLbm {
    /// Default near-inviscid relaxation frequency.
    pub const DEFAULT_OMEGA: f64 = 1.99;

    /// Hermite order of the `f` equilibrium used by the Euler solver.
    const MA_ORDER: usize = 3;

    /// Create a new Euler LBM at rest with `ρ = 1`, `T = 1/3`.
    pub fn new(nx: usize, ny: usize, omega: f64, gamma: f64) -> Self {
        let energy0 = LATTICE_CS2 / (gamma - 1.0);
        let mut sim = Self {
            nx,
            ny,
            f: vec![vec![vec![0.0f64; 9]; nx]; ny],
            g: vec![vec![vec![0.0f64; 9]; nx]; ny],
            rho: vec![vec![1.0f64; nx]; ny],
            u: vec![vec![[0.0f64; 2]; nx]; ny],
            energy: vec![vec![energy0; nx]; ny],
            temp: vec![vec![LATTICE_CS2; nx]; ny],
            force: vec![vec![[0.0f64; 2]; nx]; ny],
            cs: LATTICE_CS2.sqrt(),
            omega,
            gamma,
        };
        sim.initialize();
        sim
    }

    /// Create a new Euler LBM with [`Self::DEFAULT_OMEGA`].
    pub fn with_default_omega(nx: usize, ny: usize, gamma: f64) -> Self {
        Self::new(nx, ny, Self::DEFAULT_OMEGA, gamma)
    }

    /// Kinematic shear viscosity `ν = cs²(1/ω − 1/2)` of the scheme.
    pub fn viscosity(&self) -> f64 {
        LATTICE_CS2 * (1.0 / self.omega - 0.5)
    }

    fn params(&self) -> DdfParams {
        DdfParams {
            omega_f: self.omega,
            omega_g: self.omega,
            gamma: self.gamma,
            ma_order: Self::MA_ORDER,
        }
    }

    fn state(&mut self) -> DdfState<'_> {
        DdfState {
            nx: self.nx,
            ny: self.ny,
            f: &mut self.f,
            g: &mut self.g,
            rho: &mut self.rho,
            u: &mut self.u,
            temp: &mut self.temp,
            force: &mut self.force,
        }
    }

    fn refresh_energy(&mut self) {
        self.energy = self
            .g
            .iter()
            .map(|row| row.iter().map(|cell| cell.iter().sum()).collect())
            .collect();
    }

    /// Initialize `f` and `g` at equilibrium from the current `rho`, `u` and
    /// `temp` fields.
    pub fn initialize(&mut self) {
        let params = self.params();
        ddf::initialize(&mut self.state(), params);
        self.refresh_energy();
    }

    /// BGK collision `f ← f − ω(f − f_eq) + S_F`, `g ← g − ω(g − g_eq)`.
    pub fn collide(&mut self) {
        let params = self.params();
        ddf::collide(&mut self.state(), params);
    }

    /// Periodic streaming of both populations.
    pub fn stream(&mut self) {
        ddf::stream_periodic(&mut self.f, self.nx, self.ny);
        ddf::stream_periodic(&mut self.g, self.nx, self.ny);
    }

    /// Update `ρ`, `u`, `T`, `ρE` and the pressure-correction force.
    pub fn update_macros(&mut self) {
        let gamma = self.gamma;
        ddf::update_macros(&mut self.state(), gamma);
        self.refresh_energy();
    }

    /// Full Euler step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }

    /// Conserved totals `(mass, momentum, total energy)` over the domain.
    pub fn conserved_totals(&self) -> (f64, [f64; 2], f64) {
        ddf::totals(&self.f, &self.g, &self.rho, &self.u)
    }

    /// Compute Roe flux at interface between left and right states.
    pub fn roe_flux_x(
        &self,
        rho_l: f64,
        u_l: [f64; 2],
        p_l: f64,
        rho_r: f64,
        u_r: [f64; 2],
        p_r: f64,
    ) -> [f64; 4] {
        euler_flux_roe(rho_l, u_l, p_l, rho_r, u_r, p_r)
    }
}

/// Shock-capturing LBM with artificial viscosity.
///
/// Adds bulk viscosity term and artificial dissipation near shocks.
#[derive(Debug, Clone)]
pub struct ShockCaptureLbm {
    /// Underlying compressible LBM.
    pub lbm: CompressibleLbm,
    /// Artificial viscosity coefficient.
    pub nu_art: f64,
    /// Shock sensor threshold.
    pub shock_threshold: f64,
    /// Shock indicator field.
    pub shock_flag: Vec<Vec<bool>>,
}

impl ShockCaptureLbm {
    /// Create a new shock-capturing LBM.
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        gamma: f64,
        nu_art: f64,
        shock_threshold: f64,
    ) -> Self {
        let lbm = CompressibleLbm::new(nx, ny, omega, gamma, 2);
        let shock_flag = vec![vec![false; nx]; ny];
        Self {
            lbm,
            nu_art,
            shock_threshold,
            shock_flag,
        }
    }

    /// Detect shocks using pressure gradient sensor.
    pub fn detect_shocks(&mut self) {
        let nx = self.lbm.nx;
        let ny = self.lbm.ny;
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let rho_c = self.lbm.rho[y][x];
                let rho_r = self.lbm.rho[y][x + 1];
                let rho_l = self.lbm.rho[y][x - 1];
                let sensor = ((rho_r - 2.0 * rho_c + rho_l).abs()) / (rho_c + 1e-10);
                self.shock_flag[y][x] = sensor > self.shock_threshold;
            }
        }
    }

    /// Apply artificial viscosity at shock cells.
    pub fn apply_artificial_viscosity(&mut self) {
        let nx = self.lbm.nx;
        let ny = self.lbm.ny;
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                if self.shock_flag[y][x] {
                    // Add diffusion: smear distribution functions
                    for q in 0..9 {
                        let lap = self.lbm.f[y][x + 1][q] - 2.0 * self.lbm.f[y][x][q]
                            + self.lbm.f[y][x - 1][q];
                        self.lbm.f[y][x][q] += self.nu_art * lap;
                    }
                }
            }
        }
    }

    /// Full time step with shock capturing.
    pub fn step(&mut self) {
        self.lbm.step();
        self.detect_shocks();
        self.apply_artificial_viscosity();
    }

    /// Count number of shock cells.
    pub fn shock_cell_count(&self) -> usize {
        self.shock_flag
            .iter()
            .flat_map(|row| row.iter())
            .filter(|&&b| b)
            .count()
    }
}

/// Compressible Navier-Stokes via double distribution (density + energy).
///
/// Uses the same double-distribution model as [`CompressibleLbm`] (third
/// order `f` equilibrium at `cs² = 1/3`, total-energy population `g`,
/// pressure-correction force `F = −∇(ρT − ρcs²)`) with independent
/// relaxation frequencies for momentum (`omega_f`) and energy (`omega_g`).
/// The temperature is recovered as `T = (γ−1)(Σg/ρ − |u|²/2)`.  The domain is
/// periodic.
///
/// Stability note: the adiabatic sound speed `sqrt(γT)` must stay well below
/// one lattice cell per step, so temperatures should remain of the order of
/// the reference value `T = cs² = 1/3`.
#[derive(Debug, Clone)]
pub struct NsCompressibleLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Density/momentum distribution f.
    pub f: Vec<Vec<Vec<f64>>>,
    /// Total-energy distribution g.
    pub g: Vec<Vec<Vec<f64>>>,
    /// Relaxation for f.
    pub omega_f: f64,
    /// Relaxation for g.
    pub omega_g: f64,
    /// Density field.
    pub rho: Vec<Vec<f64>>,
    /// Velocity field (Guo-corrected).
    pub u: Vec<Vec<[f64; 2]>>,
    /// Temperature field (`p = ρT`).
    pub temp: Vec<Vec<f64>>,
    /// Pressure-correction force used in the last collision.
    pub force: Vec<Vec<[f64; 2]>>,
    /// Nominal Prandtl number supplied by the caller.  It is stored for
    /// reference only: the transport coefficients of the scheme are set by
    /// `omega_f` (viscosity) and `omega_g` (energy diffusion).
    pub pr: f64,
    /// Gamma.
    pub gamma: f64,
}

impl NsCompressibleLbm {
    /// Hermite order of the `f` equilibrium.
    const MA_ORDER: usize = 3;

    /// Create a new NS compressible LBM at rest with `ρ = 1` and the
    /// reference temperature `T = cs² = 1/3`.
    pub fn new(nx: usize, ny: usize, omega_f: f64, omega_g: f64, pr: f64, gamma: f64) -> Self {
        let mut sim = Self {
            nx,
            ny,
            f: vec![vec![vec![0.0f64; 9]; nx]; ny],
            g: vec![vec![vec![0.0f64; 9]; nx]; ny],
            omega_f,
            omega_g,
            rho: vec![vec![1.0f64; nx]; ny],
            u: vec![vec![[0.0f64; 2]; nx]; ny],
            temp: vec![vec![LATTICE_CS2; nx]; ny],
            force: vec![vec![[0.0f64; 2]; nx]; ny],
            pr,
            gamma,
        };
        sim.initialize();
        sim
    }

    fn params(&self) -> DdfParams {
        DdfParams {
            omega_f: self.omega_f,
            omega_g: self.omega_g,
            gamma: self.gamma,
            ma_order: Self::MA_ORDER,
        }
    }

    fn state(&mut self) -> DdfState<'_> {
        DdfState {
            nx: self.nx,
            ny: self.ny,
            f: &mut self.f,
            g: &mut self.g,
            rho: &mut self.rho,
            u: &mut self.u,
            temp: &mut self.temp,
            force: &mut self.force,
        }
    }

    /// Initialize `f` and `g` at equilibrium from the current `rho`, `u` and
    /// `temp` fields.
    pub fn initialize(&mut self) {
        let params = self.params();
        ddf::initialize(&mut self.state(), params);
    }

    /// Collision step for both f (Guo-forced) and g.
    pub fn collide(&mut self) {
        let params = self.params();
        ddf::collide(&mut self.state(), params);
    }

    /// Periodic streaming step.
    pub fn stream(&mut self) {
        ddf::stream_periodic(&mut self.f, self.nx, self.ny);
        ddf::stream_periodic(&mut self.g, self.nx, self.ny);
    }

    /// Update `ρ`, `u`, `T` and the pressure-correction force.
    pub fn update_macros(&mut self) {
        let gamma = self.gamma;
        ddf::update_macros(&mut self.state(), gamma);
    }

    /// Full NS compressible step.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
        self.update_macros();
    }

    /// Conserved totals `(mass, momentum, total energy)` over the domain.
    pub fn conserved_totals(&self) -> (f64, [f64; 2], f64) {
        ddf::totals(&self.f, &self.g, &self.rho, &self.u)
    }
}

/// Mach number correction for extended equilibrium.
///
/// `order = 2` gives the standard second-order Hermite D2Q9 equilibrium;
/// `order >= 3` adds the third-order Hermite terms
/// `((c·u)³ − 3cs²(c·u)|u|²)/(6cs⁶)`, which on D2Q9 at `cs² = 1/3` restore
/// the lattice-supported third moments `Σ f_eq c_x c_y² = ρ(cs² u_x +
/// u_x u_y²)` without changing mass, momentum or the pressure tensor.
/// `cs` must be the lattice sound speed `sqrt(1/3)` for the D2Q9 weights to be
/// consistent.
#[derive(Debug, Clone)]
pub struct MaCorrection {
    /// Maximum Mach number order (2 or 3).
    pub order: usize,
    /// Reference speed of sound.
    pub cs: f64,
}

impl MaCorrection {
    /// Create a new Mach correction calculator.
    pub fn new(order: usize, cs: f64) -> Self {
        Self { order, cs }
    }

    /// Compute corrected equilibrium at given Ma.
    pub fn corrected_feq(&self, rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let cs2 = self.cs * self.cs;
        let mut feq = [0.0f64; 9];
        let usq = ux * ux + uy * uy;
        for q in 0..9 {
            let eu = f64::from(EX9[q]) * ux + f64::from(EY9[q]) * uy;
            let base = 1.0 + eu / cs2 + 0.5 * (eu * eu / (cs2 * cs2) - usq / cs2);
            let correction = if self.order >= 3 {
                eu * eu * eu / (6.0 * cs2 * cs2 * cs2) - eu * usq / (2.0 * cs2 * cs2)
            } else {
                0.0
            };
            feq[q] = W9[q] * rho * (base + correction);
        }
        feq
    }

    /// Estimate Mach correction error at given Ma.
    pub fn truncation_error(&self, ma: f64) -> f64 {
        ma.powi(self.order as i32 + 1)
    }
}

/// 1D Riemann problem (Sod shock tube test).
///
/// Implements contact discontinuity and shock wave evolution.
#[derive(Debug, Clone)]
pub struct ShockWave1D {
    /// Number of cells.
    pub n: usize,
    /// Density profile.
    pub rho: Vec<f64>,
    /// Velocity profile.
    pub vel: Vec<f64>,
    /// Pressure profile.
    pub pressure: Vec<f64>,
    /// Energy profile.
    pub energy: Vec<f64>,
    /// Grid spacing dx.
    pub dx: f64,
    /// Time step dt.
    pub dt: f64,
    /// Adiabatic index gamma.
    pub gamma: f64,
    /// Current time.
    pub time: f64,
}

impl ShockWave1D {
    /// Create a new Sod shock tube problem.
    pub fn new_sod(n: usize, dx: f64, dt: f64, gamma: f64) -> Self {
        let mut rho = vec![0.0f64; n];
        let mut pressure = vec![0.0f64; n];
        let vel = vec![0.0f64; n];
        let mid = n / 2;
        for i in 0..n {
            if i < mid {
                rho[i] = 1.0;
                pressure[i] = 1.0;
            } else {
                rho[i] = 0.125;
                pressure[i] = 0.1;
            }
        }
        let energy: Vec<f64> = (0..n)
            .map(|i| pressure[i] / (gamma - 1.0) + 0.5 * rho[i] * vel[i] * vel[i])
            .collect();
        Self {
            n,
            rho,
            vel,
            pressure,
            energy,
            dx,
            dt,
            gamma,
            time: 0.0,
        }
    }

    /// Compute local speed of sound.
    pub fn sound_speed(&self, i: usize) -> f64 {
        let p = self.pressure[i].max(0.0);
        let rho = self.rho[i].max(1e-10);
        (self.gamma * p / rho).sqrt()
    }

    /// Compute CFL number.
    pub fn cfl(&self) -> f64 {
        let max_speed: f64 = (0..self.n)
            .map(|i| self.vel[i].abs() + self.sound_speed(i))
            .fold(0.0f64, f64::max);
        max_speed * self.dt / self.dx
    }

    /// First-order Godunov step (upwind).
    pub fn step_godunov(&mut self) {
        let n = self.n;
        let dt = self.dt;
        let dx = self.dx;
        let g = self.gamma;
        let mut rho_new = self.rho.clone();
        let mut vel_new = self.vel.clone();
        let mut e_new = self.energy.clone();

        for i in 1..n - 1 {
            // Left and right states
            let rho_l = self.rho[i - 1];
            let rho_r = self.rho[i];
            let u_l = self.vel[i - 1];
            let u_r = self.vel[i];
            let p_l = (g - 1.0) * (self.energy[i - 1] - 0.5 * rho_l * u_l * u_l);
            let p_r = (g - 1.0) * (self.energy[i] - 0.5 * rho_r * u_r * u_r);

            // Simple upwind flux
            let flux_rho_l = if u_l > 0.0 { rho_l * u_l } else { rho_r * u_r };
            let flux_mom_l = if u_l > 0.0 {
                rho_l * u_l * u_l + p_l
            } else {
                rho_r * u_r * u_r + p_r
            };
            let flux_e_l = if u_l > 0.0 {
                (self.energy[i - 1] + p_l) * u_l
            } else {
                (self.energy[i] + p_r) * u_r
            };

            let rho_rr = self.rho[(i + 1).min(n - 1)];
            let u_rr = self.vel[(i + 1).min(n - 1)];
            let p_rr = (g - 1.0) * (self.energy[(i + 1).min(n - 1)] - 0.5 * rho_rr * u_rr * u_rr);
            let flux_rho_r = if u_r > 0.0 {
                rho_r * u_r
            } else {
                rho_rr * u_rr
            };
            let flux_mom_r = if u_r > 0.0 {
                rho_r * u_r * u_r + p_r
            } else {
                rho_rr * u_rr * u_rr + p_rr
            };
            let flux_e_r = if u_r > 0.0 {
                (self.energy[i] + p_r) * u_r
            } else {
                (self.energy[(i + 1).min(n - 1)] + p_rr) * u_rr
            };

            rho_new[i] = self.rho[i] - dt / dx * (flux_rho_r - flux_rho_l);
            let mom = self.rho[i] * self.vel[i] - dt / dx * (flux_mom_r - flux_mom_l);
            e_new[i] = self.energy[i] - dt / dx * (flux_e_r - flux_e_l);
            rho_new[i] = rho_new[i].max(1e-10);
            vel_new[i] = mom / rho_new[i];
            e_new[i] = e_new[i].max(1e-10);
        }
        self.rho = rho_new;
        self.vel = vel_new;
        self.energy = e_new;
        self.pressure = (0..n)
            .map(|i| {
                let p =
                    (g - 1.0) * (self.energy[i] - 0.5 * self.rho[i] * self.vel[i] * self.vel[i]);
                p.max(0.0)
            })
            .collect();
        self.time += dt;
    }

    /// Run for n steps.
    pub fn run(&mut self, nsteps: usize) {
        for _ in 0..nsteps {
            self.step_godunov();
        }
    }
}

/// Supersonic inlet boundary condition.
///
/// Specifies inflow Mach number, angle of attack, density, and pressure.
#[derive(Debug, Clone)]
pub struct SupersonicInlet {
    /// Inflow Mach number.
    pub mach: f64,
    /// Angle of attack in radians.
    pub angle: f64,
    /// Inflow density.
    pub rho_in: f64,
    /// Inflow pressure.
    pub p_in: f64,
    /// Adiabatic index.
    pub gamma: f64,
    /// Speed of sound.
    pub cs: f64,
}

impl SupersonicInlet {
    /// Create a new supersonic inlet.
    pub fn new(mach: f64, angle: f64, rho_in: f64, p_in: f64, gamma: f64) -> Self {
        let cs = (gamma * p_in / rho_in).sqrt();
        Self {
            mach,
            angle,
            rho_in,
            p_in,
            gamma,
            cs,
        }
    }

    /// Inflow velocity vector \[ux, uy\].
    pub fn velocity(&self) -> [f64; 2] {
        let u_mag = self.mach * self.cs;
        [u_mag * self.angle.cos(), u_mag * self.angle.sin()]
    }

    /// Apply inlet BC to the D2Q9 density/momentum populations at column x=0.
    ///
    /// The populations are set to the third-order equilibrium at the lattice
    /// sound speed `cs² = 1/3` (the only value consistent with the D2Q9
    /// weights), so they carry exactly `ρ_in` and `ρ_in u_in`.  The inflow
    /// pressure/temperature is imposed through the energy population of the
    /// double-distribution solvers, not through `f`.
    pub fn apply(&self, f: &mut [Vec<Vec<f64>>], ny: usize) {
        let [ux, uy] = self.velocity();
        let feq = equilibrium_compressible(self.rho_in, ux, uy, LATTICE_CS2);
        for row in f.iter_mut().take(ny) {
            if let Some(cell) = row.first_mut() {
                *cell = feq.to_vec();
            }
        }
    }

    /// Total pressure from isentropic relations.
    pub fn total_pressure(&self) -> f64 {
        self.p_in
            * (1.0 + (self.gamma - 1.0) / 2.0 * self.mach * self.mach)
                .powf(self.gamma / (self.gamma - 1.0))
    }

    /// Total temperature ratio T0/T.
    pub fn total_temperature_ratio(&self) -> f64 {
        1.0 + (self.gamma - 1.0) / 2.0 * self.mach * self.mach
    }
}

/// Compressible flow statistics.
///
/// Computes Mach number, pressure ratio, density ratio, and entropy production.
#[derive(Debug, Clone)]
pub struct CompressibleStats {
    /// Density field snapshot.
    pub rho: Vec<Vec<f64>>,
    /// Velocity field snapshot.
    pub u: Vec<Vec<[f64; 2]>>,
    /// Pressure field snapshot.
    pub p: Vec<Vec<f64>>,
    /// Speed of sound field.
    pub cs_field: Vec<Vec<f64>>,
    /// Gamma.
    pub gamma: f64,
    /// Reference density.
    pub rho_ref: f64,
    /// Reference pressure.
    pub p_ref: f64,
}

impl CompressibleStats {
    /// Create new stats from snapshot.
    pub fn new(rho: Vec<Vec<f64>>, u: Vec<Vec<[f64; 2]>>, p: Vec<Vec<f64>>, gamma: f64) -> Self {
        let ny = rho.len();
        let nx = if ny > 0 { rho[0].len() } else { 0 };
        let rho_ref = 1.0;
        let p_ref = 1.0;
        let cs_field = (0..ny)
            .map(|y| {
                (0..nx)
                    .map(|x| (gamma * p[y][x] / rho[y][x].max(1e-10)).sqrt())
                    .collect()
            })
            .collect();
        Self {
            rho,
            u,
            p,
            cs_field,
            gamma,
            rho_ref,
            p_ref,
        }
    }

    /// Local Mach number at (x, y).
    pub fn mach_at(&self, x: usize, y: usize) -> f64 {
        let [ux, uy] = self.u[y][x];
        let u_mag = (ux * ux + uy * uy).sqrt();
        let cs = self.cs_field[y][x];
        u_mag / cs.max(1e-30)
    }

    /// Max Mach number in domain.
    pub fn max_mach(&self) -> f64 {
        let ny = self.rho.len();
        let nx = if ny > 0 { self.rho[0].len() } else { 0 };
        let mut ma_max = 0.0f64;
        for y in 0..ny {
            for x in 0..nx {
                ma_max = ma_max.max(self.mach_at(x, y));
            }
        }
        ma_max
    }

    /// Pressure ratio p/p_ref.
    pub fn pressure_ratio(&self, x: usize, y: usize) -> f64 {
        self.p[y][x] / self.p_ref
    }

    /// Density ratio rho/rho_ref.
    pub fn density_ratio(&self, x: usize, y: usize) -> f64 {
        self.rho[y][x] / self.rho_ref
    }

    /// Average entropy production (simplified).
    pub fn avg_entropy_production(&self) -> f64 {
        let ny = self.rho.len();
        let nx = if ny > 0 { self.rho[0].len() } else { 0 };
        let mut sum = 0.0;
        for y in 0..ny {
            for x in 0..nx {
                let s =
                    (self.p[y][x] / self.p_ref) / (self.rho[y][x] / self.rho_ref).powf(self.gamma);
                sum += s.ln().abs();
            }
        }
        sum / (nx * ny).max(1) as f64
    }
}

/// Roe flux splitting for Euler equations.
///
/// Returns approximate Riemann flux \[F_rho, F_rhou, F_rhov, F_E\].
pub fn euler_flux_roe(
    rho_l: f64,
    u_l: [f64; 2],
    p_l: f64,
    rho_r: f64,
    u_r: [f64; 2],
    p_r: f64,
) -> [f64; 4] {
    let gamma = 1.4;
    let e_l = p_l / (gamma - 1.0) + 0.5 * rho_l * (u_l[0] * u_l[0] + u_l[1] * u_l[1]);
    let e_r = p_r / (gamma - 1.0) + 0.5 * rho_r * (u_r[0] * u_r[0] + u_r[1] * u_r[1]);
    let h_l = (e_l + p_l) / rho_l;
    let h_r = (e_r + p_r) / rho_r;

    // Roe averages
    let sqrt_rho_l = rho_l.sqrt();
    let sqrt_rho_r = rho_r.sqrt();
    let denom = sqrt_rho_l + sqrt_rho_r;
    let u_roe = (sqrt_rho_l * u_l[0] + sqrt_rho_r * u_r[0]) / denom;
    let v_roe = (sqrt_rho_l * u_l[1] + sqrt_rho_r * u_r[1]) / denom;
    let h_roe = (sqrt_rho_l * h_l + sqrt_rho_r * h_r) / denom;
    let cs_roe = ((gamma - 1.0) * (h_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe)))
        .max(0.0)
        .sqrt();

    // Average flux
    let f_l = [
        rho_l * u_l[0],
        rho_l * u_l[0] * u_l[0] + p_l,
        rho_l * u_l[0] * u_l[1],
        (e_l + p_l) * u_l[0],
    ];
    let f_r = [
        rho_r * u_r[0],
        rho_r * u_r[0] * u_r[0] + p_r,
        rho_r * u_r[0] * u_r[1],
        (e_r + p_r) * u_r[0],
    ];

    // Lax-Friedrichs upwinding (simplified Roe)
    let lambda = u_roe.abs() + cs_roe;
    [
        0.5 * (f_l[0] + f_r[0]) - 0.5 * lambda * (rho_r - rho_l),
        0.5 * (f_l[1] + f_r[1]) - 0.5 * lambda * (rho_r * u_r[0] - rho_l * u_l[0]),
        0.5 * (f_l[2] + f_r[2]) - 0.5 * lambda * (rho_r * u_r[1] - rho_l * u_l[1]),
        0.5 * (f_l[3] + f_r[3]) - 0.5 * lambda * (e_r - e_l),
    ]
}

// ============================================================
// Tests
// ============================================================
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compressible_lbm_creation() {
        let sim = CompressibleLbm::new(8, 8, 1.7, 1.4, 2);
        assert_eq!(sim.nx, 8);
        assert_eq!(sim.ny, 8);
    }

    #[test]
    fn test_compressible_lbm_step() {
        let mut sim = CompressibleLbm::new(8, 8, 1.7, 1.4, 2);
        sim.step();
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                assert!(sim.rho[y][x] > 0.0);
                assert!(sim.rho[y][x].is_finite());
            }
        }
    }

    #[test]
    fn test_equilibrium_compressible_sum_rho() {
        let feq = equilibrium_compressible(1.0, 0.1, 0.05, 1.0 / 3.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_local_mach_zero_velocity() {
        let ma = local_mach(0.0, 0.0, 1.0);
        assert!(ma.abs() < 1e-15);
    }

    #[test]
    fn test_local_mach_unit() {
        let ma = local_mach(1.0, 0.0, 1.0);
        assert!((ma - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_entropy_generation_zero_shear() {
        let s = entropy_generation(1.0, 1.0, 0.0, 0.0, 1e-3);
        assert!(s.abs() < 1e-15);
    }

    #[test]
    fn test_euler_lbm_creation() {
        let sim = EulerLbm::new(8, 8, EulerLbm::DEFAULT_OMEGA, 1.4);
        assert_eq!(sim.nx, 8);
    }

    #[test]
    fn test_euler_lbm_step() {
        let mut sim = EulerLbm::new(8, 8, EulerLbm::DEFAULT_OMEGA, 1.4);
        sim.step();
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                assert!(sim.rho[y][x] > 0.0);
            }
        }
    }

    #[test]
    fn test_shock_capture_creation() {
        let sim = ShockCaptureLbm::new(8, 8, 1.7, 1.4, 0.01, 0.1);
        assert_eq!(sim.lbm.nx, 8);
    }

    #[test]
    fn test_shock_capture_step() {
        let mut sim = ShockCaptureLbm::new(8, 8, 1.7, 1.4, 0.01, 0.1);
        sim.step();
        assert!(sim.shock_cell_count() <= sim.lbm.nx * sim.lbm.ny);
    }

    #[test]
    fn test_ns_compressible_step() {
        let mut sim = NsCompressibleLbm::new(8, 8, 1.7, 1.7, 0.71, 1.4);
        sim.step();
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                assert!(sim.rho[y][x] > 0.0);
            }
        }
    }

    #[test]
    fn test_ma_correction_creation() {
        let mc = MaCorrection::new(2, 1.0 / 3.0_f64.sqrt());
        assert_eq!(mc.order, 2);
    }

    #[test]
    fn test_ma_correction_feq_sum() {
        let mc = MaCorrection::new(2, 1.0 / 3.0_f64.sqrt());
        let feq = mc.corrected_feq(1.0, 0.05, 0.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_shock_wave_1d_creation() {
        let sw = ShockWave1D::new_sod(100, 0.01, 1e-4, 1.4);
        assert_eq!(sw.n, 100);
    }

    #[test]
    fn test_shock_wave_1d_step() {
        let mut sw = ShockWave1D::new_sod(100, 0.01, 1e-5, 1.4);
        sw.run(10);
        assert!(sw.time > 0.0);
    }

    #[test]
    fn test_shock_wave_cfl() {
        let sw = ShockWave1D::new_sod(100, 0.01, 1e-5, 1.4);
        let cfl = sw.cfl();
        assert!(cfl > 0.0);
    }

    #[test]
    fn test_rayleigh_taylor_creation() {
        let rt = RayleighTaylorLbm::new(16, 32, 0.1, 1e-4).expect("valid parameters");
        assert_eq!(rt.nx, 16);
    }

    #[test]
    fn test_rayleigh_taylor_growth_rate() {
        let rt = RayleighTaylorLbm::new(16, 32, 0.1, 1e-4).expect("valid parameters");
        let gr = rt.growth_rate(1.0);
        assert!(gr > 0.0);
    }

    #[test]
    fn test_supersonic_inlet_velocity() {
        let inlet = SupersonicInlet::new(2.0, 0.0, 1.0, 1.0 / 1.4, 1.4);
        let v = inlet.velocity();
        assert!(v[0] > 0.0);
        assert!(v[1].abs() < 1e-10);
    }

    #[test]
    fn test_supersonic_inlet_total_pressure() {
        let inlet = SupersonicInlet::new(2.0, 0.0, 1.0, 1.0, 1.4);
        let pt = inlet.total_pressure();
        assert!(pt > inlet.p_in); // total > static
    }

    #[test]
    fn test_compressible_stats_mach() {
        let n = 4;
        let rho = vec![vec![1.0f64; n]; n];
        let u = vec![vec![[0.5f64, 0.0f64]; n]; n];
        let p = vec![vec![1.0 / 3.0_f64; n]; n];
        let stats = CompressibleStats::new(rho, u, p, 1.4);
        let ma = stats.mach_at(0, 0);
        assert!(ma > 0.0);
    }

    #[test]
    fn test_detonation_cj_velocity() {
        let dcj = DetonationLbm::cj_velocity(1.4, 10.0);
        assert!(dcj > 0.0);
    }

    #[test]
    fn test_detonation_step() {
        let mut det = DetonationLbm::new(50, 0.01, 1e-5, 1.4, 10.0, 5.0, 1e6);
        det.step();
        for &rho in &det.rho {
            assert!(rho > 0.0);
        }
    }

    #[test]
    fn test_euler_flux_roe_symmetry() {
        let rho = 1.0;
        let u = [0.1f64, 0.0f64];
        let p = 1.0 / 3.0;
        let flux = euler_flux_roe(rho, u, p, rho, u, p);
        // Same state on both sides: flux should be exact physical flux
        let expected_rho_flux = rho * u[0];
        assert!((flux[0] - expected_rho_flux).abs() < 1e-10);
    }

    #[test]
    fn test_max_mach_compressible() {
        let mut sim = CompressibleLbm::new(8, 8, 1.7, 1.4, 2);
        sim.u[0][0] = [0.1, 0.0];
        let ma = sim.max_mach();
        assert!(ma >= 0.0);
    }

    #[test]
    fn test_rt_bubble_velocity() {
        let rt = RayleighTaylorLbm::new(16, 32, 0.1, 1e-4).expect("valid parameters");
        let vb = rt.bubble_velocity();
        assert!(vb > 0.0);
    }
}
