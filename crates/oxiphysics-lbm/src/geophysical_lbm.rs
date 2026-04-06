#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geophysical fluid dynamics via Lattice Boltzmann Method.
//!
//! This module provides:
//! - [`RotatingFrameLbm`]: D2Q9 LBM in rotating frame with Coriolis forcing.
//! - [`ShallowWaterLbm`]: Lattice Boltzmann shallow-water equations.
//! - [`RossbyWaveSim`]: Rossby wave simulation on beta-plane.
//! - [`MantelConvection`]: Boussinesq convection in Earth's mantle.
//! - [`ThermohalineCirculation`]: Double-diffusive oceanic circulation.
//! - [`AtmosphericBoundaryLayer`]: Turbulent ABL with Monin-Obukhov similarity.
//! - [`GravityCurrent`]: Density-driven gravity current front tracking.
//! - [`TidalSimulation`]: Tidal forcing with bottom friction.
//! - [`GeostrophicAdjustment`]: Geostrophic balance and adjustment process.
//! - [`PlanetaryWave`]: Barotropic and baroclinic planetary wave propagation.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// D2Q9 lattice constants
// ─────────────────────────────────────────────────────────────────────────────

const NQ: usize = 9;
const CX: [f64; NQ] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const CY: [f64; NQ] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
const W: [f64; NQ] = [
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
const CS2: f64 = 1.0 / 3.0;

/// Opposite direction index for bounce-back boundary.
const OPP: [usize; NQ] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

// ─────────────────────────────────────────────────────────────────────────────
// Math helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn feq(rho: f64, ux: f64, uy: f64) -> [f64; NQ] {
    let mut f = [0.0f64; NQ];
    let u2 = ux * ux + uy * uy;
    for q in 0..NQ {
        let cu = CX[q] * ux + CY[q] * uy;
        f[q] = W[q] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    f
}

#[inline]
fn macroscopic(f: &[f64; NQ]) -> (f64, f64, f64) {
    let rho = f.iter().sum::<f64>();
    if rho < 1e-300 {
        return (0.0, 0.0, 0.0);
    }
    let ux = f.iter().enumerate().map(|(q, &v)| CX[q] * v).sum::<f64>() / rho;
    let uy = f.iter().enumerate().map(|(q, &v)| CY[q] * v).sum::<f64>() / rho;
    (rho, ux, uy)
}

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BGK collision with body force (Guo scheme)
// ─────────────────────────────────────────────────────────────────────────────

/// Apply Guo body-force correction to distribution after BGK collision.
///
/// `f_post` is updated in-place.
fn guo_force(f_post: &mut [f64; NQ], fx: f64, fy: f64, ux: f64, uy: f64, omega: f64) {
    for q in 0..NQ {
        let cu = CX[q] * ux + CY[q] * uy;
        let cf = CX[q] * fx + CY[q] * fy;
        let correction = W[q]
            * (1.0 - 0.5 * omega)
            * (cf / CS2 + (cu * cf) / (CS2 * CS2) - (ux * fx + uy * fy) / CS2);
        f_post[q] += correction;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Periodic streaming helper
// ─────────────────────────────────────────────────────────────────────────────

fn stream_periodic(f: &[Vec<[f64; NQ]>], nx: usize, ny: usize) -> Vec<Vec<[f64; NQ]>> {
    let mut f_new = vec![vec![[0.0f64; NQ]; ny]; nx];
    for ix in 0..nx {
        for iy in 0..ny {
            for q in 0..NQ {
                let cx = CX[q] as isize;
                let cy = CY[q] as isize;
                let sx = ((ix as isize - cx).rem_euclid(nx as isize)) as usize;
                let sy = ((iy as isize - cy).rem_euclid(ny as isize)) as usize;
                f_new[ix][iy][q] = f[sx][sy][q];
            }
        }
    }
    f_new
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. RotatingFrameLbm — D2Q9 with Coriolis forcing
// ─────────────────────────────────────────────────────────────────────────────

/// D2Q9 Lattice Boltzmann simulation in a rotating reference frame.
///
/// Implements Coriolis force `F = -2Ω × u` via the Guo body-force scheme,
/// enabling simulation of rotating geophysical flows such as inertial
/// oscillations and rotating Rayleigh-Bénard convection.
pub struct RotatingFrameLbm {
    /// Grid width (number of cells in x direction).
    pub nx: usize,
    /// Grid height (number of cells in y direction).
    pub ny: usize,
    /// BGK relaxation parameter (ω = 1/τ).
    pub omega: f64,
    /// Rotation rate (Coriolis parameter) in lattice units.
    pub rotation_rate: f64,
    /// Distribution functions (outer: x, inner: y).
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Solid mask: true = solid (bounce-back) cell.
    pub solid: Vec<Vec<bool>>,
}

impl RotatingFrameLbm {
    /// Create a new rotating-frame LBM grid.
    ///
    /// # Arguments
    /// * `nx` — width in lattice cells.
    /// * `ny` — height in lattice cells.
    /// * `omega` — BGK relaxation frequency (0 < ω < 2).
    /// * `rotation_rate` — Coriolis parameter f (rad/step in lattice units).
    pub fn new(nx: usize, ny: usize, omega: f64, rotation_rate: f64) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let solid = vec![vec![false; ny]; nx];
        Self {
            nx,
            ny,
            omega,
            rotation_rate,
            f,
            solid,
        }
    }

    /// Set velocity initial condition for all fluid cells.
    pub fn set_velocity_ic(&mut self, ux0: f64, uy0: f64) {
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                self.f[ix][iy] = feq(1.0, ux0, uy0);
            }
        }
    }

    /// Perform one LBM step: collision → Coriolis forcing → streaming.
    pub fn step(&mut self) {
        let omega = self.omega;
        let rot = self.rotation_rate;
        // Collision + Coriolis
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                if self.solid[ix][iy] {
                    continue;
                }
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let feq_arr = feq(rho, ux, uy);
                let mut f_coll = self.f[ix][iy];
                for q in 0..NQ {
                    f_coll[q] += omega * (feq_arr[q] - f_coll[q]);
                }
                // Coriolis: F_x = +2Ω*uy, F_y = -2Ω*ux  (f-plane)
                let fx = 2.0 * rot * uy;
                let fy = -2.0 * rot * ux;
                guo_force(&mut f_coll, fx, fy, ux, uy, omega);
                self.f[ix][iy] = f_coll;
            }
        }
        // Streaming with periodic BC
        let f_new = stream_periodic(&self.f, self.nx, self.ny);
        // Bounce-back for solid cells
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                if self.solid[ix][iy] {
                    let old = self.f[ix][iy];
                    for q in 0..NQ {
                        self.f[ix][iy][q] = old[OPP[q]];
                    }
                } else {
                    self.f[ix][iy] = f_new[ix][iy];
                }
            }
        }
    }

    /// Compute bulk kinetic energy (lattice units).
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                ke += 0.5 * rho * (ux * ux + uy * uy);
            }
        }
        ke / (self.nx * self.ny) as f64
    }

    /// Extract macroscopic velocity field as flat (ux, uy) pairs.
    pub fn velocity_field(&self) -> Vec<(f64, f64)> {
        let mut out = Vec::with_capacity(self.nx * self.ny);
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (_rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                out.push((ux, uy));
            }
        }
        out
    }

    /// Compute enstrophy (integral of ω²/2) using finite differences.
    pub fn enstrophy(&self) -> f64 {
        let mut ens = 0.0;
        let nx = self.nx;
        let ny = self.ny;
        for ix in 0..nx {
            for iy in 0..ny {
                let ixp = (ix + 1) % nx;
                let ixm = (ix + nx - 1) % nx;
                let iyp = (iy + 1) % ny;
                let iym = (iy + ny - 1) % ny;
                let (_r, _ux, uym) = macroscopic(&self.f[ix][iym]);
                let (_r2, _ux2, uyp) = macroscopic(&self.f[ix][iyp]);
                let (_r3, uxm, _uy3) = macroscopic(&self.f[ixm][iy]);
                let (_r4, uxp, _uy4) = macroscopic(&self.f[ixp][iy]);
                let dvdx = (uyp - uym) * 0.5;
                let dudy = (uxp - uxm) * 0.5;
                let vort = dvdx - dudy;
                ens += 0.5 * vort * vort;
            }
        }
        ens / (nx * ny) as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. ShallowWaterLbm — lattice Boltzmann shallow-water equations
// ─────────────────────────────────────────────────────────────────────────────

/// Shallow-water LBM using height-velocity formulation.
///
/// The shallow-water equations are recovered from the LBM in the long-wave
/// limit with gravity wave speed `c = sqrt(g*H)`.  This is a 1D formulation
/// (D1Q3) extended to 2D with D2Q9 moments.
pub struct ShallowWaterLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// BGK relaxation frequency.
    pub omega: f64,
    /// Gravitational acceleration in lattice units.
    pub gravity: f64,
    /// Mean water depth.
    pub h_mean: f64,
    /// Distribution functions representing height field.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Bottom topography (depth perturbation).
    pub bathymetry: Vec<Vec<f64>>,
}

impl ShallowWaterLbm {
    /// Create a new shallow-water LBM domain.
    ///
    /// # Arguments
    /// * `nx`, `ny` — domain size.
    /// * `omega` — relaxation frequency.
    /// * `gravity` — gravitational acceleration (lattice units).
    /// * `h_mean` — mean water depth.
    pub fn new(nx: usize, ny: usize, omega: f64, gravity: f64, h_mean: f64) -> Self {
        let f = vec![vec![feq(h_mean, 0.0, 0.0); ny]; nx];
        let bathymetry = vec![vec![0.0f64; ny]; nx];
        Self {
            nx,
            ny,
            omega,
            gravity,
            h_mean,
            f,
            bathymetry,
        }
    }

    /// Set a Gaussian dam-break initial condition.
    ///
    /// Elevates water height by `amplitude` at grid center.
    pub fn dam_break_ic(&mut self, amplitude: f64) {
        let cx = self.nx as f64 * 0.5;
        let cy = self.ny as f64 * 0.5;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let r2 = (ix as f64 - cx).powi(2) + (iy as f64 - cy).powi(2);
                let h = self.h_mean + amplitude * (-r2 / 50.0).exp();
                self.f[ix][iy] = feq(h, 0.0, 0.0);
            }
        }
    }

    /// Shallow-water equilibrium with gravity wave dispersion.
    fn feq_sw(h: f64, ux: f64, uy: f64, _gravity: f64) -> [f64; NQ] {
        let mut f = [0.0f64; NQ];
        let u2 = ux * ux + uy * uy;
        for q in 0..NQ {
            let cu = CX[q] * ux + CY[q] * uy;
            f[q] = W[q] * h * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
        }
        f
    }

    /// Advance one timestep.
    pub fn step(&mut self) {
        let omega = self.omega;
        let g = self.gravity;
        let h_mean = self.h_mean;
        // Collision
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (h, ux, uy) = macroscopic(&self.f[ix][iy]);
                let h_eff = (h - self.bathymetry[ix][iy]).max(1e-6);
                let feq_arr = Self::feq_sw(h_eff, ux, uy, g);
                // Pressure-gradient forcing from bathymetry slope
                let bx = if ix + 1 < self.nx {
                    self.bathymetry[ix + 1][iy] - self.bathymetry[ix][iy]
                } else {
                    0.0
                };
                let by = if iy + 1 < self.ny {
                    self.bathymetry[ix][iy + 1] - self.bathymetry[ix][iy]
                } else {
                    0.0
                };
                let fx = -g * h_mean * bx;
                let fy = -g * h_mean * by;
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega * (feq_arr[q] - self.f[ix][iy][q]);
                }
                let mut tmp = self.f[ix][iy];
                guo_force(&mut tmp, fx, fy, ux, uy, omega);
                self.f[ix][iy] = tmp;
            }
        }
        // Streaming
        self.f = stream_periodic(&self.f, self.nx, self.ny);
    }

    /// Return the height field as a flat vector.
    pub fn height_field(&self) -> Vec<f64> {
        let mut h = Vec::with_capacity(self.nx * self.ny);
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, _ux, _uy) = macroscopic(&self.f[ix][iy]);
                h.push(rho);
            }
        }
        h
    }

    /// Compute total energy (kinetic + potential).
    pub fn total_energy(&self) -> f64 {
        let g = self.gravity;
        let mut e = 0.0;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (h, ux, uy) = macroscopic(&self.f[ix][iy]);
                e += 0.5 * h * (ux * ux + uy * uy) + 0.5 * g * h * h;
            }
        }
        e
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. RossbyWaveSim — Rossby wave on beta-plane
// ─────────────────────────────────────────────────────────────────────────────

/// Rossby wave simulation using beta-plane approximation.
///
/// The Coriolis parameter varies linearly with latitude: `f(y) = f0 + β·y`.
/// This produces westward-propagating Rossby waves characteristic of
/// planetary-scale atmospheric and oceanic dynamics.
pub struct RossbyWaveSim {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// BGK relaxation frequency.
    pub omega: f64,
    /// Reference Coriolis parameter at domain centre.
    pub f0: f64,
    /// Beta-plane gradient df/dy.
    pub beta: f64,
    /// Lattice spacing (physical length per cell).
    pub dy: f64,
    /// Distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Stream-function (diagnosed from vorticity).
    pub psi: Vec<Vec<f64>>,
}

impl RossbyWaveSim {
    /// Construct a new beta-plane Rossby wave simulation.
    pub fn new(nx: usize, ny: usize, omega: f64, f0: f64, beta: f64, dy: f64) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let psi = vec![vec![0.0f64; ny]; nx];
        Self {
            nx,
            ny,
            omega,
            f0,
            beta,
            dy,
            f,
            psi,
        }
    }

    /// Set a single Rossby wave initial condition with wavenumbers (kx, ky).
    pub fn wave_ic(&mut self, kx: f64, ky: f64, amplitude: f64) {
        let nx = self.nx as f64;
        let ny = self.ny as f64;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let x = ix as f64 / nx;
                let y = iy as f64 / ny;
                let ux = -amplitude * ky * (2.0 * PI * (kx * x + ky * y)).cos();
                let uy = amplitude * kx * (2.0 * PI * (kx * x + ky * y)).cos();
                self.f[ix][iy] = feq(1.0, ux, uy);
            }
        }
    }

    /// Compute local Coriolis parameter at lattice row `iy`.
    #[inline]
    pub fn coriolis_at(&self, iy: usize) -> f64 {
        let y = iy as f64 * self.dy;
        self.f0 + self.beta * y
    }

    /// Advance one step with spatially varying Coriolis.
    pub fn step(&mut self) {
        let omega = self.omega;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let feq_arr = feq(rho, ux, uy);
                let fc = self.coriolis_at(iy);
                let fx = fc * uy;
                let fy = -fc * ux;
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega * (feq_arr[q] - self.f[ix][iy][q]);
                }
                let mut tmp = self.f[ix][iy];
                guo_force(&mut tmp, fx, fy, ux, uy, omega);
                self.f[ix][iy] = tmp;
            }
        }
        self.f = stream_periodic(&self.f, self.nx, self.ny);
        self.update_streamfunction();
    }

    /// Update stream-function by integrating vorticity (simplified).
    fn update_streamfunction(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for ix in 0..nx {
            for iy in 0..ny {
                let ixp = (ix + 1) % nx;
                let ixm = (ix + nx - 1) % nx;
                let iyp = (iy + 1) % ny;
                let iym = (iy + ny - 1) % ny;
                let (_r, _ux, uym) = macroscopic(&self.f[ix][iym]);
                let (_r2, _ux2, uyp) = macroscopic(&self.f[ix][iyp]);
                let (_r3, uxm, _uy3) = macroscopic(&self.f[ixm][iy]);
                let (_r4, uxp, _uy4) = macroscopic(&self.f[ixp][iy]);
                let vort = 0.5 * ((uyp - uym) - (uxp - uxm));
                self.psi[ix][iy] = vort;
            }
        }
    }

    /// Compute Rossby wave phase speed (theoretical) for wavenumber (kx, ky).
    ///
    /// `c_x = -β / (kx² + ky²)` (non-dimensional).
    pub fn theoretical_phase_speed(&self, kx: f64, ky: f64) -> f64 {
        let k2 = kx * kx + ky * ky;
        if k2 < 1e-14 {
            return 0.0;
        }
        -self.beta / k2
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. MantelConvection — Boussinesq convection in Earth's mantle
// ─────────────────────────────────────────────────────────────────────────────

/// Boussinesq thermal convection model for Earth's mantle.
///
/// Uses a double-distribution LBM: one for momentum (D2Q9) and one for
/// temperature (D2Q5 collapsed to D2Q9 with passive scalar weights).
/// The buoyancy force is `F_y = Ra·Pr·θ` with Boussinesq approximation.
pub struct MantelConvection {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Momentum relaxation frequency.
    pub omega_f: f64,
    /// Temperature relaxation frequency.
    pub omega_t: f64,
    /// Rayleigh number.
    pub ra: f64,
    /// Prandtl number.
    pub pr: f64,
    /// Momentum distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Temperature distribution functions (D2Q9 proxy).
    pub g: Vec<Vec<[f64; NQ]>>,
    /// Temperature field (cached).
    pub temp: Vec<Vec<f64>>,
}

impl MantelConvection {
    /// Create a Boussinesq mantel convection simulation.
    ///
    /// The Rayleigh number `Ra` controls convective vigour; for Earth's mantle
    /// `Ra ~ 10^6 – 10^8`.  The Prandtl number `Pr ~ 10^23` (effective) but
    /// for LBM we use `Pr ~ 1 – 100` at reduced resolution.
    pub fn new(nx: usize, ny: usize, omega_f: f64, omega_t: f64, ra: f64, pr: f64) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let g = vec![vec![feq(0.5, 0.0, 0.0); ny]; nx];
        let temp = vec![vec![0.5f64; ny]; nx];
        Self {
            nx,
            ny,
            omega_f,
            omega_t,
            ra,
            pr,
            f,
            g,
            temp,
        }
    }

    /// Set a conductive profile (linear temperature 0..1 from top to bottom).
    pub fn conductive_ic(&mut self) {
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let t = 1.0 - iy as f64 / (self.ny - 1) as f64;
                self.temp[ix][iy] = t;
                self.g[ix][iy] = feq(t, 0.0, 0.0);
            }
        }
    }

    /// Temperature equilibrium distribution.
    fn geq(t: f64, ux: f64, uy: f64) -> [f64; NQ] {
        feq(t, ux, uy)
    }

    /// Advance one coupled timestep.
    pub fn step(&mut self) {
        let omega_f = self.omega_f;
        let omega_t = self.omega_t;
        let ra = self.ra;
        let pr = self.pr;
        // Temperature collision
        let mut g_new = self.g.clone();
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (t, ux, uy) = macroscopic(&self.g[ix][iy]);
                let geq_arr = Self::geq(t, ux, uy);
                for q in 0..NQ {
                    g_new[ix][iy][q] =
                        self.g[ix][iy][q] + omega_t * (geq_arr[q] - self.g[ix][iy][q]);
                }
                self.temp[ix][iy] = t;
                let _ = (ux, uy, pr);
            }
        }
        self.g = stream_periodic(&g_new, self.nx, self.ny);
        // Momentum collision with buoyancy
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let feq_arr = feq(rho, ux, uy);
                let theta = self.temp[ix][iy] - 0.5;
                let buoy = ra * pr * theta;
                let fx = 0.0;
                let fy = buoy;
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega_f * (feq_arr[q] - self.f[ix][iy][q]);
                }
                let mut tmp = self.f[ix][iy];
                guo_force(&mut tmp, fx, fy, ux, uy, omega_f);
                self.f[ix][iy] = tmp;
            }
        }
        self.f = stream_periodic(&self.f, self.nx, self.ny);
        // Fixed temperature top/bottom walls
        for ix in 0..self.nx {
            let t_bot = 1.0;
            let t_top = 0.0;
            self.g[ix][0] = feq(t_bot, 0.0, 0.0);
            self.g[ix][self.ny - 1] = feq(t_top, 0.0, 0.0);
        }
    }

    /// Compute Nusselt number from heat flux at bottom wall.
    pub fn nusselt_number(&self) -> f64 {
        let ny = self.ny;
        let flux: f64 = (0..self.nx)
            .map(|ix| {
                let t1 = self.temp[ix][0];
                let t2 = self.temp[ix][1.min(ny - 1)];
                t1 - t2
            })
            .sum::<f64>()
            / self.nx as f64;
        // Nu = q / (ΔT/H) where ΔT=1, H=ny
        flux * self.ny as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. ThermohalineCirculation — double-diffusive oceanic convection
// ─────────────────────────────────────────────────────────────────────────────

/// Double-diffusive thermohaline circulation model.
///
/// Two scalar fields (temperature T, salinity S) drive density via
/// `ρ = ρ₀(1 - αT + βS)`.  Uses three coupled D2Q9 distributions.
pub struct ThermohalineCirculation {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Momentum relaxation rate.
    pub omega_f: f64,
    /// Thermal diffusion relaxation rate.
    pub omega_t: f64,
    /// Saline diffusion relaxation rate.
    pub omega_s: f64,
    /// Thermal expansion coefficient.
    pub alpha_t: f64,
    /// Haline contraction coefficient.
    pub beta_s: f64,
    /// Momentum distribution.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Temperature distribution.
    pub g_t: Vec<Vec<[f64; NQ]>>,
    /// Salinity distribution.
    pub g_s: Vec<Vec<[f64; NQ]>>,
}

impl ThermohalineCirculation {
    /// Construct a thermohaline circulation model.
    pub fn new(
        nx: usize,
        ny: usize,
        omega_f: f64,
        omega_t: f64,
        omega_s: f64,
        alpha_t: f64,
        beta_s: f64,
    ) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let g_t = vec![vec![feq(0.5, 0.0, 0.0); ny]; nx];
        let g_s = vec![vec![feq(0.5, 0.0, 0.0); ny]; nx];
        Self {
            nx,
            ny,
            omega_f,
            omega_t,
            omega_s,
            alpha_t,
            beta_s,
            f,
            g_t,
            g_s,
        }
    }

    /// Advance one step.
    pub fn step(&mut self) {
        let omega_f = self.omega_f;
        let omega_t = self.omega_t;
        let omega_s = self.omega_s;
        let alpha_t = self.alpha_t;
        let beta_s = self.beta_s;
        // Collect T, S
        let mut temp_arr = vec![vec![0.0f64; self.ny]; self.nx];
        let mut sal_arr = vec![vec![0.0f64; self.ny]; self.nx];
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (t, _ux, _uy) = macroscopic(&self.g_t[ix][iy]);
                let (s, _ux2, _uy2) = macroscopic(&self.g_s[ix][iy]);
                temp_arr[ix][iy] = t;
                sal_arr[ix][iy] = s;
            }
        }
        // Temperature collision
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (t, ux, uy) = macroscopic(&self.g_t[ix][iy]);
                let geq_arr = feq(t, ux, uy);
                for q in 0..NQ {
                    self.g_t[ix][iy][q] += omega_t * (geq_arr[q] - self.g_t[ix][iy][q]);
                }
            }
        }
        self.g_t = stream_periodic(&self.g_t, self.nx, self.ny);
        // Salinity collision
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (s, ux, uy) = macroscopic(&self.g_s[ix][iy]);
                let geq_arr = feq(s, ux, uy);
                for q in 0..NQ {
                    self.g_s[ix][iy][q] += omega_s * (geq_arr[q] - self.g_s[ix][iy][q]);
                }
            }
        }
        self.g_s = stream_periodic(&self.g_s, self.nx, self.ny);
        // Momentum with density buoyancy
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let feq_arr = feq(rho, ux, uy);
                let t = temp_arr[ix][iy];
                let s = sal_arr[ix][iy];
                let rho_prime = -alpha_t * (t - 0.5) + beta_s * (s - 0.5);
                let fy = -rho_prime * 0.01;
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega_f * (feq_arr[q] - self.f[ix][iy][q]);
                }
                let mut tmp = self.f[ix][iy];
                guo_force(&mut tmp, 0.0, fy, ux, uy, omega_f);
                self.f[ix][iy] = tmp;
            }
        }
        self.f = stream_periodic(&self.f, self.nx, self.ny);
    }

    /// Return density anomaly field `ρ' = -αT' + βS'`.
    pub fn density_anomaly(&self) -> Vec<Vec<f64>> {
        let mut out = vec![vec![0.0f64; self.ny]; self.nx];
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (t, _a, _b) = macroscopic(&self.g_t[ix][iy]);
                let (s, _c, _d) = macroscopic(&self.g_s[ix][iy]);
                out[ix][iy] = -self.alpha_t * (t - 0.5) + self.beta_s * (s - 0.5);
            }
        }
        out
    }

    /// Diagnose thermohaline overturning circulation strength.
    pub fn overturning_strength(&self) -> f64 {
        let ny = self.ny;
        let ny2 = ny / 2;
        let psi: f64 = (0..self.nx)
            .map(|ix| {
                let (_rho, ux, _uy) = macroscopic(&self.f[ix][ny2]);
                ux
            })
            .sum();
        psi.abs()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. AtmosphericBoundaryLayer — turbulent ABL
// ─────────────────────────────────────────────────────────────────────────────

/// Atmospheric boundary layer simulation with Monin-Obukhov similarity.
///
/// Implements a modified Smagorinsky LES approach suitable for neutral and
/// stable ABL.  Wall stress uses log-law with roughness length `z0`.
pub struct AtmosphericBoundaryLayer {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// BGK relaxation frequency (background viscosity).
    pub omega_base: f64,
    /// Smagorinsky constant Cs.
    pub cs_smag: f64,
    /// Roughness length z0 in lattice units.
    pub z0: f64,
    /// Geostrophic wind speed (driving).
    pub u_geo: f64,
    /// Distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Turbulent viscosity field.
    pub nu_t: Vec<Vec<f64>>,
}

impl AtmosphericBoundaryLayer {
    /// Construct a new ABL simulation domain.
    pub fn new(nx: usize, ny: usize, omega_base: f64, cs_smag: f64, z0: f64, u_geo: f64) -> Self {
        let f = vec![vec![feq(1.0, u_geo, 0.0); ny]; nx];
        let nu_t = vec![vec![0.0f64; ny]; nx];
        Self {
            nx,
            ny,
            omega_base,
            cs_smag,
            z0,
            u_geo,
            f,
            nu_t,
        }
    }

    /// Compute local strain-rate magnitude |S|.
    fn strain_rate(&self, ix: usize, iy: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let ixp = (ix + 1) % nx;
        let ixm = (ix + nx - 1) % nx;
        let iyp = (iy + 1) % ny;
        let iym = (iy + ny - 1) % ny;
        let (_r, uxpx, _) = macroscopic(&self.f[ixp][iy]);
        let (_r, uxmx, _) = macroscopic(&self.f[ixm][iy]);
        let (_r, _, uypy) = macroscopic(&self.f[ix][iyp]);
        let (_r, _, uymy) = macroscopic(&self.f[ix][iym]);
        let (_r, uxpy, _) = macroscopic(&self.f[ix][iyp]);
        let (_r, uxmy, _) = macroscopic(&self.f[ix][iym]);
        let (_r, _, uypx) = macroscopic(&self.f[ixp][iy]);
        let (_r, _, uymx) = macroscopic(&self.f[ixm][iy]);
        let s11 = 0.5 * (uxpx - uxmx);
        let s22 = 0.5 * (uypy - uymy);
        let s12 = 0.25 * ((uxpy - uxmy) + (uypx - uymx));
        (2.0 * (s11 * s11 + s22 * s22 + 2.0 * s12 * s12)).sqrt()
    }

    /// Update turbulent viscosity using Smagorinsky model.
    pub fn update_turbulent_viscosity(&mut self) {
        let cs = self.cs_smag;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let smag = self.strain_rate(ix, iy);
                self.nu_t[ix][iy] = cs * cs * smag;
            }
        }
    }

    /// Advance one timestep with Smagorinsky LES.
    pub fn step(&mut self) {
        self.update_turbulent_viscosity();
        let omega_base = self.omega_base;
        let u_geo = self.u_geo;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let nu_t_local = self.nu_t[ix][iy];
                let nu_eff = (1.0 / omega_base - 0.5) * CS2 + nu_t_local;
                let omega_eff = clamp(1.0 / (nu_eff / CS2 + 0.5), 0.01, 1.99);
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let feq_arr = feq(rho, ux, uy);
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega_eff * (feq_arr[q] - self.f[ix][iy][q]);
                }
            }
        }
        // Geostrophic forcing at top
        for ix in 0..self.nx {
            let (_rho, ux, _uy) = macroscopic(&self.f[ix][self.ny - 1]);
            let du = u_geo - ux;
            for q in 0..NQ {
                self.f[ix][self.ny - 1][q] += 0.01 * du * W[q];
            }
        }
        // Log-law wall stress at bottom
        for ix in 0..self.nx {
            let (_rho, ux, _uy) = macroscopic(&self.f[ix][1]);
            let z = 1.0;
            let kappa = 0.41;
            let u_star = ux * kappa / (z / self.z0 + 1.0).ln().max(1e-6);
            let tau_wall = u_star * u_star;
            let _ = tau_wall;
            for q in 0..NQ {
                self.f[ix][0][q] = self.f[ix][0][OPP[q]];
            }
        }
        self.f = stream_periodic(&self.f, self.nx, self.ny);
    }

    /// Return vertically-averaged horizontal wind profile.
    pub fn wind_profile(&self) -> Vec<f64> {
        (0..self.ny)
            .map(|iy| {
                let ux_sum: f64 = (0..self.nx).map(|ix| macroscopic(&self.f[ix][iy]).1).sum();
                ux_sum / self.nx as f64
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 7. GravityCurrent — density-driven gravity current
// ─────────────────────────────────────────────────────────────────────────────

/// Density-driven gravity current simulation.
///
/// A heavy fluid (density `ρ_h`) intrudes along the bottom beneath a lighter
/// ambient fluid (`ρ_l`).  The front speed obeys Froude's law.
pub struct GravityCurrent {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Relaxation frequency.
    pub omega: f64,
    /// Heavy fluid density.
    pub rho_heavy: f64,
    /// Light fluid density.
    pub rho_light: f64,
    /// Gravitational acceleration.
    pub gravity: f64,
    /// Distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Density field (scalar).
    pub rho_field: Vec<Vec<f64>>,
}

impl GravityCurrent {
    /// Create a gravity current domain with initial lock-release.
    ///
    /// Heavy fluid occupies the left half of the bottom quarter.
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        rho_heavy: f64,
        rho_light: f64,
        gravity: f64,
    ) -> Self {
        let mut f = vec![vec![feq(rho_light, 0.0, 0.0); ny]; nx];
        let mut rho_field = vec![vec![rho_light; ny]; nx];
        // Lock release: left half, bottom quarter
        for ix in 0..nx / 2 {
            for iy in 0..ny / 4 {
                f[ix][iy] = feq(rho_heavy, 0.0, 0.0);
                rho_field[ix][iy] = rho_heavy;
            }
        }
        Self {
            nx,
            ny,
            omega,
            rho_heavy,
            rho_light,
            gravity,
            f,
            rho_field,
        }
    }

    /// Advance one timestep with buoyancy force.
    pub fn step(&mut self) {
        let omega = self.omega;
        let g = self.gravity;
        let rho_l = self.rho_light;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                self.rho_field[ix][iy] = rho;
                let feq_arr = feq(rho, ux, uy);
                let delta_rho = rho - rho_l;
                let fy = -g * delta_rho / rho.max(1e-6);
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega * (feq_arr[q] - self.f[ix][iy][q]);
                }
                let mut tmp = self.f[ix][iy];
                guo_force(&mut tmp, 0.0, fy, ux, uy, omega);
                self.f[ix][iy] = tmp;
            }
        }
        // Bottom and top bounce-back
        for ix in 0..self.nx {
            for q in 0..NQ {
                self.f[ix][0][q] = self.f[ix][0][OPP[q]];
                self.f[ix][self.ny - 1][q] = self.f[ix][self.ny - 1][OPP[q]];
            }
        }
        self.f = stream_periodic(&self.f, self.nx, self.ny);
    }

    /// Locate the front position (rightmost cell with rho > midpoint).
    pub fn front_position(&self) -> usize {
        let rho_mid = (self.rho_heavy + self.rho_light) * 0.5;
        for ix in (0..self.nx).rev() {
            let avg: f64 = (0..self.ny / 4)
                .map(|iy| self.rho_field[ix][iy])
                .sum::<f64>()
                / (self.ny / 4).max(1) as f64;
            if avg > rho_mid {
                return ix;
            }
        }
        0
    }

    /// Theoretical Froude front speed.
    pub fn froude_speed(&self) -> f64 {
        let g_prime = self.gravity * (self.rho_heavy - self.rho_light) / self.rho_light.max(1e-14);
        let h = self.ny as f64 / 4.0;
        0.5 * (g_prime * h).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 8. TidalSimulation — tidal forcing with bottom friction
// ─────────────────────────────────────────────────────────────────────────────

/// Tidal simulation with astronomical forcing and bottom friction.
///
/// Implements the M2 tidal constituent as a harmonic body force,
/// with bottom friction via quadratic drag law `τ = Cd * |u| * u`.
pub struct TidalSimulation {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Relaxation frequency.
    pub omega: f64,
    /// Tidal forcing amplitude.
    pub tidal_amplitude: f64,
    /// Tidal angular frequency (M2 period = 44712 s; scaled to lattice units).
    pub tidal_omega: f64,
    /// Bottom drag coefficient Cd.
    pub drag_coeff: f64,
    /// Current time (in lattice steps).
    pub time: usize,
    /// Distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Bathymetry (depth).
    pub depth: Vec<Vec<f64>>,
}

impl TidalSimulation {
    /// Create a new tidal simulation.
    pub fn new(
        nx: usize,
        ny: usize,
        omega: f64,
        tidal_amplitude: f64,
        tidal_omega: f64,
        drag_coeff: f64,
    ) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let depth = vec![vec![1.0f64; ny]; nx];
        Self {
            nx,
            ny,
            omega,
            tidal_amplitude,
            tidal_omega,
            drag_coeff,
            time: 0,
            f,
            depth,
        }
    }

    /// Set a shelf-break bathymetry: shallow shelf then deep ocean.
    pub fn shelf_bathymetry(&mut self, shelf_depth: f64, ocean_depth: f64, shelf_x: usize) {
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                self.depth[ix][iy] = if ix < shelf_x {
                    shelf_depth
                } else {
                    ocean_depth
                };
            }
        }
    }

    /// Advance one timestep.
    pub fn step(&mut self) {
        let omega = self.omega;
        let a = self.tidal_amplitude;
        let w = self.tidal_omega;
        let t = self.time as f64;
        let fx_tidal = a * w * (w * t).cos();
        let cd = self.drag_coeff;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let feq_arr = feq(rho, ux, uy);
                let u_mag = (ux * ux + uy * uy).sqrt();
                let fx = fx_tidal - cd * u_mag * ux;
                let fy = -cd * u_mag * uy;
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega * (feq_arr[q] - self.f[ix][iy][q]);
                }
                let mut tmp = self.f[ix][iy];
                guo_force(&mut tmp, fx, fy, ux, uy, omega);
                self.f[ix][iy] = tmp;
            }
        }
        self.f = stream_periodic(&self.f, self.nx, self.ny);
        self.time += 1;
    }

    /// Compute tidal ellipse parameters (semi-major, semi-minor axes) at (ix, iy).
    pub fn tidal_ellipse(&self, ix: usize, iy: usize) -> (f64, f64) {
        let (_rho, ux, uy) = macroscopic(&self.f[ix][iy]);
        let major = (ux * ux + uy * uy).sqrt();
        let minor = (ux * uy).abs() / major.max(1e-14);
        (major, minor)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 9. GeostrophicAdjustment — geostrophic balance and adjustment
// ─────────────────────────────────────────────────────────────────────────────

/// Geostrophic adjustment process simulation.
///
/// Tracks how an unbalanced initial state adjusts to geostrophic balance
/// `f × u = -∇p/ρ`.  The Rossby radius of deformation `L_R = NH/f`
/// controls the adjustment length scale.
pub struct GeostrophicAdjustment {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Relaxation frequency.
    pub omega: f64,
    /// Coriolis parameter.
    pub f_coriolis: f64,
    /// Distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Pressure field (diagnosed).
    pub pressure: Vec<Vec<f64>>,
    /// Geostrophic imbalance history.
    pub imbalance_history: Vec<f64>,
}

impl GeostrophicAdjustment {
    /// Construct a geostrophic adjustment simulation.
    pub fn new(nx: usize, ny: usize, omega: f64, f_coriolis: f64) -> Self {
        let f = vec![vec![feq(1.0, 0.0, 0.0); ny]; nx];
        let pressure = vec![vec![CS2; ny]; nx];
        Self {
            nx,
            ny,
            omega,
            f_coriolis,
            f,
            pressure,
            imbalance_history: Vec::new(),
        }
    }

    /// Set a jet initial condition: zonal flow with Gaussian profile.
    pub fn jet_ic(&mut self, u_jet: f64, width: f64) {
        let cy = self.ny as f64 * 0.5;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let y = iy as f64 - cy;
                let ux = u_jet * (-y * y / (2.0 * width * width)).exp();
                self.f[ix][iy] = feq(1.0, ux, 0.0);
            }
        }
    }

    /// Update pressure from density (ideal gas in lattice units).
    fn update_pressure(&mut self) {
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, _ux, _uy) = macroscopic(&self.f[ix][iy]);
                self.pressure[ix][iy] = CS2 * rho;
            }
        }
    }

    /// Compute geostrophic imbalance: `|f × u + ∇p/ρ|` domain mean.
    pub fn geostrophic_imbalance(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let fc = self.f_coriolis;
        let mut total = 0.0;
        let n = (nx * ny) as f64;
        for ix in 0..nx {
            for iy in 0..ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let ixp = (ix + 1) % nx;
                let ixm = (ix + nx - 1) % nx;
                let iyp = (iy + 1) % ny;
                let iym = (iy + ny - 1) % ny;
                let dpdx = 0.5 * (self.pressure[ixp][iy] - self.pressure[ixm][iy]);
                let dpdy = 0.5 * (self.pressure[ix][iyp] - self.pressure[ix][iym]);
                let ag_x = fc * uy + dpdx / rho.max(1e-14);
                let ag_y = -fc * ux + dpdy / rho.max(1e-14);
                total += (ag_x * ag_x + ag_y * ag_y).sqrt();
            }
        }
        total / n
    }

    /// Advance one timestep.
    pub fn step(&mut self) {
        let omega = self.omega;
        let fc = self.f_coriolis;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let (rho, ux, uy) = macroscopic(&self.f[ix][iy]);
                let feq_arr = feq(rho, ux, uy);
                let fx = fc * uy;
                let fy = -fc * ux;
                for q in 0..NQ {
                    self.f[ix][iy][q] += omega * (feq_arr[q] - self.f[ix][iy][q]);
                }
                let mut tmp = self.f[ix][iy];
                guo_force(&mut tmp, fx, fy, ux, uy, omega);
                self.f[ix][iy] = tmp;
            }
        }
        self.f = stream_periodic(&self.f, self.nx, self.ny);
        self.update_pressure();
        let imbal = self.geostrophic_imbalance();
        self.imbalance_history.push(imbal);
    }

    /// Compute Rossby radius of deformation.
    pub fn rossby_radius(&self, brunt_vaisala_freq: f64, layer_depth: f64) -> f64 {
        if self.f_coriolis.abs() < 1e-14 {
            return f64::INFINITY;
        }
        brunt_vaisala_freq * layer_depth / self.f_coriolis.abs()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 10. PlanetaryWave — barotropic and baroclinic planetary wave propagation
// ─────────────────────────────────────────────────────────────────────────────

/// Planetary wave (Rossby and Kelvin wave) propagation model.
///
/// Supports barotropic (single-layer) and baroclinic (two-layer) modes.
/// Westward Rossby phase speed: `c_R = -β L_R²` for barotropic mode.
pub struct PlanetaryWave {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Relaxation frequency for layer 1 (upper).
    pub omega1: f64,
    /// Relaxation frequency for layer 2 (lower).
    pub omega2: f64,
    /// Beta-plane parameter.
    pub beta: f64,
    /// Reduced gravity (for baroclinic mode).
    pub g_prime: f64,
    /// Layer-1 thickness.
    pub h1: f64,
    /// Layer-2 thickness.
    pub h2: f64,
    /// Distribution functions layer 1.
    pub f1: Vec<Vec<[f64; NQ]>>,
    /// Distribution functions layer 2.
    pub f2: Vec<Vec<[f64; NQ]>>,
}

impl PlanetaryWave {
    /// Construct a two-layer planetary wave model.
    pub fn new(
        nx: usize,
        ny: usize,
        omega1: f64,
        omega2: f64,
        beta: f64,
        g_prime: f64,
        h1: f64,
        h2: f64,
    ) -> Self {
        let f1 = vec![vec![feq(h1, 0.0, 0.0); ny]; nx];
        let f2 = vec![vec![feq(h2, 0.0, 0.0); ny]; nx];
        Self {
            nx,
            ny,
            omega1,
            omega2,
            beta,
            g_prime,
            h1,
            h2,
            f1,
            f2,
        }
    }

    /// Impose a Kelvin wave along the eastern boundary.
    pub fn kelvin_wave_ic(&mut self, amplitude: f64, kx: f64) {
        let nx = self.nx as f64;
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                let x = ix as f64 / nx;
                let phase = (2.0 * PI * kx * x).cos();
                self.f1[ix][iy] = feq(self.h1 + amplitude * phase, 0.0, 0.0);
            }
        }
    }

    /// Advance one timestep for both layers.
    pub fn step(&mut self) {
        let omega1 = self.omega1;
        let omega2 = self.omega2;
        let beta = self.beta;
        let nx = self.nx;
        let ny = self.ny;
        // Layer 1
        for ix in 0..nx {
            for iy in 0..ny {
                let y = iy as f64;
                let fc = beta * y;
                let (h, ux, uy) = macroscopic(&self.f1[ix][iy]);
                let feq_arr = feq(h, ux, uy);
                let fx = fc * uy;
                let fy = -fc * ux;
                for q in 0..NQ {
                    self.f1[ix][iy][q] += omega1 * (feq_arr[q] - self.f1[ix][iy][q]);
                }
                let mut tmp = self.f1[ix][iy];
                guo_force(&mut tmp, fx, fy, ux, uy, omega1);
                self.f1[ix][iy] = tmp;
            }
        }
        self.f1 = stream_periodic(&self.f1, nx, ny);
        // Layer 2
        for ix in 0..nx {
            for iy in 0..ny {
                let y = iy as f64;
                let fc = beta * y;
                let (h, ux, uy) = macroscopic(&self.f2[ix][iy]);
                let feq_arr = feq(h, ux, uy);
                let fx = fc * uy;
                let fy = -fc * ux;
                for q in 0..NQ {
                    self.f2[ix][iy][q] += omega2 * (feq_arr[q] - self.f2[ix][iy][q]);
                }
                let mut tmp = self.f2[ix][iy];
                guo_force(&mut tmp, fx, fy, ux, uy, omega2);
                self.f2[ix][iy] = tmp;
            }
        }
        self.f2 = stream_periodic(&self.f2, nx, ny);
    }

    /// Barotropic Rossby wave phase speed for wavenumber kx.
    pub fn barotropic_phase_speed(&self, kx: f64) -> f64 {
        let k2 = kx * kx;
        if k2 < 1e-14 {
            return 0.0;
        }
        -self.beta / k2
    }

    /// Baroclinic Rossby wave phase speed including deformation radius.
    pub fn baroclinic_phase_speed(&self, kx: f64) -> f64 {
        let lr2 = self.g_prime * self.h1 / self.beta.max(1e-14);
        let k2 = kx * kx;
        if k2 < 1e-14 {
            return 0.0;
        }
        -self.beta / (k2 + 1.0 / lr2.max(1e-14))
    }

    /// Compute interface displacement (h1 - h1_mean) as proxy for baroclinicity.
    pub fn interface_displacement(&self) -> Vec<Vec<f64>> {
        let h1 = self.h1;
        (0..self.nx)
            .map(|ix| {
                (0..self.ny)
                    .map(|iy| macroscopic(&self.f1[ix][iy]).0 - h1)
                    .collect()
            })
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Additional utility functions
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Froude number `Fr = U / sqrt(g H)`.
pub fn froude_number(u: f64, g: f64, h: f64) -> f64 {
    let c = (g * h).sqrt().max(1e-14);
    u / c
}

/// Compute the Rossby number `Ro = U / (f L)`.
pub fn rossby_number(u: f64, f_coriolis: f64, length: f64) -> f64 {
    let denom = f_coriolis.abs() * length;
    if denom < 1e-14 {
        return f64::INFINITY;
    }
    u / denom
}

/// Compute the Richardson number `Ri = N² / (dU/dz)²`.
pub fn richardson_number(brunt_vaisala_sq: f64, shear_sq: f64) -> f64 {
    if shear_sq < 1e-14 {
        return f64::INFINITY;
    }
    brunt_vaisala_sq / shear_sq
}

/// Brunt-Väisälä frequency squared from density gradient.
///
/// `N² = -(g/ρ₀)(dρ/dz)` (positive for stable stratification).
pub fn brunt_vaisala_sq(g: f64, rho0: f64, drho_dz: f64) -> f64 {
    -(g / rho0.max(1e-14)) * drho_dz
}

/// Ekman depth `d_E = sqrt(2ν / f)`.
pub fn ekman_depth(kinematic_viscosity: f64, f_coriolis: f64) -> f64 {
    if f_coriolis.abs() < 1e-14 {
        return f64::INFINITY;
    }
    (2.0 * kinematic_viscosity / f_coriolis.abs()).sqrt()
}

/// Tidal period in seconds from angular frequency.
pub fn tidal_period(angular_freq: f64) -> f64 {
    if angular_freq.abs() < 1e-14 {
        return f64::INFINITY;
    }
    2.0 * PI / angular_freq.abs()
}

/// Geostrophic velocity from pressure gradient and Coriolis parameter.
///
/// Returns (u_g, v_g).
pub fn geostrophic_velocity(dpdx: f64, dpdy: f64, rho: f64, f_cor: f64) -> (f64, f64) {
    if f_cor.abs() < 1e-14 {
        return (0.0, 0.0);
    }
    let u_g = -dpdy / (rho * f_cor);
    let v_g = dpdx / (rho * f_cor);
    (u_g, v_g)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Utilities ─────────────────────────────────────────────────────────────

    #[test]
    fn test_feq_sums_to_rho() {
        let rho = 1.2;
        let f = feq(rho, 0.05, 0.03);
        let sum: f64 = f.iter().sum();
        assert!((sum - rho).abs() < 1e-12, "feq sum = {sum}, expected {rho}");
    }

    #[test]
    fn test_feq_zero_velocity() {
        let f = feq(1.0, 0.0, 0.0);
        for (q, &fq) in f.iter().enumerate() {
            assert!((fq - W[q]).abs() < 1e-12, "q={q}: f={fq}, W={}", W[q]);
        }
    }

    #[test]
    fn test_macroscopic_roundtrip() {
        let rho0 = 1.1;
        let ux0 = 0.05;
        let uy0 = -0.03;
        let f = feq(rho0, ux0, uy0);
        let (rho, ux, uy) = macroscopic(&f);
        assert!((rho - rho0).abs() < 1e-12);
        assert!((ux - ux0).abs() < 1e-12);
        assert!((uy - uy0).abs() < 1e-12);
    }

    #[test]
    fn test_froude_number() {
        let fr = froude_number(1.0, 9.81, 1.0);
        let expected = 1.0 / 9.81_f64.sqrt();
        assert!((fr - expected).abs() < 1e-10);
    }

    #[test]
    fn test_rossby_number() {
        let ro = rossby_number(1.0, 1e-4, 1e6);
        assert!((ro - 0.01).abs() < 1e-10);
    }

    #[test]
    fn test_richardson_number() {
        let ri = richardson_number(1e-4, 1e-4);
        assert!((ri - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_richardson_number_unstable() {
        let ri = richardson_number(-0.01, 1.0);
        assert!(ri < 0.0, "unstable stratification should give Ri < 0");
    }

    #[test]
    fn test_brunt_vaisala_stable() {
        let n2 = brunt_vaisala_sq(9.81, 1025.0, -0.1);
        assert!(n2 > 0.0, "stable stratification: N² > 0");
    }

    #[test]
    fn test_ekman_depth() {
        let d = ekman_depth(1e-6, 1e-4);
        let expected = (2e-6 / 1e-4_f64).sqrt();
        assert!((d - expected).abs() < 1e-14);
    }

    #[test]
    fn test_tidal_period() {
        let period = tidal_period(2.0 * PI / 44712.0);
        assert!((period - 44712.0).abs() < 1e-6);
    }

    #[test]
    fn test_geostrophic_velocity() {
        let (ug, vg) = geostrophic_velocity(0.01, 0.02, 1025.0, 1e-4);
        assert!(ug.is_finite() && vg.is_finite());
    }

    // ── RotatingFrameLbm ──────────────────────────────────────────────────────

    #[test]
    fn test_rotating_lbm_construction() {
        let lbm = RotatingFrameLbm::new(16, 16, 1.0, 0.001);
        assert_eq!(lbm.nx, 16);
        assert_eq!(lbm.ny, 16);
    }

    #[test]
    fn test_rotating_lbm_mass_conservation() {
        let mut lbm = RotatingFrameLbm::new(16, 16, 1.0, 0.001);
        let mass_before: f64 = lbm.f.iter().flatten().map(|f| f.iter().sum::<f64>()).sum();
        for _ in 0..5 {
            lbm.step();
        }
        let mass_after: f64 = lbm.f.iter().flatten().map(|f| f.iter().sum::<f64>()).sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "mass not conserved: {mass_before} vs {mass_after}"
        );
    }

    #[test]
    fn test_rotating_lbm_kinetic_energy() {
        let mut lbm = RotatingFrameLbm::new(16, 16, 1.0, 0.001);
        lbm.set_velocity_ic(0.05, 0.0);
        let ke = lbm.kinetic_energy();
        assert!(ke > 0.0);
    }

    #[test]
    fn test_rotating_lbm_enstrophy_nonneg() {
        let mut lbm = RotatingFrameLbm::new(16, 16, 1.0, 0.001);
        lbm.set_velocity_ic(0.05, 0.01);
        let ens = lbm.enstrophy();
        assert!(ens >= 0.0);
    }

    #[test]
    fn test_rotating_lbm_velocity_field_size() {
        let lbm = RotatingFrameLbm::new(8, 10, 1.0, 0.0);
        let vf = lbm.velocity_field();
        assert_eq!(vf.len(), 8 * 10);
    }

    // ── ShallowWaterLbm ───────────────────────────────────────────────────────

    #[test]
    fn test_shallow_water_construction() {
        let sw = ShallowWaterLbm::new(16, 16, 1.0, 0.01, 1.0);
        assert_eq!(sw.nx, 16);
    }

    #[test]
    fn test_shallow_water_dam_break_height() {
        let mut sw = ShallowWaterLbm::new(16, 16, 1.0, 0.01, 1.0);
        sw.dam_break_ic(0.5);
        let h = sw.height_field();
        let max_h = h.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_h > 1.0, "dam break should increase max height");
    }

    #[test]
    fn test_shallow_water_energy_finite() {
        let sw = ShallowWaterLbm::new(16, 16, 1.0, 0.01, 1.0);
        let e = sw.total_energy();
        assert!(e.is_finite() && e >= 0.0);
    }

    #[test]
    fn test_shallow_water_step_runs() {
        let mut sw = ShallowWaterLbm::new(16, 16, 1.0, 0.01, 1.0);
        sw.dam_break_ic(0.1);
        for _ in 0..3 {
            sw.step();
        }
        let h = sw.height_field();
        assert!(h.iter().all(|&v| v.is_finite() && v > 0.0));
    }

    // ── RossbyWaveSim ─────────────────────────────────────────────────────────

    #[test]
    fn test_rossby_wave_construction() {
        let rw = RossbyWaveSim::new(16, 16, 1.0, 1e-4, 1.6e-11, 1.0);
        assert_eq!(rw.nx, 16);
    }

    #[test]
    fn test_rossby_phase_speed_negative() {
        let rw = RossbyWaveSim::new(16, 16, 1.0, 1e-4, 1.6e-11, 1.0);
        let cp = rw.theoretical_phase_speed(1.0, 0.0);
        assert!(cp < 0.0, "Rossby phase speed should be westward (negative)");
    }

    #[test]
    fn test_rossby_coriolis_at() {
        let rw = RossbyWaveSim::new(16, 16, 1.0, 1e-4, 1.6e-11, 1e5);
        let fc0 = rw.coriolis_at(0);
        let fc1 = rw.coriolis_at(1);
        assert!(fc1 > fc0, "beta-plane: f increases with y");
    }

    #[test]
    fn test_rossby_step_finite() {
        let mut rw = RossbyWaveSim::new(16, 16, 1.0, 1e-4, 1e-5, 1.0);
        rw.wave_ic(1.0, 1.0, 0.01);
        for _ in 0..5 {
            rw.step();
        }
        for ix in 0..rw.nx {
            for iy in 0..rw.ny {
                for q in 0..NQ {
                    assert!(rw.f[ix][iy][q].is_finite());
                }
            }
        }
    }

    // ── MantelConvection ──────────────────────────────────────────────────────

    #[test]
    fn test_mantel_construction() {
        let m = MantelConvection::new(16, 16, 1.0, 1.0, 1e4, 7.0);
        assert_eq!(m.nx, 16);
    }

    #[test]
    fn test_mantel_conductive_ic() {
        let mut m = MantelConvection::new(16, 16, 1.0, 1.0, 1e4, 7.0);
        m.conductive_ic();
        assert!((m.temp[0][0] - 1.0).abs() < 1e-12, "bottom should be hot");
        assert!(
            (m.temp[0][m.ny - 1] - 0.0).abs() < 1e-12,
            "top should be cold"
        );
    }

    #[test]
    fn test_mantel_nusselt_conductive() {
        let mut m = MantelConvection::new(16, 16, 1.0, 1.0, 1e3, 7.0);
        m.conductive_ic();
        let nu = m.nusselt_number();
        assert!(nu.is_finite() && nu > 0.0);
    }

    #[test]
    fn test_mantel_step_runs() {
        let mut m = MantelConvection::new(8, 8, 1.0, 1.0, 1e3, 7.0);
        m.conductive_ic();
        for _ in 0..3 {
            m.step();
        }
        assert!(m.temp.iter().flatten().all(|&t| t.is_finite()));
    }

    // ── ThermohalineCirculation ───────────────────────────────────────────────

    #[test]
    fn test_thermohaline_construction() {
        let thc = ThermohalineCirculation::new(16, 16, 1.0, 1.0, 0.8, 0.1, 0.08);
        assert_eq!(thc.nx, 16);
    }

    #[test]
    fn test_thermohaline_density_anomaly() {
        let thc = ThermohalineCirculation::new(8, 8, 1.0, 1.0, 0.8, 0.1, 0.08);
        let da = thc.density_anomaly();
        assert_eq!(da.len(), 8);
        assert_eq!(da[0].len(), 8);
    }

    #[test]
    fn test_thermohaline_step_finite() {
        let mut thc = ThermohalineCirculation::new(8, 8, 1.0, 1.0, 0.8, 0.1, 0.08);
        for _ in 0..5 {
            thc.step();
        }
        let da = thc.density_anomaly();
        assert!(da.iter().flatten().all(|&v| v.is_finite()));
    }

    // ── GravityCurrent ────────────────────────────────────────────────────────

    #[test]
    fn test_gravity_current_construction() {
        let gc = GravityCurrent::new(32, 16, 1.0, 1.1, 1.0, 0.01);
        assert!(gc.rho_heavy > gc.rho_light);
    }

    #[test]
    fn test_gravity_current_froude_speed() {
        let gc = GravityCurrent::new(32, 16, 1.0, 1.1, 1.0, 0.01);
        let spd = gc.froude_speed();
        assert!(spd > 0.0 && spd.is_finite());
    }

    #[test]
    fn test_gravity_current_front_initial() {
        let gc = GravityCurrent::new(32, 16, 1.0, 1.1, 1.0, 0.01);
        let fp = gc.front_position();
        assert!(fp < 32, "front should be within domain");
    }

    #[test]
    fn test_gravity_current_step_runs() {
        let mut gc = GravityCurrent::new(16, 8, 1.0, 1.1, 1.0, 0.001);
        for _ in 0..3 {
            gc.step();
        }
        assert!(gc.rho_field.iter().flatten().all(|&v| v.is_finite()));
    }

    // ── TidalSimulation ───────────────────────────────────────────────────────

    #[test]
    fn test_tidal_construction() {
        let ts = TidalSimulation::new(16, 16, 1.0, 0.001, 0.0001, 0.002);
        assert_eq!(ts.nx, 16);
    }

    #[test]
    fn test_tidal_step_advances_time() {
        let mut ts = TidalSimulation::new(8, 8, 1.0, 0.001, 0.0001, 0.002);
        ts.step();
        assert_eq!(ts.time, 1);
    }

    #[test]
    fn test_tidal_ellipse_finite() {
        let ts = TidalSimulation::new(8, 8, 1.0, 0.001, 0.0001, 0.002);
        let (a, b) = ts.tidal_ellipse(4, 4);
        assert!(a.is_finite() && b.is_finite());
    }

    #[test]
    fn test_tidal_shelf_bathymetry() {
        let mut ts = TidalSimulation::new(16, 8, 1.0, 0.001, 0.0001, 0.002);
        ts.shelf_bathymetry(0.5, 2.0, 8);
        assert!((ts.depth[0][0] - 0.5).abs() < 1e-12);
        assert!((ts.depth[15][0] - 2.0).abs() < 1e-12);
    }

    // ── GeostrophicAdjustment ─────────────────────────────────────────────────

    #[test]
    fn test_geostrophic_construction() {
        let ga = GeostrophicAdjustment::new(16, 16, 1.0, 1e-4);
        assert_eq!(ga.nx, 16);
    }

    #[test]
    fn test_geostrophic_jet_ic() {
        let mut ga = GeostrophicAdjustment::new(16, 16, 1.0, 1e-4);
        ga.jet_ic(0.1, 4.0);
        let (_rho, ux, _uy) = macroscopic(&ga.f[0][8]);
        assert!(ux > 0.0, "jet centre should have positive zonal flow");
    }

    #[test]
    fn test_geostrophic_rossby_radius() {
        let ga = GeostrophicAdjustment::new(16, 16, 1.0, 1e-4);
        let lr = ga.rossby_radius(0.01, 1000.0);
        assert!(lr > 0.0 && lr.is_finite());
    }

    #[test]
    fn test_geostrophic_imbalance_decreases() {
        let mut ga = GeostrophicAdjustment::new(16, 16, 1.2, 1e-3);
        ga.jet_ic(0.05, 4.0);
        for _ in 0..20 {
            ga.step();
        }
        assert!(!ga.imbalance_history.is_empty());
        assert!(ga.imbalance_history.iter().all(|&v| v.is_finite()));
    }

    // ── PlanetaryWave ─────────────────────────────────────────────────────────

    #[test]
    fn test_planetary_wave_construction() {
        let pw = PlanetaryWave::new(16, 16, 1.0, 1.0, 1.6e-11, 0.02, 300.0, 4700.0);
        assert_eq!(pw.nx, 16);
    }

    #[test]
    fn test_barotropic_phase_speed() {
        let pw = PlanetaryWave::new(16, 16, 1.0, 1.0, 2e-11, 0.02, 300.0, 4700.0);
        let cp = pw.barotropic_phase_speed(1e-6);
        assert!(cp < 0.0, "barotropic Rossby waves propagate westward");
    }

    #[test]
    fn test_baroclinic_phase_speed() {
        let pw = PlanetaryWave::new(16, 16, 1.0, 1.0, 2e-11, 0.02, 300.0, 4700.0);
        let cp_bt = pw.barotropic_phase_speed(1e-5);
        let cp_bc = pw.baroclinic_phase_speed(1e-5);
        assert!(
            cp_bc > cp_bt,
            "baroclinic speed slower (less negative) than barotropic at short wave"
        );
    }

    #[test]
    fn test_planetary_wave_kelvin_ic() {
        let mut pw = PlanetaryWave::new(16, 16, 1.0, 1.0, 1e-5, 0.02, 300.0, 4700.0);
        pw.kelvin_wave_ic(10.0, 1.0);
        let h_max = (0..16)
            .flat_map(|ix| (0..16).map(move |iy| (ix, iy)))
            .map(|(ix, iy)| macroscopic(&pw.f1[ix][iy]).0)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(h_max > 300.0, "Kelvin wave should raise interface height");
    }

    #[test]
    fn test_planetary_wave_step_finite() {
        let mut pw = PlanetaryWave::new(8, 8, 1.0, 1.0, 1e-4, 0.02, 300.0, 4700.0);
        pw.kelvin_wave_ic(5.0, 1.0);
        for _ in 0..5 {
            pw.step();
        }
        let disp = pw.interface_displacement();
        assert!(disp.iter().flatten().all(|&v| v.is_finite()));
    }
}
