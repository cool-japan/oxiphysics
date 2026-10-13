// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neutral atmospheric boundary layer with a log-law wall model.

use super::{CS2, MIRROR_Y, NQ, clamp, feq, guo_force, macroscopic, stream_walled_y};

/// von Kármán constant.
const VON_KARMAN: f64 = 0.41;

/// Height of the first cell centre above the wall plane (the wall plane lies
/// half-way between the bottom row and the reflected ghost row).
const FIRST_CELL_HEIGHT: f64 = 0.5;

/// Neutral atmospheric boundary layer with a rough-wall log-law stress model.
///
/// A Smagorinsky large-eddy closure sets the local viscosity.  There is no
/// thermal stratification, so this is a *neutral* boundary layer (no
/// Monin–Obukhov stability correction).  The vertical coordinate is `y`; the
/// domain is periodic in x with free-slip (specular) walls at the bottom and
/// top, located half a lattice spacing beyond the first and last row.
///
/// * **Wall model** — the bottom row feels the rough-wall stress of the
///   logarithmic law, `τ_w = ρ u*²` with `u* = κ |u| / ln(1 + z₁/z₀)` (`z₁ = ½`
///   the height of the first cell centre, `z₀` the roughness length), as the
///   body force density `F_x = −τ_w sign(u_x) / Δz` (`Δz = 1`).  A smaller
///   roughness length gives a smaller drag.
/// * **Geostrophic forcing** — the top row is relaxed towards the geostrophic
///   wind by the body force density `F_x = r ρ (u_geo − u_x)` with the
///   relaxation rate `r =` [`relax_rate`](Self::relax_rate).
///
/// Both forces are added with the Guo scheme, so they change the momentum and
/// leave the mass untouched.
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
    /// Relaxation rate `r` of the top row towards the geostrophic wind
    /// (per time step); `0.01` unless changed.
    pub relax_rate: f64,
    /// Distribution functions.
    pub f: Vec<Vec<[f64; NQ]>>,
    /// Turbulent viscosity field.
    pub nu_t: Vec<Vec<f64>>,
}

impl AtmosphericBoundaryLayer {
    /// Construct a new ABL simulation domain with the flow initialised to the
    /// geostrophic wind.
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
            relax_rate: 0.01,
            f,
            nu_t,
        }
    }

    /// Velocity `(u_x, u_y)` at row `iy`; rows beyond a free-slip wall mirror
    /// the adjacent fluid row (`u_x` kept, `u_y` reversed).
    fn velocity_with_mirror(&self, ix: usize, iy: isize) -> (f64, f64) {
        if iy < 0 {
            let (_, ux, uy) = macroscopic(&self.f[ix][0]);
            (ux, -uy)
        } else if iy as usize >= self.ny {
            let (_, ux, uy) = macroscopic(&self.f[ix][self.ny - 1]);
            (ux, -uy)
        } else {
            let (_, ux, uy) = macroscopic(&self.f[ix][iy as usize]);
            (ux, uy)
        }
    }

    /// Compute local strain-rate magnitude |S| (central differences, periodic
    /// in x, mirrored at the free-slip walls).
    fn strain_rate(&self, ix: usize, iy: usize) -> f64 {
        let nx = self.nx;
        let ixp = (ix + 1) % nx;
        let ixm = (ix + nx - 1) % nx;
        let row = iy as isize;
        let (uxpx, uypx) = self.velocity_with_mirror(ixp, row);
        let (uxmx, uymx) = self.velocity_with_mirror(ixm, row);
        let (uxpy, uypy) = self.velocity_with_mirror(ix, row + 1);
        let (uxmy, uymy) = self.velocity_with_mirror(ix, row - 1);
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

    /// Friction velocity `u* = κ |u| / ln(1 + z₁/z₀)` of the logarithmic law
    /// for the first-cell velocity `u`.
    pub fn friction_velocity_for(&self, u: f64) -> f64 {
        let z0 = self.z0.max(1e-12);
        VON_KARMAN * u.abs() / (1.0 + FIRST_CELL_HEIGHT / z0).ln().max(1e-6)
    }

    /// Mean friction velocity over the bottom row of the current state.
    pub fn friction_velocity(&self) -> f64 {
        if self.nx == 0 || self.ny == 0 {
            return 0.0;
        }
        let sum: f64 = (0..self.nx)
            .map(|ix| self.friction_velocity_for(macroscopic(&self.f[ix][0]).1))
            .sum();
        sum / self.nx as f64
    }

    /// Total mass `Σ ρ` over the grid (conserved to rounding).
    pub fn total_mass(&self) -> f64 {
        self.f
            .iter()
            .flatten()
            .map(|cell| cell.iter().sum::<f64>())
            .sum()
    }

    /// Advance one timestep with Smagorinsky LES, the log-law wall stress and
    /// the geostrophic relaxation, then stream with free-slip walls.
    pub fn step(&mut self) {
        self.update_turbulent_viscosity();
        let (nx, ny) = (self.nx, self.ny);
        let omega_base = self.omega_base;
        let u_geo = self.u_geo;
        let relax_rate = self.relax_rate;
        for ix in 0..nx {
            for iy in 0..ny {
                let nu_t_local = self.nu_t[ix][iy];
                let nu_eff = (1.0 / omega_base - 0.5) * CS2 + nu_t_local;
                let omega_eff = clamp(1.0 / (nu_eff / CS2 + 0.5), 0.01, 1.99);
                let (rho, ux0, uy0) = macroscopic(&self.f[ix][iy]);
                // Body force density along x: geostrophic relaxation on the top
                // row, log-law wall stress on the bottom row.
                let mut fx = 0.0;
                if iy + 1 == ny {
                    fx += relax_rate * rho * (u_geo - ux0);
                }
                if iy == 0 {
                    let u_star = self.friction_velocity_for(ux0);
                    fx -= rho * u_star * u_star * ux0.signum();
                }
                let ux = if rho > 1e-12 {
                    ux0 + 0.5 * fx / rho
                } else {
                    0.0
                };
                let feq_arr = feq(rho, ux, uy0);
                let cell = &mut self.f[ix][iy];
                for (f_q, feq_q) in cell.iter_mut().zip(feq_arr.iter()) {
                    *f_q += omega_eff * (*feq_q - *f_q);
                }
                guo_force(cell, fx, 0.0, ux, uy0, omega_eff);
            }
        }
        self.f = stream_walled_y(&self.f, nx, ny, &MIRROR_Y);
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
