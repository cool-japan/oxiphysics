// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shared double-distribution-function (DDF) core for the compressible D2Q9
//! solvers ([`super::CompressibleLbm`], [`super::EulerLbm`] and
//! [`super::NsCompressibleLbm`]).
//!
//! # Model
//!
//! * The density/momentum population `f` lives on the standard D2Q9 lattice
//!   with the lattice (isothermal) sound speed `cs² = 1/3`, which is the only
//!   value for which the D2Q9 weights reproduce the Hermite moments.  Its
//!   equilibrium is the Hermite expansion of order 2 or 3 (see
//!   [`super::MaCorrection`]), so `Σ f_eq = ρ` and `Σ f_eq c = ρu` exactly.
//! * The total-energy population `g` carries `ρE = ρ(|u|²/2 + T/(γ−1))`
//!   (ideal gas with `R = 1` in lattice units, `p = ρT`).  Its equilibrium is
//!   `g_eq_i = f_eq_i E + w_i p (c_i·u)/cs²`, giving `Σ g_eq = ρE` and
//!   `Σ g_eq c = (ρE + p) u`, i.e. the Euler energy flux.
//! * The thermodynamic pressure `ρT` is coupled back into the momentum
//!   equation through a Guo forcing term `F = −∇(ρT − ρ cs²)` evaluated with
//!   periodic central differences.  The lattice already supplies `−∇(ρ cs²)`,
//!   so the macroscopic momentum equation carries the full `−∇(ρT)`.  On a
//!   periodic grid the discrete forces telescope, so total momentum is
//!   conserved to round-off.
//!
//! Mass (`Σf`) and total energy (`Σg`) are conserved exactly by both the
//! collision and the periodic streaming, because the temperature is always
//! recovered from `Σg` with the same velocity that is used to build `g_eq`.

use super::{EX9, EY9, LATTICE_CS2, MaCorrection, W9};

/// D2Q9 population field indexed `[y][x][q]`.
pub(super) type Populations = Vec<Vec<Vec<f64>>>;

/// Relaxation and thermodynamic parameters of one DDF collision.
#[derive(Debug, Clone, Copy)]
pub(super) struct DdfParams {
    /// Relaxation frequency of the density/momentum population `f`.
    pub omega_f: f64,
    /// Relaxation frequency of the total-energy population `g`.
    pub omega_g: f64,
    /// Adiabatic exponent γ.
    pub gamma: f64,
    /// Hermite order of the `f` equilibrium (2 or ≥3).
    pub ma_order: usize,
}

/// Mutable view of all fields of a DDF solver.
pub(super) struct DdfState<'a> {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Density/momentum populations.
    pub f: &'a mut Populations,
    /// Total-energy populations.
    pub g: &'a mut Populations,
    /// Density field.
    pub rho: &'a mut Vec<Vec<f64>>,
    /// Velocity field.
    pub u: &'a mut Vec<Vec<[f64; 2]>>,
    /// Temperature field (`p = ρT`).
    pub temp: &'a mut Vec<Vec<f64>>,
    /// Guo pressure-correction force field.
    pub force: &'a mut Vec<Vec<[f64; 2]>>,
}

/// Total-energy equilibrium for the DDF model.
///
/// `g_eq_i = f_eq_i E + w_i p (c_i·u)/cs²` with `E = |u|²/2 + T/(γ−1)` and
/// `p = ρT`.  For any `f_eq` with `Σ f_eq = ρ` and `Σ f_eq c = ρu` this gives
/// `Σ g_eq = ρE` and `Σ g_eq c = (ρE + p) u`.
pub fn energy_equilibrium(
    feq: &[f64; 9],
    rho: f64,
    ux: f64,
    uy: f64,
    temp: f64,
    gamma: f64,
) -> [f64; 9] {
    let e_tot = 0.5 * (ux * ux + uy * uy) + temp / (gamma - 1.0);
    let p = rho * temp;
    let mut geq = [0.0f64; 9];
    for (q, g) in geq.iter_mut().enumerate() {
        let cu = f64::from(EX9[q]) * ux + f64::from(EY9[q]) * uy;
        *g = feq[q] * e_tot + W9[q] * p * cu / LATTICE_CS2;
    }
    geq
}

/// Temperature recovered from the total energy density.
///
/// `T = (γ−1)(ρE/ρ − |u|²/2)`.
pub(super) fn temperature_from_energy(rho: f64, energy: f64, u: [f64; 2], gamma: f64) -> f64 {
    (gamma - 1.0) * (energy / rho - 0.5 * (u[0] * u[0] + u[1] * u[1]))
}

/// Guo force `F = −∇(ρT − ρ cs²)` with periodic second-order central
/// differences.  Returns a field indexed `[y][x]`.
pub(super) fn pressure_correction_force(
    rho: &[Vec<f64>],
    temp: &[Vec<f64>],
    nx: usize,
    ny: usize,
) -> Vec<Vec<[f64; 2]>> {
    let excess = |x: usize, y: usize| -> f64 { rho[y][x] * (temp[y][x] - LATTICE_CS2) };
    let mut force = vec![vec![[0.0f64; 2]; nx]; ny];
    for (y, row) in force.iter_mut().enumerate() {
        let yp = (y + 1) % ny;
        let ym = (y + ny - 1) % ny;
        for (x, fxy) in row.iter_mut().enumerate() {
            let xp = (x + 1) % nx;
            let xm = (x + nx - 1) % nx;
            fxy[0] = -0.5 * (excess(xp, y) - excess(xm, y));
            fxy[1] = -0.5 * (excess(x, yp) - excess(x, ym));
        }
    }
    force
}

/// Guo source term `S_i = (1 − ω/2) w_i [(c_i − u)/cs² + (c_i·u) c_i/cs⁴]·F`.
fn guo_source(q: usize, u: [f64; 2], force: [f64; 2], omega: f64) -> f64 {
    let cx = f64::from(EX9[q]);
    let cy = f64::from(EY9[q]);
    let cu = cx * u[0] + cy * u[1];
    let inv_cs2 = 1.0 / LATTICE_CS2;
    let inv_cs4 = inv_cs2 * inv_cs2;
    let term_x = (cx - u[0]) * inv_cs2 + cu * cx * inv_cs4;
    let term_y = (cy - u[1]) * inv_cs2 + cu * cy * inv_cs4;
    (1.0 - 0.5 * omega) * W9[q] * (term_x * force[0] + term_y * force[1])
}

/// Build `f` and `g` at equilibrium from the current `rho`, `u` and `temp`
/// fields, and compute the matching pressure-correction force.
///
/// `f` is built with the velocity `u − F/(2ρ)` so that the Guo-corrected
/// velocity recovered by [`update_macros`] equals the requested `u`.
pub(super) fn initialize(state: &mut DdfState<'_>, params: DdfParams) {
    *state.force = pressure_correction_force(state.rho, state.temp, state.nx, state.ny);
    let corr = MaCorrection::new(params.ma_order, LATTICE_CS2.sqrt());
    for y in 0..state.ny {
        for x in 0..state.nx {
            let rho = state.rho[y][x];
            let [ux, uy] = state.u[y][x];
            let [fx, fy] = state.force[y][x];
            let (sx, sy) = if rho > 0.0 {
                (ux - 0.5 * fx / rho, uy - 0.5 * fy / rho)
            } else {
                (ux, uy)
            };
            let f_shift = corr.corrected_feq(rho, sx, sy);
            let feq = corr.corrected_feq(rho, ux, uy);
            let geq = energy_equilibrium(&feq, rho, ux, uy, state.temp[y][x], params.gamma);
            state.f[y][x] = f_shift.to_vec();
            state.g[y][x] = geq.to_vec();
        }
    }
}

/// BGK collision of `f` (with Guo forcing) and `g`.
///
/// Uses the force stored in `state.force`, which must be the one used to
/// define the current velocity field (see [`update_macros`]).
pub(super) fn collide(state: &mut DdfState<'_>, params: DdfParams) {
    let corr = MaCorrection::new(params.ma_order, LATTICE_CS2.sqrt());
    let omega_f = params.omega_f;
    let omega_g = params.omega_g;
    for y in 0..state.ny {
        for x in 0..state.nx {
            let rho = state.rho[y][x];
            let u = state.u[y][x];
            let force = state.force[y][x];
            let feq = corr.corrected_feq(rho, u[0], u[1]);
            let geq = energy_equilibrium(&feq, rho, u[0], u[1], state.temp[y][x], params.gamma);
            let fc = &mut state.f[y][x];
            for (q, fq) in fc.iter_mut().enumerate() {
                *fq = *fq * (1.0 - omega_f) + omega_f * feq[q] + guo_source(q, u, force, omega_f);
            }
            let gc = &mut state.g[y][x];
            for (q, gq) in gc.iter_mut().enumerate() {
                *gq = *gq * (1.0 - omega_g) + omega_g * geq[q];
            }
        }
    }
}

/// Periodic D2Q9 streaming of one population field.
pub(super) fn stream_periodic(pops: &mut Populations, nx: usize, ny: usize) {
    let mut streamed = vec![vec![vec![0.0f64; 9]; nx]; ny];
    for (y, row) in pops.iter().enumerate() {
        for (x, cell) in row.iter().enumerate() {
            for (q, &value) in cell.iter().enumerate() {
                let xn = (x as i64 + i64::from(EX9[q])).rem_euclid(nx as i64) as usize;
                let yn = (y as i64 + i64::from(EY9[q])).rem_euclid(ny as i64) as usize;
                streamed[yn][xn][q] = value;
            }
        }
    }
    *pops = streamed;
}

/// Recover `ρ`, `u`, `T` and the pressure-correction force from `f` and `g`.
///
/// The velocity depends on the force (Guo: `ρu = Σ f c + F/2`), the force
/// depends on the temperature and the temperature depends on the velocity.
/// This is resolved with one predictor–corrector pass:
/// `u* = Σfc/ρ → T*(u*) → F(ρ, T*) → u = (Σfc + F/2)/ρ → T(u)`.
/// The final temperature uses the final velocity, so `Σ g = ρE(u, T)` holds
/// exactly and the next collision conserves energy.
pub(super) fn update_macros(state: &mut DdfState<'_>, gamma: f64) {
    let nx = state.nx;
    let ny = state.ny;
    // Raw moments: [Σf c_x, Σf c_y, Σg].
    let mut raw = vec![vec![[0.0f64; 3]; nx]; ny];
    for (y, raw_row) in raw.iter_mut().enumerate() {
        for (x, raw_cell) in raw_row.iter_mut().enumerate() {
            let mut rho = 0.0;
            let mut mx = 0.0;
            let mut my = 0.0;
            for (q, &fq) in state.f[y][x].iter().enumerate() {
                rho += fq;
                mx += f64::from(EX9[q]) * fq;
                my += f64::from(EY9[q]) * fq;
            }
            let energy: f64 = state.g[y][x].iter().sum();
            *raw_cell = [mx, my, energy];
            state.rho[y][x] = rho;
            if rho > 0.0 {
                let u_star = [mx / rho, my / rho];
                state.temp[y][x] = temperature_from_energy(rho, energy, u_star, gamma);
            }
        }
    }
    *state.force = pressure_correction_force(state.rho, state.temp, nx, ny);
    for (y, raw_row) in raw.iter().enumerate() {
        for (x, &[mx, my, energy]) in raw_row.iter().enumerate() {
            let rho = state.rho[y][x];
            if rho > 0.0 {
                let [fx, fy] = state.force[y][x];
                let u = [(mx + 0.5 * fx) / rho, (my + 0.5 * fy) / rho];
                state.u[y][x] = u;
                state.temp[y][x] = temperature_from_energy(rho, energy, u, gamma);
            }
        }
    }
}

/// Conserved totals `(Σρ, Σρu, Σρe_total)` of a DDF state, where the energy
/// total is `Σ g`.
pub(super) fn totals(
    f: &Populations,
    g: &Populations,
    rho: &[Vec<f64>],
    u: &[Vec<[f64; 2]>],
) -> (f64, [f64; 2], f64) {
    let mass: f64 = f.iter().flatten().flatten().sum();
    let energy: f64 = g.iter().flatten().flatten().sum();
    let mut momentum = [0.0f64; 2];
    for (rho_row, u_row) in rho.iter().zip(u.iter()) {
        for (&r, v) in rho_row.iter().zip(u_row.iter()) {
            momentum[0] += r * v[0];
            momentum[1] += r * v[1];
        }
    }
    (mass, momentum, energy)
}
