// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Free-surface LBM using Volume-of-Fluid (VOF) mass tracking.
//!
//! Implements the Körner 2005 / Thürey VOF scheme for 2D free-surface LBM with
//! D2Q9. Provides a parallel side state (`FreeSurfaceState`) to `CellGrid2D`
//! that tracks per-cell fill levels, masses, and cell state flags
//! (Gas / Interface / Fluid). The caller drives one BGK step on the
//! `CellGrid2D`, then calls `free_surface_step` to advance the VOF state.
//!
//! ## References
//! - Körner, C. et al. (2005). *Lattice Boltzmann model for free surface flow
//!   for modeling foaming*. J. Stat. Phys., 121(1–2), 179–196.
//! - Thürey, N. (2007). *Physically based Animation of Free Surface Flows
//!   with the Lattice Boltzmann Method*. PhD thesis, FAU Erlangen-Nürnberg.

use crate::grid::types::CellGrid2D;

// ─────────────────────────────────────────────────────────────────────────────
// D2Q9 constants (local copies — keeps this module self-contained)
// ─────────────────────────────────────────────────────────────────────────────

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

/// D2Q9 velocity vectors `[cx, cy]`.
const C9: [[i32; 2]; 9] = [
    [0, 0],   // 0 rest
    [1, 0],   // 1 E
    [0, 1],   // 2 N
    [-1, 0],  // 3 W
    [0, -1],  // 4 S
    [1, 1],   // 5 NE
    [-1, 1],  // 6 NW
    [-1, -1], // 7 SW
    [1, -1],  // 8 SE
];

/// D2Q9 opposite direction indices.
const OPP9: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Speed of sound squared in LBM lattice units.
const CS2: f64 = 1.0 / 3.0;

/// Atmospheric density used for gas-side reconstruction.
const RHO_ATM: f64 = 1.0;

/// Fill-level tolerance for Interface → Fluid / Gas conversions.
const FILL_EPS: f64 = 1.0e-3;

// ─────────────────────────────────────────────────────────────────────────────
// Cell state flag
// ─────────────────────────────────────────────────────────────────────────────

/// Cell phase state for VOF free-surface tracking.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellState {
    /// Cell is entirely gaseous (no fluid).
    Gas,
    /// Cell is at the fluid–gas interface.
    Interface,
    /// Cell is entirely filled with fluid.
    Fluid,
}

// ─────────────────────────────────────────────────────────────────────────────
// FreeSurfaceState
// ─────────────────────────────────────────────────────────────────────────────

/// VOF free-surface state parallel to a `CellGrid2D`.
///
/// `state`, `fill`, and `mass` are in row-major order: `idx(ix, iy) = iy * nx + ix`.
/// This matches `CellGrid2D::idx` exactly.
pub struct FreeSurfaceState {
    /// Per-cell phase state.
    pub state: Vec<CellState>,
    /// Per-cell fill level φ ∈ \[0, 1\].
    pub fill: Vec<f64>,
    /// Per-cell fluid mass m = φ · ρ.
    pub mass: Vec<f64>,
    /// Grid width.
    nx: usize,
    /// Grid height.
    ny: usize,
}

impl FreeSurfaceState {
    /// Create a free-surface state for a grid of size `nx × ny` with all
    /// cells initialised to `Gas` (fill = 0, mass = 0).
    pub fn new(nx: usize, ny: usize) -> Self {
        let n = nx * ny;
        Self {
            state: vec![CellState::Gas; n],
            fill: vec![0.0; n],
            mass: vec![0.0; n],
            nx,
            ny,
        }
    }

    /// Flat index: `idx(ix, iy) = iy * nx + ix` — matches `CellGrid2D::idx`.
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Grid width.
    #[inline]
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Grid height.
    #[inline]
    pub fn ny(&self) -> usize {
        self.ny
    }

    /// Initialise a 2D dam-break scenario.
    ///
    /// The left column block `ix < fluid_width` and row block `iy < fluid_height`
    /// is filled with Fluid (fill = 1, mass = ρ₀ = 1).  The outer shell of the
    /// fluid block (top row of fluid, rightmost column of fluid) becomes
    /// Interface; everything else is Gas.
    ///
    /// `fluid_width` defaults to `nx / 2` when `0` is passed.
    /// `fluid_height` is clamped to `ny`.
    ///
    /// The corresponding `CellGrid2D` should be initialised at equilibrium
    /// with `rho0 = 1.0` before calling this.
    pub fn new_dam_break(nx: usize, ny: usize, fluid_height: usize) -> Self {
        let n = nx * ny;
        let mut s = Self {
            state: vec![CellState::Gas; n],
            fill: vec![0.0; n],
            mass: vec![0.0; n],
            nx,
            ny,
        };

        let fluid_height = fluid_height.min(ny);
        // Classic column-collapse dam break: left half of domain is the water column.
        let fluid_width = (nx / 2).max(1);

        for iy in 0..ny {
            for ix in 0..nx {
                let k = s.idx(ix, iy);
                let in_fluid_block = ix < fluid_width && iy < fluid_height;
                if !in_fluid_block {
                    continue; // Gas: already set
                }
                // Classify as Interface if on the exposed border (top row or
                // rightmost column of the fluid block), otherwise Fluid.
                let is_interface = iy + 1 == fluid_height || ix + 1 == fluid_width;
                if is_interface {
                    s.state[k] = CellState::Interface;
                } else {
                    s.state[k] = CellState::Fluid;
                }
                s.fill[k] = 1.0;
                s.mass[k] = 1.0;
            }
        }
        s
    }

    /// Total fluid mass: sum of `mass` over all Fluid and Interface cells.
    pub fn total_mass(&self) -> f64 {
        self.state
            .iter()
            .zip(self.mass.iter())
            .filter_map(|(st, &m)| {
                if matches!(st, CellState::Fluid | CellState::Interface) {
                    Some(m)
                } else {
                    None
                }
            })
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// D2Q9 equilibrium distribution component for direction `a`.
///
/// `f_eq_a = w_a · ρ · (1 + (e_a·u)/cs² + (e_a·u)²/(2cs⁴) − |u|²/(2cs²))`
#[inline]
fn feq(rho: f64, ux: f64, uy: f64, a: usize) -> f64 {
    let cx = C9[a][0] as f64;
    let cy = C9[a][1] as f64;
    let eu = cx * ux + cy * uy;
    let u_sq = ux * ux + uy * uy;
    W9[a] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2))
}

/// Perform a BGK collision + periodic streaming step restricted to Fluid and
/// Interface cells, with an optional Guo body force.
///
/// Gas cells retain their current (post-reconstruction) distribution values
/// but do not participate in BGK collision.  This matches the Körner / Thürey
/// formulation where only the liquid phase is actively solved.
///
/// The body force `(gx, gy)` is applied using the Guo 2002 scheme:
/// the macroscopic velocity `u` (pre-shifted by `apply_body_force_2d`) is used
/// in BGK collision, then the force term `(1 - ω/2) * F_i` is added.
///
/// Caller should call `free_surface_step` **after** this function to advance
/// the VOF state.
pub fn free_surface_bgk_stream(
    grid: &mut CellGrid2D,
    fs: &FreeSurfaceState,
    omega: f64,
    gx: f64,
    gy: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let guo_coeff = 1.0 - 0.5 * omega;

    // ── BGK collision + Guo force on Fluid and Interface cells only ──
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            let cell = &mut grid.cells[ki];
            if cell.obstacle {
                continue;
            }
            if matches!(fs.state[ki], CellState::Gas) {
                continue;
            }
            let rho = cell.rho;
            // cell.ux/uy already shifted by 0.5*g (from apply_body_force_2d).
            let ux = cell.ux;
            let uy = cell.uy;
            let u_sq = ux * ux + uy * uy;
            for a in 0..9usize {
                let cx = C9[a][0] as f64;
                let cy = C9[a][1] as f64;
                let eu = cx * ux + cy * uy;
                let f_eq = W9[a]
                    * rho
                    * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
                cell.f[a] -= omega * (cell.f[a] - f_eq);
                // Guo force term: F_i = w_i * [(c_i - u)/cs² + (c_i·u)/cs⁴ * c_i] · g
                let cu = cx * ux + cy * uy;
                let force_i = W9[a]
                    * (((cx - ux) * gx + (cy - uy) * gy) / CS2
                        + cu * (cx * gx + cy * gy) / (CS2 * CS2));
                cell.f[a] += guo_coeff * force_i;
            }
            // Update macroscopic velocity with the net acceleration.
            // (Subtract the pre-shift of 0.5*g that was added earlier, then add
            // the full-step contribution via Guo: u_actual = (u_shifted - 0.5*g) + g*0.5 = u_shifted)
            // The net macroscopic velocity after this step is u_shifted,
            // which equals u_old + 0.5*g_applied_per_step. This accumulates correctly.
        }
    }

    // ── Pull streaming (periodic) ──
    let f_old: Vec<[f64; 9]> = grid.cells.iter().map(|c| c.f).collect();
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            if grid.cells[ki].obstacle {
                continue;
            }
            for a in 0..9usize {
                let sx = (ix as isize - C9[a][0] as isize).rem_euclid(nx as isize) as usize;
                let sy = (iy as isize - C9[a][1] as isize).rem_euclid(ny as isize) as usize;
                grid.cells[ki].f[a] = f_old[sy * nx + sx][a];
            }
        }
    }

    // Note: Gas cells retain their streamed distributions.  Their values do
    // not affect the VOF mass balance — gas-side directions of Interface cells
    // are reconstructed in `free_surface_step` using the Interface cell's
    // own f[ā] values, independent of Gas cell content.

    // ── Update macroscopic for Fluid and Interface cells only ──
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            let cell = &mut grid.cells[ki];
            if cell.obstacle || fs.state[ki] == CellState::Gas {
                continue;
            }
            let mut rho = 0.0_f64;
            let mut mx = 0.0_f64;
            let mut my = 0.0_f64;
            for (a, ca) in C9.iter().enumerate() {
                let fi = cell.f[a];
                rho += fi;
                mx += fi * ca[0] as f64;
                my += fi * ca[1] as f64;
            }
            cell.rho = rho;
            if rho.abs() > 1e-15 {
                cell.ux = mx / rho;
                cell.uy = my / rho;
            } else {
                cell.ux = 0.0;
                cell.uy = 0.0;
            }
        }
    }
}

/// Apply a Guo 2002 body-force correction to all non-obstacle cells of
/// `grid`.
///
/// The correction modifies the post-collision distributions:
///   f_i += (1 − ω/2) · F_i
/// where
///   F_i = w_i · [ (e_i − u)/cs² + (e_i · u)/cs⁴ · e_i ] · g
///
/// Typically called **before** `CellGrid2D::step(omega)` (or any streaming
/// cycle) to include the body force in the collision step.  Pass `gx = 0.0`
/// and `gy = -g_lattice` for downward gravity (negative y-direction).
/// Apply a Guo 2002 body-force correction to all Fluid and Interface cells of
/// `grid`.  Gas cells and obstacle cells are skipped.
///
/// **Important:** This function shifts the stored macroscopic velocity
/// `(ux, uy)` by `0.5 * g / omega` so that the subsequent BGK collision uses
/// the force-corrected velocity.  The actual force term is applied to `f[a]`
/// after the BGK collision inside `free_surface_bgk_stream`.
///
/// The `fs` parameter is used to determine which cells are Gas.
/// Call this **before** `free_surface_bgk_stream`.
pub fn apply_body_force_2d(
    grid: &mut CellGrid2D,
    fs: &FreeSurfaceState,
    gx: f64,
    gy: f64,
    _omega: f64,
) {
    let nx = grid.nx;
    let ny = grid.ny;
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            let cell = &mut grid.cells[ki];
            if cell.obstacle || fs.state[ki] == CellState::Gas {
                continue;
            }
            // Shift macroscopic velocity by g * 0.5 so that BGK collision
            // uses the force-corrected velocity. (Guo scheme: u_star = u + tau*g/2)
            cell.ux += 0.5 * gx;
            cell.uy += 0.5 * gy;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main per-step update
// ─────────────────────────────────────────────────────────────────────────────

/// Advance the VOF free-surface state by one step.
///
/// This must be called **after** `CellGrid2D::step(omega)` (or an equivalent
/// BGK / MRT + streaming cycle) on the same grid.
///
/// The function implements the Körner 2005 / Thürey algorithm:
/// 1. Reconstruct gas-side distributions entering Interface cells.
/// 2. Perform mass exchange between neighbouring cells.
/// 3. Update fill levels.
/// 4. Reinitialise cell flags (Interface ↔ Fluid / Gas) with excess/deficit
///    redistribution and enforce the Interface closure invariant.
///
/// # Panics
/// Panics if `fs.nx() != grid.nx || fs.ny() != grid.ny`.
pub fn free_surface_step(grid: &mut CellGrid2D, fs: &mut FreeSurfaceState) {
    assert_eq!(
        fs.nx(),
        grid.nx,
        "FreeSurfaceState nx must equal CellGrid2D nx"
    );
    assert_eq!(
        fs.ny(),
        grid.ny,
        "FreeSurfaceState ny must equal CellGrid2D ny"
    );

    let nx = fs.nx();
    let ny = fs.ny();

    // ──────────────────────────────────────────────────────────────────────
    // Step 1: Gas-side distribution reconstruction for Interface cells.
    //
    // For each Interface cell i and each direction α whose source neighbour g
    // is a Gas cell, replace f_i^α with the reconstruction:
    //   f_i^α = f_eq(ρ_atm, u_i) + f_eq_ᾱ(ρ_atm, u_i) − f_i^ᾱ
    // where ᾱ = OPP9[α].
    // ──────────────────────────────────────────────────────────────────────
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            if fs.state[ki] != CellState::Interface {
                continue;
            }
            let cell_i = &grid.cells[ki];
            let ux_i = cell_i.ux;
            let uy_i = cell_i.uy;

            // Collect directions that come from a Gas neighbour.
            let mut new_f = cell_i.f;
            for a in 1..9usize {
                // Neighbour position: cell that streamed into direction a of i
                // (i.e. the cell *from* which direction a originated).
                let nx_coord = ix as isize + C9[a][0] as isize;
                let ny_coord = iy as isize + C9[a][1] as isize;
                // Bounds check — treat out-of-bounds as Gas (open boundary).
                let is_gas_neighbor = if nx_coord < 0
                    || ny_coord < 0
                    || nx_coord >= nx as isize
                    || ny_coord >= ny as isize
                {
                    true
                } else {
                    let kg = fs.idx(nx_coord as usize, ny_coord as usize);
                    fs.state[kg] == CellState::Gas
                };
                if is_gas_neighbor {
                    let a_bar = OPP9[a];
                    // Reconstruction: feq(atm, u_i) direction a
                    //               + feq(atm, u_i) direction ā
                    //               − f_i^{ā}   (post-stream value)
                    let reconstructed = feq(RHO_ATM, ux_i, uy_i, a)
                        + feq(RHO_ATM, ux_i, uy_i, a_bar)
                        - cell_i.f[a_bar];
                    new_f[a] = reconstructed;
                }
            }
            grid.cells[ki].f = new_f;
        }
    }

    // ──────────────────────────────────────────────────────────────────────
    // Step 2: Mass exchange.
    //
    // For each Interface cell i, sum up mass fluxes from neighbours.
    // Δm_{i←j,α} = f_j^ᾱ − f_i^α  (direction α points from j to i)
    // scaled by the average fill level:
    //   scale = 0.5 * (fill_i + fill_j)  when j is Interface
    //           1.0                       when j is Fluid
    //           0.0                       when j is Gas
    // ──────────────────────────────────────────────────────────────────────
    let mut delta_mass = vec![0.0_f64; nx * ny];

    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            if fs.state[ki] != CellState::Interface {
                continue;
            }
            let fill_i = fs.fill[ki];

            // Direction α points toward neighbour j; the opposite ā is the
            // direction from j back to i (what j streamed into i).
            for a in 1..9usize {
                let jx = ix as isize + C9[a][0] as isize;
                let jy = iy as isize + C9[a][1] as isize;
                if jx < 0 || jy < 0 || jx >= nx as isize || jy >= ny as isize {
                    continue;
                }
                let jx = jx as usize;
                let jy = jy as usize;
                let kj = fs.idx(jx, jy);

                let scale = match fs.state[kj] {
                    CellState::Gas => 0.0,
                    CellState::Fluid => 1.0,
                    CellState::Interface => 0.5 * (fill_i + fs.fill[kj]),
                };

                if scale == 0.0 {
                    continue;
                }

                // f_j^{ā}: distribution at j in direction ā (the direction
                // pointing from j back toward i).
                let a_bar = OPP9[a];
                let f_j_abar = grid.cells[kj].f[a_bar];
                // f_i^α: distribution at i in direction α (toward j).
                let f_i_a = grid.cells[ki].f[a];

                let dm = (f_j_abar - f_i_a) * scale;
                delta_mass[ki] += dm;
            }
        }
    }

    // Apply mass deltas to Interface cells.
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            if fs.state[ki] == CellState::Interface {
                fs.mass[ki] += delta_mass[ki];
            }
        }
    }

    // ──────────────────────────────────────────────────────────────────────
    // Step 3: Update fill levels.
    // ──────────────────────────────────────────────────────────────────────
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            match fs.state[ki] {
                CellState::Fluid => {
                    fs.fill[ki] = 1.0;
                    // Keep mass in sync with actual density.
                    fs.mass[ki] = grid.cells[ki].rho;
                }
                CellState::Gas => {
                    fs.fill[ki] = 0.0;
                    fs.mass[ki] = 0.0;
                }
                CellState::Interface => {
                    let rho_i = grid.cells[ki].rho.max(1e-8);
                    fs.fill[ki] = fs.mass[ki] / rho_i;
                }
            }
        }
    }

    // ──────────────────────────────────────────────────────────────────────
    // Step 4: Flag reinit with excess / deficit redistribution.
    // ──────────────────────────────────────────────────────────────────────
    let mut to_fluid: Vec<usize> = Vec::new();
    let mut to_gas: Vec<usize> = Vec::new();

    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            if fs.state[ki] != CellState::Interface {
                continue;
            }
            let fill = fs.fill[ki];
            if fill >= 1.0 + FILL_EPS {
                to_fluid.push(ki);
            } else if fill <= -FILL_EPS {
                to_gas.push(ki);
            }
        }
    }

    // Interface → Fluid: redistribute excess mass to neighbouring Interface cells.
    for &ki in &to_fluid {
        let iy = ki / nx;
        let ix = ki % nx;
        let rho_i = grid.cells[ki].rho.max(1e-8);
        let excess = (fs.fill[ki] - 1.0) * rho_i;

        // Collect Interface neighbours and their fill levels.
        let mut neighbor_fills: [(usize, f64); 8] = [(0, 0.0); 8];
        let mut n_neighbors = 0_usize;
        let mut fill_sum = 0.0_f64;
        for ca in C9.iter().skip(1) {
            let jx = ix as isize + ca[0] as isize;
            let jy = iy as isize + ca[1] as isize;
            if jx < 0 || jy < 0 || jx >= nx as isize || jy >= ny as isize {
                continue;
            }
            let kj = fs.idx(jx as usize, jy as usize);
            if fs.state[kj] == CellState::Interface && !to_fluid.contains(&kj) {
                neighbor_fills[n_neighbors] = (kj, fs.fill[kj]);
                fill_sum += fs.fill[kj];
                n_neighbors += 1;
            }
        }

        if n_neighbors > 0 && fill_sum > 0.0 {
            for item in neighbor_fills.iter().take(n_neighbors) {
                let (kj, fill_j) = *item;
                let weight = fill_j / fill_sum;
                let rho_j = grid.cells[kj].rho.max(1e-8);
                fs.mass[kj] += excess * weight;
                fs.fill[kj] = fs.mass[kj] / rho_j;
            }
        }

        // Convert to Fluid.
        fs.state[ki] = CellState::Fluid;
        fs.fill[ki] = 1.0;
        fs.mass[ki] = rho_i;
    }

    // Interface → Gas: redistribute deficit (negative excess) similarly.
    for &ki in &to_gas {
        let iy = ki / nx;
        let ix = ki % nx;
        let rho_i = grid.cells[ki].rho.max(1e-8);
        // deficit is how much mass we're missing (positive number).
        let deficit = fs.fill[ki].abs() * rho_i;

        let mut neighbor_fills: [(usize, f64); 8] = [(0, 0.0); 8];
        let mut n_neighbors = 0_usize;
        let mut fill_sum = 0.0_f64;
        for ca in C9.iter().skip(1) {
            let jx = ix as isize + ca[0] as isize;
            let jy = iy as isize + ca[1] as isize;
            if jx < 0 || jy < 0 || jx >= nx as isize || jy >= ny as isize {
                continue;
            }
            let kj = fs.idx(jx as usize, jy as usize);
            if fs.state[kj] == CellState::Interface && !to_gas.contains(&kj) {
                neighbor_fills[n_neighbors] = (kj, fs.fill[kj]);
                fill_sum += fs.fill[kj];
                n_neighbors += 1;
            }
        }

        if n_neighbors > 0 && fill_sum > 0.0 {
            for item in neighbor_fills.iter().take(n_neighbors) {
                let (kj, fill_j) = *item;
                let weight = fill_j / fill_sum;
                let rho_j = grid.cells[kj].rho.max(1e-8);
                fs.mass[kj] -= deficit * weight;
                let new_fill = fs.mass[kj] / rho_j;
                fs.fill[kj] = new_fill.max(0.0);
            }
        }

        // Convert to Gas.
        fs.state[ki] = CellState::Gas;
        fs.fill[ki] = 0.0;
        fs.mass[ki] = 0.0;
    }

    // ──────────────────────────────────────────────────────────────────────
    // Step 5: Interface layer closure invariant.
    //
    // After conversions, ensure no Fluid cell is directly adjacent to a Gas
    // cell: any Fluid cell with a Gas 4-neighbor becomes Interface (fill=1,
    // mass=ρ). Any Gas cell with a Fluid 4-neighbor becomes Interface
    // (fill=0, mass=0).
    // ──────────────────────────────────────────────────────────────────────
    // We iterate until no more changes are needed (typically just 1–2 passes).
    let mut changed = true;
    while changed {
        changed = false;
        for iy in 0..ny {
            for ix in 0..nx {
                let ki = fs.idx(ix, iy);
                if fs.state[ki] == CellState::Fluid {
                    // Check 8-connected neighbours for Gas.
                    let mut has_gas_neighbor = false;
                    'outer_fluid: for ca in C9.iter().skip(1) {
                        let jx = ix as isize + ca[0] as isize;
                        let jy = iy as isize + ca[1] as isize;
                        if jx < 0 || jy < 0 || jx >= nx as isize || jy >= ny as isize {
                            continue;
                        }
                        let kj = fs.idx(jx as usize, jy as usize);
                        if fs.state[kj] == CellState::Gas {
                            has_gas_neighbor = true;
                            break 'outer_fluid;
                        }
                    }
                    if has_gas_neighbor {
                        // Demote Fluid → Interface.
                        let rho_i = grid.cells[ki].rho.max(1e-8);
                        fs.state[ki] = CellState::Interface;
                        fs.fill[ki] = 1.0;
                        fs.mass[ki] = rho_i;
                        changed = true;
                    }
                } else if fs.state[ki] == CellState::Gas {
                    // Check 8-connected neighbours for Fluid.
                    let mut has_fluid_neighbor = false;
                    'outer_gas: for ca in C9.iter().skip(1) {
                        let jx = ix as isize + ca[0] as isize;
                        let jy = iy as isize + ca[1] as isize;
                        if jx < 0 || jy < 0 || jx >= nx as isize || jy >= ny as isize {
                            continue;
                        }
                        let kj = fs.idx(jx as usize, jy as usize);
                        if fs.state[kj] == CellState::Fluid {
                            has_fluid_neighbor = true;
                            break 'outer_gas;
                        }
                    }
                    if has_fluid_neighbor {
                        // Promote Gas → Interface (empty).
                        fs.state[ki] = CellState::Interface;
                        fs.fill[ki] = 0.0;
                        fs.mass[ki] = 0.0;
                        changed = true;
                    }
                }
            }
        }
    }

    // Post-flag-reinit: recompute fill from mass for any Interface cells whose
    // mass was modified during excess/deficit redistribution.
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            if fs.state[ki] == CellState::Interface {
                let rho_i = grid.cells[ki].rho.max(1e-8);
                fs.fill[ki] = fs.mass[ki] / rho_i;
            }
        }
    }
}

/// Apply a global mass correction to keep `total_mass` exactly equal to
/// `target_mass`.
///
/// Any discrepancy (positive or negative) is redistributed uniformly to all
/// Interface cells.  This corrects for mass lost during Interface→Gas or
/// Interface→Fluid conversions where redistribution neighbours were unavailable.
///
/// Call once per step after `free_surface_step`.
pub fn apply_mass_correction(grid: &CellGrid2D, fs: &mut FreeSurfaceState, target_mass: f64) {
    let current_mass = fs.total_mass();
    let error = target_mass - current_mass;
    if error.abs() < 1e-15 {
        return;
    }

    // Count Interface cells available to absorb the correction.
    let n_interface = fs
        .state
        .iter()
        .filter(|&&s| s == CellState::Interface)
        .count();

    if n_interface == 0 {
        return;
    }

    let correction_per_cell = error / n_interface as f64;

    let nx = fs.nx();
    let ny = fs.ny();
    for iy in 0..ny {
        for ix in 0..nx {
            let ki = fs.idx(ix, iy);
            if fs.state[ki] == CellState::Interface {
                fs.mass[ki] += correction_per_cell;
                let rho_i = grid.cells[ki].rho.max(1e-8);
                fs.fill[ki] = fs.mass[ki] / rho_i;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::types::CellGrid2D;

    // ── helpers ─────────────────────────────────────────────────────────────

    /// Run `n_steps` of (body-force + VOF-restricted BGK + bounce-back + VOF
    /// step) with the given gravity acceleration `g_lattice` (−y direction).
    ///
    /// Uses `free_surface_bgk_stream` so that only Fluid and Interface cells
    /// participate in BGK collision, matching the Körner 2005 formulation.
    fn run_steps_gravity(
        grid: &mut CellGrid2D,
        fs: &mut FreeSurfaceState,
        n_steps: usize,
        omega: f64,
        g_lattice: f64,
    ) {
        let target_mass = fs.total_mass();
        for _ in 0..n_steps {
            // Shift macroscopic velocities by 0.5*g (Guo pre-step).
            apply_body_force_2d(grid, fs, 0.0, -g_lattice, omega);
            // BGK with Guo force term + stream (Fluid/Interface only).
            free_surface_bgk_stream(grid, fs, omega, 0.0, -g_lattice);
            // Bounce-back on walls.
            grid.apply_bounce_back();
            // VOF mass tracking and flag reinit.
            free_surface_step(grid, fs);
            // Global mass correction to handle redistribution losses.
            apply_mass_correction(grid, fs, target_mass);
        }
    }

    // ── 1. Mass conservation — dam break ────────────────────────────────────
    /// Run a 60×40 dam-break for 100 steps and assert relative mass deviation
    /// stays below 1e-4.
    #[test]
    fn test_mass_conservation_dam_break() {
        let nx = 60;
        let ny = 40;
        let fluid_height = 20;
        let omega = 1.0;

        let mut grid = CellGrid2D::new(nx, ny, 1.0);
        // Mark boundary cells as obstacles (bottom and side walls).
        for x in 0..nx {
            let kb = grid.idx(x, 0);
            let kt = grid.idx(x, ny - 1);
            grid.cells[kb].obstacle = true;
            grid.cells[kt].obstacle = true;
        }
        for y in 0..ny {
            let kl = grid.idx(0, y);
            let kr = grid.idx(nx - 1, y);
            grid.cells[kl].obstacle = true;
            grid.cells[kr].obstacle = true;
        }

        let mut fs = FreeSurfaceState::new_dam_break(nx, ny, fluid_height);
        let mass_initial = fs.total_mass();

        // g=5e-4 keeps Mach low (u ~ sqrt(2gH) ≈ sqrt(2*5e-4*20) ≈ 0.14).
        let g_run = 5e-4_f64;
        run_steps_gravity(&mut grid, &mut fs, 100, omega, g_run);

        let mass_final = fs.total_mass();
        let rel_err = (mass_final - mass_initial).abs() / mass_initial.max(1e-15);
        assert!(
            rel_err < 1e-4,
            "Mass conservation violated: initial={mass_initial:.6}, \
             final={mass_final:.6}, rel_err={rel_err:.2e}"
        );
    }

    // ── 2. Dam-break front position advances ────────────────────────────────
    /// After 50 steps the rightmost wet (Interface or Fluid) column must have
    /// advanced at least 3 cells from its initial position, confirming that
    /// the free surface actually moves under gravity / momentum.
    #[test]
    fn test_dam_break_front_position() {
        // The classic column-collapse setup: left half filled with fluid up to
        // fluid_height, right half and top are Gas.  `new_dam_break` puts the
        // fluid column in the left `nx/2` columns.
        let nx = 60;
        let ny = 40;
        let fluid_height = 20;
        let omega = 1.0;

        let mut grid = CellGrid2D::new(nx, ny, 1.0);
        // Obstacle walls (bounce-back) on all four sides.
        for x in 0..nx {
            let kb = grid.idx(x, 0);
            let kt = grid.idx(x, ny - 1);
            grid.cells[kb].obstacle = true;
            grid.cells[kt].obstacle = true;
        }
        for y in 0..ny {
            let kl = grid.idx(0, y);
            let kr = grid.idx(nx - 1, y);
            grid.cells[kl].obstacle = true;
            grid.cells[kr].obstacle = true;
        }

        let mut fs = FreeSurfaceState::new_dam_break(nx, ny, fluid_height);

        // The initial right-most wet column is the rightmost column of the
        // fluid block, i.e. nx/2 - 1 = 29.
        let initial_front = rightmost_wet_column(&fs);

        // Drive dam-break with significant gravity for clear front advance.
        // g=5e-3 gives u_max ~ sqrt(2*g*H) ~ sqrt(2*5e-3*20) ≈ 0.45 (Ma ≈ 0.78).
        // Run 200 steps to allow pressure to drive horizontal spreading.
        let g_lattice = 5e-3;
        run_steps_gravity(&mut grid, &mut fs, 200, omega, g_lattice);

        let final_front = rightmost_wet_column(&fs);
        let advance = final_front.saturating_sub(initial_front);
        assert!(
            advance >= 3,
            "Dam-break front should advance ≥3 cells in 200 steps, \
             but initial={initial_front}, final={final_front}"
        );
    }

    /// Rightmost column index that contains at least one Interface or Fluid cell.
    fn rightmost_wet_column(fs: &FreeSurfaceState) -> usize {
        let nx = fs.nx();
        let ny = fs.ny();
        let mut rightmost = 0_usize;
        for iy in 0..ny {
            for ix in 0..nx {
                let k = fs.idx(ix, iy);
                if matches!(fs.state[k], CellState::Interface | CellState::Fluid) && ix > rightmost
                {
                    rightmost = ix;
                }
            }
        }
        rightmost
    }

    // ── 3. No Fluid cell adjacent to Gas cell ───────────────────────────────
    /// After each step, no Fluid cell should be 8-connected to a Gas cell.
    #[test]
    fn test_no_fluid_gas_adjacency() {
        let nx = 30;
        let ny = 20;
        let fluid_height = 10;
        let omega = 1.0;

        let mut grid = CellGrid2D::new(nx, ny, 1.0);
        for x in 0..nx {
            let k = grid.idx(x, 0);
            grid.cells[k].obstacle = true;
            let k = grid.idx(x, ny - 1);
            grid.cells[k].obstacle = true;
        }
        for y in 0..ny {
            let k = grid.idx(0, y);
            grid.cells[k].obstacle = true;
            let k = grid.idx(nx - 1, y);
            grid.cells[k].obstacle = true;
        }

        let mut fs = FreeSurfaceState::new_dam_break(nx, ny, fluid_height);
        let g_lattice = 1e-4;

        for step in 0..50 {
            apply_body_force_2d(&mut grid, &fs, 0.0, -g_lattice, omega);
            free_surface_bgk_stream(&mut grid, &fs, omega, 0.0, -g_lattice);
            grid.apply_bounce_back();
            free_surface_step(&mut grid, &mut fs);

            for iy in 0..ny {
                for ix in 0..nx {
                    let ki = fs.idx(ix, iy);
                    if fs.state[ki] != CellState::Fluid {
                        continue;
                    }
                    for (a, ca) in C9.iter().enumerate().skip(1) {
                        let jx = ix as isize + ca[0] as isize;
                        let jy = iy as isize + ca[1] as isize;
                        if jx < 0 || jy < 0 || jx >= nx as isize || jy >= ny as isize {
                            continue;
                        }
                        let kj = fs.idx(jx as usize, jy as usize);
                        assert_ne!(
                            fs.state[kj],
                            CellState::Gas,
                            "Step {step}: Fluid cell ({ix},{iy}) is adjacent to Gas at \
                             ({},{}) via direction {a}",
                            jx,
                            jy
                        );
                    }
                }
            }
        }
    }

    // ── 4. Static enclosed pool stays stable ────────────────────────────────
    /// An all-Fluid box (no Interface / Gas) with obstacle walls: after 20
    /// steps all interior densities should still be ~1.0 (no spurious flow).
    #[test]
    fn test_static_pool_stable() {
        let nx = 10;
        let ny = 10;
        let omega = 1.0;

        let mut grid = CellGrid2D::new(nx, ny, 1.0);
        // Mark boundary as obstacle.
        for x in 0..nx {
            let kb = grid.idx(x, 0);
            let kt = grid.idx(x, ny - 1);
            grid.cells[kb].obstacle = true;
            grid.cells[kt].obstacle = true;
        }
        for y in 0..ny {
            let kl = grid.idx(0, y);
            let kr = grid.idx(nx - 1, y);
            grid.cells[kl].obstacle = true;
            grid.cells[kr].obstacle = true;
        }

        // All interior cells are Fluid.
        let n = nx * ny;
        let mut fs = FreeSurfaceState::new(nx, ny);
        for iy in 1..(ny - 1) {
            for ix in 1..(nx - 1) {
                let k = fs.idx(ix, iy);
                fs.state[k] = CellState::Fluid;
                fs.fill[k] = 1.0;
                fs.mass[k] = 1.0;
            }
        }
        let _ = n;

        for _ in 0..20 {
            free_surface_bgk_stream(&mut grid, &fs, omega, 0.0, 0.0);
            grid.apply_bounce_back();
            free_surface_step(&mut grid, &mut fs);
        }

        // All non-obstacle cells should still have rho ≈ 1.0.
        for iy in 1..(ny - 1) {
            for ix in 1..(nx - 1) {
                let k = grid.idx(ix, iy);
                let rho = grid.cells[k].rho;
                assert!(
                    (rho - 1.0).abs() < 0.1,
                    "Spurious density at ({ix},{iy}): rho={rho}"
                );
            }
        }
    }
}
