// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Phase 21.3 validation: D3Q19 LBM lid-driven cubic cavity at Re=100.
//!
//! Reproduces the Ghia, Ghia & Shin (1982) Table I `u/u_lid` profile along
//! the vertical centerline within a 10% relative tolerance (combined with a
//! small absolute floor to handle the near-zero reference value at
//! y/L ≈ 0.7344 without hypersensitive rejection).
//!
//! The test is expensive: 33³ = 35 937 cells × 19 directions × up to
//! 20 000 time-steps. It is gated behind `#[ignore]` and has to be launched
//! explicitly with
//!
//! ```bash
//! cargo test -p oxiphysics --test validation_lbm -- --ignored
//! ```

#![allow(missing_docs)]

use oxiphysics::lbm::{LbmGrid3D, lattice::LatticeType, stream_3d};

// ---------------------------------------------------------------------------
// Ghia, Ghia & Shin (1982), Table I — Re = 100.
// (y/L, u/u_lid) pairs along the vertical centerline.
// Mirrors `crates/oxiphysics/tests/regression_baselines/lbm_ghia_re100_u_midline.json`.
// ---------------------------------------------------------------------------

const GHIA_RE100: [(f64, f64); 17] = [
    (1.0000, 1.00000),
    (0.9766, 0.84123),
    (0.9688, 0.78871),
    (0.9609, 0.73722),
    (0.9531, 0.68717),
    (0.8516, 0.23151),
    (0.7344, 0.00332),
    (0.6172, -0.13641),
    (0.5000, -0.20581),
    (0.4531, -0.21090),
    (0.2813, -0.15662),
    (0.1719, -0.10150),
    (0.1016, -0.06434),
    (0.0703, -0.04775),
    (0.0625, -0.04192),
    (0.0547, -0.03717),
    (0.0000, 0.00000),
];

// ---------------------------------------------------------------------------
// D3Q19 velocity table (matches `oxiphysics_lbm::lattice::D3Q19_VELOCITIES`)
// Reproduced locally to keep the lattice indexing obvious in this test.
// ---------------------------------------------------------------------------

const D3Q19_C: [[i32; 3]; 19] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, 1, 0],
    [1, -1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, 1],
    [1, 0, -1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, 1],
    [0, 1, -1],
    [0, -1, -1],
];

const D3Q19_W: [f64; 19] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

const D3Q19_OPP: [usize; 19] = [
    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15,
];

// ---------------------------------------------------------------------------
// Equilibrium distribution for D3Q19.
// ---------------------------------------------------------------------------

#[inline]
fn feq_d3q19(w: f64, rho: f64, ux: f64, uy: f64, uz: f64, cx: f64, cy: f64, cz: f64) -> f64 {
    // cs² = 1/3
    let cs2 = 1.0_f64 / 3.0;
    let eu = cx * ux + cy * uy + cz * uz;
    let u_sq = ux * ux + uy * uy + uz * uz;
    w * rho * (1.0 + eu / cs2 + eu * eu / (2.0 * cs2 * cs2) - u_sq / (2.0 * cs2))
}

// ---------------------------------------------------------------------------
// BGK collision on fluid cells only.
// ---------------------------------------------------------------------------

fn collide_fluid_bgk(grid: &mut LbmGrid3D, omega: f64, is_wall: &[bool]) {
    let nx = grid.nx;
    let ny = grid.ny;
    let nz = grid.nz;
    let n = nx * ny * nz;
    for (k, &wall) in is_wall.iter().enumerate().take(n) {
        if wall {
            continue;
        }
        // Read macroscopics (already current).
        let rho = grid.rho[k];
        let ux = grid.ux[k];
        let uy = grid.uy[k];
        let uz = grid.uz[k];
        for i in 0..19 {
            let w = D3Q19_W[i];
            let cx = D3Q19_C[i][0] as f64;
            let cy = D3Q19_C[i][1] as f64;
            let cz = D3Q19_C[i][2] as f64;
            let feq = feq_d3q19(w, rho, ux, uy, uz, cx, cy, cz);
            grid.f[i][k] -= omega * (grid.f[i][k] - feq);
        }
    }
}

// ---------------------------------------------------------------------------
// Full-way bounce-back on a list of wall cells. Swap f[i] ↔ f[opp(i)].
// ---------------------------------------------------------------------------

fn apply_bounce_back(grid: &mut LbmGrid3D, wall_cells: &[usize]) {
    for &k in wall_cells {
        let mut tmp = [0.0_f64; 19];
        for (i, tmp_i) in tmp.iter_mut().enumerate() {
            *tmp_i = grid.f[i][k];
        }
        for (i, f_i) in grid.f.iter_mut().enumerate() {
            f_i[k] = tmp[D3Q19_OPP[i]];
        }
    }
}

// ---------------------------------------------------------------------------
// Moving-wall full-way bounce-back with Ladd momentum correction.
// Correction sign matches `oxiphysics_lbm::streaming::stream_half_way_bb_moving_wall_2d`:
//    f[opp][k] = f_pre[i][k] + 2 w_i ρ (e_i · u_wall) / cs²
// Here we operate in place on the lid nodes: after standard streaming moved
// populations around, swap f[i] ↔ f[opp(i)] at lid cell and add the
// +2 w ρ e·u / cs² correction to the outgoing (post-swap) population.
// ---------------------------------------------------------------------------

fn apply_bounce_back_moving(grid: &mut LbmGrid3D, lid_cells: &[usize], u_wall: [f64; 3]) {
    let cs2 = 1.0_f64 / 3.0;
    for &k in lid_cells {
        let mut tmp = [0.0_f64; 19];
        for (i, tmp_i) in tmp.iter_mut().enumerate() {
            *tmp_i = grid.f[i][k];
        }
        let rho = tmp.iter().sum::<f64>();
        for (i, f_i) in grid.f.iter_mut().enumerate() {
            let opp = D3Q19_OPP[i];
            // e_opp direction carries the reflected population. The Ladd
            // correction is applied to the OUTGOING direction (opp), using
            // e_opp · u_wall.
            let cx = D3Q19_C[opp][0] as f64;
            let cy = D3Q19_C[opp][1] as f64;
            let cz = D3Q19_C[opp][2] as f64;
            let e_dot_u = cx * u_wall[0] + cy * u_wall[1] + cz * u_wall[2];
            let corr = 2.0 * D3Q19_W[opp] * rho * e_dot_u / cs2;
            f_i[k] = tmp[opp] + corr;
        }
    }
}

// ---------------------------------------------------------------------------
// Helper: linear interpolation of u(y) along y = 0..ny-1 at (x_mid, z_mid).
// Returns u/u_lid at an arbitrary real y ∈ [0, ny-1].
// ---------------------------------------------------------------------------

fn interp_u(ux_profile: &[f64], y_real: f64, u_lid: f64) -> f64 {
    if ux_profile.is_empty() || u_lid.abs() < 1e-30 {
        return 0.0;
    }
    let ny = ux_profile.len();
    let y_clamped = y_real.clamp(0.0, (ny - 1) as f64);
    let y_lo = y_clamped.floor() as usize;
    let y_hi = (y_lo + 1).min(ny - 1);
    let t = y_clamped - (y_lo as f64);
    let u = (1.0 - t) * ux_profile[y_lo] + t * ux_profile[y_hi];
    u / u_lid
}

// ---------------------------------------------------------------------------
// The Ghia validation test.
// ---------------------------------------------------------------------------

/// Lid-driven cubic cavity, Re = 100, D3Q19 BGK, 33³ grid.
///
/// The test is marked `#[ignore]` because 33³ × 19 × up to 20 000 steps
/// exceeds the CI budget by a wide margin. Run locally with:
///
/// ```bash
/// cargo test -p oxiphysics --test validation_lbm -- --ignored --nocapture
/// ```
#[test]
#[ignore = "runtime > 60 s: 33^3 × 19 × up to 20 000 LBM steps; run with --ignored"]
fn test_cavity_re100_centerline() {
    // ---- Grid & flow parameters -----------------------------------------
    let nx: usize = 33;
    let ny: usize = 33;
    let nz: usize = 33;
    let l_grid: f64 = (ny - 1) as f64; // 32 grid spans between y=0 and y=ny-1
    let u_lid: f64 = 0.1;
    let reynolds: f64 = 100.0;
    let nu = u_lid * l_grid / reynolds; // 0.1 * 32 / 100 = 0.032
    // BGK: ν = cs² (1/ω − 1/2) = (1/3)(1/ω − 1/2)
    // ⇒ ω = 1 / (3 ν + 0.5)
    let omega = 1.0 / (3.0 * nu + 0.5);
    // Sanity: ν = 0.032 ⇒ ω ≈ 1 / 0.596 ≈ 1.6779
    debug_assert!(
        omega > 1.6 && omega < 1.7,
        "ω out of expected band: {omega}"
    );

    // ---- Build grid & wall masks ----------------------------------------
    let mut grid = LbmGrid3D::new(nx, ny, nz, LatticeType::D3Q19);

    let mut is_wall = vec![false; nx * ny * nz];
    let mut static_walls: Vec<usize> = Vec::new();
    let mut lid_cells: Vec<usize> = Vec::new();

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let k = grid.idx(x, y, z);
                let at_bottom = y == 0;
                let at_top = y == ny - 1;
                let at_xlow = x == 0;
                let at_xhigh = x == nx - 1;
                let at_zlow = z == 0;
                let at_zhigh = z == nz - 1;
                let any_bound = at_bottom || at_top || at_xlow || at_xhigh || at_zlow || at_zhigh;
                if any_bound {
                    is_wall[k] = true;
                    if at_top && !at_xlow && !at_xhigh && !at_zlow && !at_zhigh {
                        // Interior of the top plane is the moving lid.
                        lid_cells.push(k);
                    } else {
                        // Five stationary walls plus lid edges (treat the lid
                        // edges as stationary no-slip to avoid corner
                        // inconsistency with the x and z walls).
                        static_walls.push(k);
                    }
                }
            }
        }
    }

    // ---- Initial condition: quiescent (ρ=1, u=0) ------------------------
    for k in 0..nx * ny * nz {
        grid.rho[k] = 1.0;
        grid.ux[k] = 0.0;
        grid.uy[k] = 0.0;
        grid.uz[k] = 0.0;
        for (i, f_i) in grid.f.iter_mut().enumerate() {
            f_i[k] = D3Q19_W[i]; // equilibrium at ρ=1, u=0
        }
    }

    // ---- Time-stepping loop ---------------------------------------------
    let u_wall = [u_lid, 0.0, 0.0];
    let max_steps = 20_000_usize;
    let conv_tol = 1e-5_f64;
    let conv_start_step = 500_usize;

    let x_mid = nx / 2;
    let z_mid = nz / 2;

    // Previous-centerline snapshot for convergence check.
    let mut prev_ux_line: Vec<f64> = (0..ny)
        .map(|y| grid.ux[grid.idx(x_mid, y, z_mid)])
        .collect();

    let mut converged = false;
    let mut final_step = max_steps;
    let mut final_delta = f64::INFINITY;

    for step in 0..max_steps {
        // (1) collision on fluid cells
        collide_fluid_bgk(&mut grid, omega, &is_wall);
        // (2) streaming (periodic; boundary cells will have their values
        //     fixed by bounce-back immediately after)
        stream_3d(&mut grid);
        // (3) bounce-back: stationary walls, then moving lid
        apply_bounce_back(&mut grid, &static_walls);
        apply_bounce_back_moving(&mut grid, &lid_cells, u_wall);
        // (4) refresh macroscopics
        grid.compute_macroscopic();

        // Convergence check on the centerline profile every 50 steps.
        if step >= conv_start_step && step.is_multiple_of(50) {
            let curr: Vec<f64> = (0..ny)
                .map(|y| grid.ux[grid.idx(x_mid, y, z_mid)])
                .collect();
            let delta = curr
                .iter()
                .zip(prev_ux_line.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max)
                / u_lid;
            final_delta = delta;
            final_step = step + 1;
            if delta < conv_tol {
                converged = true;
                break;
            }
            prev_ux_line = curr;
        }
    }

    if !converged {
        eprintln!(
            "[validation_lbm] WARNING: cavity did not fully converge in {max_steps} steps \
             (max |Δu|/u_lid = {final_delta:.3e} at step {final_step}); \
             proceeding with profile-shape check."
        );
    }

    // ---- Read centerline profile and compare to Ghia --------------------
    let ux_line: Vec<f64> = (0..ny)
        .map(|y| grid.ux[grid.idx(x_mid, y, z_mid)])
        .collect();

    // Combined tolerance: 10% relative OR an absolute floor of 0.01 u_lid
    // (needed because Ghia's value at y/L = 0.7344 is 0.00332 — a strict
    // relative test becomes pathologically sensitive near zero crossings).
    let rel_tol = 0.10_f64;
    let abs_floor = 0.01_f64;

    let mut max_err_rel = 0.0_f64;
    let mut worst_point = (0.0, 0.0, 0.0); // (y/L, measured, expected)
    let mut failed_pairs: Vec<(f64, f64, f64, f64)> = Vec::new();

    for &(y_over_l, ref_u) in GHIA_RE100.iter() {
        // Skip the prescribed / degenerate endpoints.
        if y_over_l <= 0.0 + 1e-12 || y_over_l >= 1.0 - 1e-12 {
            continue;
        }
        let y_real = y_over_l * l_grid;
        let measured = interp_u(&ux_line, y_real, u_lid);
        let diff = (measured - ref_u).abs();
        let tol = rel_tol * ref_u.abs() + abs_floor;
        let rel = if ref_u.abs() > 1e-12 {
            diff / ref_u.abs()
        } else {
            diff
        };
        if rel > max_err_rel {
            max_err_rel = rel;
            worst_point = (y_over_l, measured, ref_u);
        }
        if diff > tol {
            failed_pairs.push((y_over_l, measured, ref_u, diff));
        }
    }

    eprintln!(
        "[validation_lbm] converged = {converged}, steps = {final_step}, \
         final |Δu|/u_lid = {final_delta:.3e}, ω = {omega:.6}, ν = {nu:.6}, \
         worst-point y/L = {:.4}  measured = {:.5}  ref = {:.5}  rel_err = {:.3}",
        worst_point.0, worst_point.1, worst_point.2, max_err_rel,
    );

    assert!(
        failed_pairs.is_empty(),
        "Cavity Re=100 centerline: {} points exceed combined tolerance \
         (10% rel + {abs_floor} abs). Failures (y/L, measured, ref, |Δ|): {:?}. \
         max rel err = {:.3} at y/L = {:.4} (measured = {:.5}, ref = {:.5}).",
        failed_pairs.len(),
        failed_pairs,
        max_err_rel,
        worst_point.0,
        worst_point.1,
        worst_point.2,
    );
}
